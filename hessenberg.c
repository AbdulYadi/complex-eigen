#include <math.h>
#include "complex_matrix.h"

static void housegen(const matrix *m, unsigned int startrow, unsigned int stoprow, unsigned int col,
    vector *v_u, complex *v);

void hessenberg(const matrix *m, matrix *h, matrix *q)
{
////assumption: m is squared matrix
    int row, col, col2;
    unsigned int input_row, input_col;
    vector vcol_u, vcol_v, vrow_u, vrow_v;
    matrix vcol_u_vrow_v, vcol_v_vrow_u;

    matrix_dup(m, h);
    matrix_create(m->rows, m->cols, q);

    vector_create(m->cols-1/*maximum count*/, &vcol_u);
    vector_create(m->rows, &vcol_v);
    vector_create(m->cols-1/*maximum count*/, &vrow_u);
    vector_create(m->cols-1/*maximum count*/, &vrow_v);
    matrix_create(m->rows-1, m->cols-1, &vcol_u_vrow_v);
    matrix_create(m->rows, m->cols-1, &vcol_v_vrow_u);

////reduce the matrix
    for(col = 0; col < h->cols-2; col++)
    {
        vrow_v.count = 1 + (m->rows-1) - (col+1);

    ////generate the transformation
        housegen(h, col+1, h->rows-1, col, &vcol_u, &MATRIXP(h, col+1, col));
        matrix_apply_vector(&vcol_u, col+1, col, q);

    ////premultiply the transformation
        vector_conjugate(&vcol_u, &vrow_u);
        rowvector_mul_submatrix(&vrow_u, h, col+1, col+1, &vrow_v);
        vector_product_to_matrix(&vcol_u, &vrow_v, &vcol_u_vrow_v, false);
        for(row=col+1, input_row=0; row<h->rows; row++, input_row++)
        {
            for(col2=col+1, input_col=0; col2<h->cols; col2++, input_col++)
                complex_self_sub(&MATRIXP(h, row, col2), &MATRIX(vcol_u_vrow_v, input_row, input_col));
            if(row>col+1)
                complex_init(&MATRIXP(h, row, col));
        }

    ////postmultiply the transformation
        submatrix_mul_colvector(h, &vcol_u, 0, col+1, &vcol_v);
        vector_product_to_matrix(&vcol_v, &vrow_u, &vcol_v_vrow_u, false);
        for(row=0; row<h->rows; row++)
        {
            for(col2=col+1, input_col=0; col2<h->cols; col2++, input_col++)
                complex_self_sub(&MATRIXP(h, row, col2), &MATRIX(vcol_v_vrow_u, row, input_col));
        }
    }

////accumulate the transformations
    ////unit vector
    matrix_unit_vector(q, 2, q->cols-1, q->cols-2);
    for(col = q->cols-3; col>=0; col--)
    {
        vrow_v.count = vcol_u.count = 1 + (q->rows-1) - (col+1);
        matrix_to_colvector(q, col+1, q->rows-1, col, &vcol_u);
        vector_conjugate(&vcol_u, &vrow_u);
        rowvector_mul_submatrix(&vrow_u, q, col+1, col+1, &vrow_v);
        vector_product_to_matrix(&vcol_u, &vrow_v, &vcol_u_vrow_v, false);
        for(row=col+1, input_row=0; row<q->rows; row++, input_row++)
        {
            for(col2=col+1, input_col=0; col2<q->cols; col2++, input_col++)
                complex_self_sub(&MATRIXP(q, row, col2), &MATRIX(vcol_u_vrow_v, input_row, input_col));
        }
        ////unit vector
        matrix_unit_vector(q, 1, col);
    }//for(col = q->cols-3; col>=0; col--)

    vector_destroy(&vcol_u);
    vector_destroy(&vcol_v);
    vector_destroy(&vrow_u);
    vector_destroy(&vrow_v);
    matrix_destroy(&vcol_u_vrow_v);
    matrix_destroy(&vcol_v_vrow_u);
}

static void housegen(const matrix *m, unsigned int startrow, unsigned int stoprow, unsigned int col,
    vector *v_u, complex *v)
{
    complex p, pconjugate, *u0, p_over_v, tmp, u0sqrt;
    float num;

    v_u->count = 1 + stoprow - startrow;
    matrix_to_colvector(m, startrow, stoprow, col, v_u);

    complex_init(v);
    v->r = vector_norm2(v_u);
    u0 = &VECTORP(v_u, 0);

    if(v->r == 0)
    {
        u0->r = sqrt(2);
        return;
    }

    complex_init(&p);
    if( COMPLEXP_IS_ZERO(u0) )
        p.r = 1;
    else
    {
        num = complex_abs(u0);
        p.r = u0->r / num;
        p.i = -(u0->i) / num;
    }

    complex_div(&p, v, &p_over_v);
    vector_scalar_mul(v_u, &p_over_v);

    u0->r += 1;
    complex_sqrt(u0, &u0sqrt);
    vector_scalar_div(v_u, &u0sqrt);
    complex_conjugate(&p, &pconjugate);
    complex_mul(v, &pconjugate, &tmp);
    v->r = -(tmp.r);
    v->i = -(tmp.i);
}
