#include <math.h>
#include <values.h>
#include <sys/param.h>
#include "hessenberg.h"
#include "schur.h"
#include "eigen.h"

static void rightvec(matrix *t, matrix *x);

unsigned int eigen(const matrix *m, vector *eigenval, matrix *eigenvec)
{
    matrix h, q, x;
    unsigned int err, i;
    complex *tmp, denom;
    bool hit;

    if(!matrix_is_squared(m))
        return EIGEN_CODE_NON_SQUARED_INPUT;

    matrix_init(&h);
    matrix_init(&q);
    matrix_init(&x);

    hessenberg(m, &h, &q);
    err = schur(&h, &q, 100);
    if(err==SCHUR_CODE_OK)
    {
        rightvec(&h, &x);
        matrix_mul(&q, &x, eigenvec, true);

        vector_create(m->cols, eigenval);
        for(i=0; i<m->cols; i++)//get diagonal value for eigenvalue
            complex_copy(&MATRIX(h, i, i), &VECTORP(eigenval, i));

        for(int c=0; c<eigenvec->cols; c++)
        {
            hit = false;
            for(int r=eigenvec->rows-1; !hit && r>=0; r--)
            {
                tmp = &MATRIXP(eigenvec, r, c);
                if(tmp->r!=0 || tmp->i!=0)
                    hit = true;
            }
            if(hit && !(tmp->r==1 && tmp->i==0))
            {
                complex_copy(tmp, &denom);
                for(int r=0; r<eigenvec->rows; r++)
                    complex_self_div(&MATRIXP(eigenvec, r, c), &denom);
            }
        }

        err = EIGEN_CODE_OK;
    }

    matrix_destroy(&h);
    matrix_destroy(&q);
    matrix_destroy(&x);

    return err;
}

const char* eigen_err(unsigned int code)
{
    if(code>=SCHUR_CODE_MIN && code<=SCHUR_CODE_MAX)
        return schur_err(code);
    switch(code)
    {
        case EIGEN_CODE_OK:
            return "Successful";
        case EIGEN_CODE_NON_SQUARED_INPUT:
            return "Input matrix is not squared";
    }
    return "Unknown";
}

bool eigen_success(unsigned int code)
{
    return code==EIGEN_CODE_OK || code==SCHUR_CODE_OK;
}

static void rightvec(matrix *t, matrix *x)
{
    float roundunit, smallnum, bignum, dmin, t_abs, d_abs, x_abs, s, norm2;
    int col, row, col2;
    complex *c, d, tmp;

    roundunit = powf(2, -24);
    smallnum = FLT_MIN;
    bignum = FLT_MAX;

    matrix_create(t->rows, t->cols, x);
    for(col=t->cols-1; col>=0; col--)
    {
        for(row=0; row<col; row++)
            complex_mul_real(&MATRIXP(t, row, col), -1, &MATRIXP(x, row, col));

        c = &MATRIXP(x, col, col);
        c->r = 1;
        c->i = 0;

        for(row=col+1; row<x->rows; row++)
            complex_init(&MATRIXP(x, row, col));

        t_abs = complex_abs(&MATRIXP(t, col, col));
        dmin = MAX(roundunit * t_abs, smallnum);

        for(col2 = col-1; col2>=0; col2--)
        {
            complex_sub(&MATRIXP(t, col2, col2), &MATRIXP(t, col, col), &d);

            d_abs = complex_abs(&d);
            if(d_abs <= dmin)
            {
                d.r = dmin;
                d.i = 0;
                d_abs = complex_abs(&d);
            }

            x_abs = complex_abs(&MATRIXP(x, col2, col));
            if(x_abs/bignum >= d_abs)
            {
                s = d_abs/x_abs;
                for(row=0; row<=col; row++)
                    complex_self_mul_real(&MATRIXP(x, row, col), s);
            }

            complex_self_div(&MATRIXP(x, col2, col), &d);

            for(row=0; row<col2; row++)
            {
                complex_mul(&MATRIXP(x, col2, col), &MATRIXP(t, row, col2), &tmp);
                complex_self_sub(&MATRIXP(x, row, col), &tmp);
            }
        }//for(col2 = col-1; col2>=0; col2--)

        norm2 = matrix_norm2_col(x, 0/*startrow*/, col/*stoprow*/, col/*the column*/);
        for(row=0; row<=col; row++)
            complex_self_div_real(&MATRIXP(x, row, col), norm2);
    }//for(col=t->cols-1; col>=0; col--)
}
