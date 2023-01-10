#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/param.h>
#include <stdarg.h>
#include "complex_matrix.h"

void complex_init(complex *c)
{
    c->r = c->i = 0;
}

void complex_copy(const complex *src, complex *dest)
{
    dest->r = src->r;
    dest->i = src->i;
}

float complex_abs(const complex *c)
{
    return sqrt((c->r * c->r) + (c->i * c->i));
}

void complex_div(const complex *num, const complex *denom, complex *dest)
{
    float d;

    d = (denom->r * denom->r) + (denom->i * denom->i);
    dest->r = ((num->r * denom->r) + (num->i * denom->i)) / d;
    dest->i = ((num->i * denom->r) - (num->r * denom->i)) / d;
}

void complex_div_real(const complex *num, float denom, complex *dest)
{
    dest->r = num->r / denom;
    dest->i = num->i / denom;
}

void complex_mul(const complex *c1, const complex *c2, complex *dest)
{
    dest->r = (c1->r * c2->r) - (c1->i * c2->i);
    dest->i = (c1->r * c2->i) + (c1->i * c2->r);
}

void complex_mul_real(const complex *num, float mul, complex *dest)
{
    dest->r = num->r * mul;
    dest->i = num->i * mul;
}

void complex_add(const complex *c, const complex *add, complex *dest)
{
    dest->r = c->r + add->r;
    dest->i = c->i + add->i;
}

void complex_sub(const complex *c, const complex *sub, complex *dest)
{
    dest->r = c->r - sub->r;
    dest->i = c->i - sub->i;
}

void complex_sqrt(const complex *c, complex *dest)
{
    float z;

    if(c->i==0)
    {
        if(c->r<0)
        {
            dest->r = 0;
            dest->i = sqrt(c->r * -1);
        }
        else
        {
            dest->r = sqrt(c->r);
            dest->i = 0;
        }
    }
    else
    {
        z = sqrt((c->r * c->r) + (c->i * c->i));
        dest->r = sqrt((z + c->r)/2);
        dest->i = (c->i/fabsf(c->i)) * sqrt((z - c->r)/2);
    }
}

void complex_conjugate(const complex *c, complex *dest)
{
    dest->r = c->r;
    dest->i = -(c->i);
}

void complex_invert(const complex *c, complex *dest)
{
    float denom;
    denom = (c->r * c->r) +  (c->i * c->i);
    dest->r = c->r / denom;
    dest->i = -(c->i / denom);
}

void complex_self_add(complex *self, const complex *add)
{
    self->r += add->r;
    self->i += add->i;
}

void complex_self_sub(complex *self, const complex *sub)
{
    self->r -= sub->r;
    self->i -= sub->i;
}

void complex_self_div(complex *self, const complex *denom)
{
    complex tmp;
    complex_div(self, denom, &tmp);
    complex_copy(&tmp, self);
}

void complex_self_div_real(complex *self, float real)
{
    self->r /= real;
    self->i /= real;
}

void complex_self_mul_real(complex *self, float real)
{
    self->r *= real;
    self->i *= real;
}

void complex_self_sqrt(complex *c)
{
    complex dest;
    complex_sqrt(c, &dest);
    complex_copy(&dest, c);
}

void matrix_create(unsigned int rows, unsigned int cols, matrix* m)
{
    m->rows = rows;
    m->cols = cols;
    m->data = (complex*)calloc(rows * cols, sizeof(complex));
}

bool matrix_is_squared(const matrix* m)
{
    return m->rows>0 && m->rows==m->cols;
}

void matrix_init(matrix* m)
{
    m->rows = m->cols = 0;
    m->data = NULL;
}

void matrix_destroy(matrix* m)
{
    if(m->data!=NULL) {
        free(m->data);
        m->data = NULL;
    }
    return;
}

void matrix_fill(matrix* m, ...)
{
    va_list ap;
    int row, col;

    va_start(ap, m);
    for(row=0; row<m->rows; row++)
    {
        for(col = 0; col<m->cols; col++)
        {
            MATRIXP(m, row, col).r = (float)va_arg(ap, double);
            MATRIXP(m, row, col).i = 0;
        }
    }
    va_end(ap);
}

void matrix_fill_int(matrix* m, ...)
{
    va_list ap;
    int row, col;

    va_start(ap, m);
    for(row=0; row<m->rows; row++)
    {
        for(col = 0; col<m->cols; col++)
        {
            MATRIXP(m, row, col).r = (float)va_arg(ap, int);
            MATRIXP(m, row, col).i = 0;
        }
    }
    va_end(ap);
}

void matrix_copy(const matrix* src, matrix *dest)
{
    memcpy(dest->data, src->data, src->rows * src->cols * sizeof(complex));
}

void matrix_dup(const matrix* src, matrix *dest)
{
    dest->rows = src->rows;
    dest->cols = src->cols;
    dest->data = (complex*)malloc(dest->rows * dest->cols * sizeof(complex));
    matrix_copy(src, dest);
}

//matrix_to_vector: dest->count dictates the number to fetch from matrix
void matrix_to_colvector(const matrix* src, unsigned int startrow, unsigned int stoprow, unsigned int col, vector *dest)
{
    int i;
    for(i=0; startrow<src->rows && startrow<=stoprow && i<dest->count; i++, startrow++)
        complex_copy(&MATRIXP(src, startrow, col), &VECTORP(dest, i));
}

//matrix_to_vector: dest->count dictates the row count to fetch from matrix
void matrix_to_rowvector(const matrix* src, unsigned int startcol, unsigned int stopcol, unsigned int row, vector *dest)
{
    int i;
    for(i=0; startcol<src->cols && startcol<=stopcol && i<dest->count; i++)
    {
        complex_copy(&MATRIXP(src, row, startcol), &VECTORP(dest, i));
        startcol++;
    }
}

//matrix_apply_vector: src->count dictates the row count to apply to matrix
void matrix_apply_vector(const vector* src, unsigned int startrow, unsigned int col, matrix *dest)
{
    int i;
    for(i=0; startrow<dest->rows && i<src->count; i++, startrow++)
        complex_copy(&VECTORP(src, i), &MATRIXP(dest, startrow, col));
}

void matrix_mul(const matrix* m1, const matrix* m2, matrix* dest, bool create_buffer)
{
    int row, col, col2, innercols;
    complex c;

    if(create_buffer)
        matrix_create(m1->rows, m2->cols, dest);
    else {
        if( (m1->rows * m2->cols) > (dest->rows * dest->cols) )
            dest->data = (complex*)realloc(dest->data, m1->rows * m1->cols * sizeof(complex));
        dest->rows = m1->rows;
        dest->cols = m2->cols;
    }

    innercols = MIN(m1->cols, m2->rows);
    for(row=0; row<m1->rows; row++)
    {
        for(col=0; col<m2->cols; col++)
        {
            complex_init(&MATRIXP(dest, row, col));
            for(col2=0; col2<innercols; col2++) {
                complex_mul(&MATRIXP(m1, row, col2), &MATRIXP(m2, col2, col), &c);
                complex_self_add(&MATRIXP(dest, row, col), &c);
            }
        }
    }
}

void vector_product_to_matrix(const vector* vcol, const vector* vrow, matrix *dest, bool create_buffer)
{
    int rows, cols, row, col;

    rows = vcol->count;
    cols = vrow->count;
    if(create_buffer)
        matrix_create(rows, cols, dest);
    else {
        if( (rows * cols) > (dest->rows * dest->cols) )
            dest->data = (complex*)realloc(dest->data, rows * cols * sizeof(complex));
        dest->rows = rows;
        dest->cols = cols;
    }

    for(row=0; row<rows; row++)
    {
        for(col=0; col<cols; col++)
            complex_mul(&VECTORP(vcol, row), &VECTORP(vrow, col), &MATRIXP(dest, row, col));
    }
}

void matrix_unit_vector(matrix *m, int colcount, ...)
{
    va_list ap;
    int row, col, i;
    int colarr[5];
    int *cols;

    cols = colcount>5 ? (int*)malloc(colcount * sizeof(int)) : NULL;
    va_start(ap, colcount);
    for(row=0; row<m->rows; row++)
    {
        for(i = 0; i < colcount; i++)
        {
            if(row==0)
            {
                if(cols!=NULL)
                    cols[i] = va_arg(ap, int);
                else
                    colarr[i] = va_arg(ap, int);
            }
            col = cols!=NULL ? cols[i] : colarr[i];
            MATRIXP(m, row, col).r = row==col ? 1 : 0;
            MATRIXP(m, row, col).i = 0;
        }
    }
    va_end(ap);
    if(cols!=NULL)
        free(cols);
}

void rowvector_mul_submatrix(const vector* vrow, const matrix* m, unsigned int matrix_startcol, unsigned int matrix_startrow, vector *vrow_dest)
{
    unsigned int stopcol, stoprow, submatrixcol, submatrixrow;
    unsigned int output_index, input_index;
    complex *cdest, c;

    stopcol = matrix_startcol + vrow->count;
    stoprow = matrix_startrow + vrow->count;
    for(submatrixcol=matrix_startcol, output_index=0; submatrixcol<stopcol; submatrixcol++, output_index++)
    {
        cdest = &VECTORP(vrow_dest, output_index);
        complex_init(cdest);
        for(submatrixrow=matrix_startrow, input_index=0; submatrixrow<stoprow; submatrixrow++, input_index++)
        {
            complex_mul(&VECTORP(vrow, input_index), &MATRIXP(m, submatrixrow, submatrixcol), &c);
            complex_self_add(cdest, &c);
        }
    }
}

void submatrix_mul_colvector(const matrix* m, const vector* vcol, unsigned int matrix_startrow, unsigned int matrix_startcol, vector *vcol_dest)
{
    unsigned int submatrixrow, submatrixcol, stoprow, stopcol;
    unsigned int output_index, input_index;
    complex *cdest, c;

    stoprow = matrix_startrow + vcol_dest->count;
    stopcol = matrix_startcol + vcol->count;
    for(submatrixrow=matrix_startrow, output_index=0; submatrixrow<stoprow; submatrixrow++, output_index++)
    {
        cdest = &VECTORP(vcol_dest, output_index);
        complex_init(cdest);
        for(submatrixcol=matrix_startcol, input_index=0; submatrixcol<stopcol; submatrixcol++, input_index++)
        {
            complex_mul(&MATRIXP(m, submatrixrow, submatrixcol), &VECTORP(vcol, input_index), &c);
            complex_self_add(cdest, &c);
        }
    }
}

float matrix_norm2_col(const matrix *m, unsigned int startrow, unsigned int stoprow, unsigned int col)
{
    unsigned int row;
    float norm = 0.0;

    for(row=startrow; row<=stoprow; row++)
    {
        complex* c = &MATRIXP(m, row, col);
        norm += (c->r * c->r) + (c->i * c->i);
    }
    return sqrt(norm);
}

void vector_create(unsigned int count, vector* v)
{
    v->count = count;
    v->data = (complex*)calloc(count, sizeof(complex));
}

void vector_init(vector* v)
{
    v->count = 0;
    v->data = NULL;
}

void vector_destroy(vector* v)
{
    if(v->data!=NULL) {
        free(v->data);
        v->data = NULL;
    }
    return;
}

float vector_norm2(const vector *v)
{
    int i;
    float norm = 0.0;
    for(i=0; i<v->count; i++) {
        complex* c = &v->data[i];
        norm += (c->r * c->r) + (c->i * c->i);
    }
    return sqrt(norm);
}

void vector_copy(const vector* src, vector *dest)
{
    memcpy(dest->data, src->data, src->count * sizeof(complex));
}

void vector_dup(const vector* src, vector *dest)
{
    vector_create(src->count, dest);
    vector_copy(src, dest);
}

void vector_scalar_mul(vector* v, const complex *c)
{
    complex tmp;
    int i;

    for(i=0; i<v->count; i++) {
        complex* vc = &v->data[i];
        complex_mul(c, vc, &tmp);
        complex_copy(&tmp, vc);
    }
}

void vector_scalar_div(vector* v, const complex *denom)
{
    complex tmp;
    int i;

    for(i=0; i<v->count; i++) {
        complex* vc = &v->data[i];
        complex_div(vc, denom, &tmp);
        complex_copy(&tmp, vc);
    }
}

void vector_conjugate(const vector* v, vector* dest)
{
    int i;

    dest->count = v->count;
    for(i=0; i<v->count; i++)
        complex_conjugate(&VECTORP(v,i), &VECTORP(dest,i));
}
