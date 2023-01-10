#ifndef _COMPLEX_MATRIX_
#define _COMPLEX_MATRIX_

#include <stdbool.h>

typedef struct complex {
    float r;//real
    float i;//imaginary
} complex;

#define COMPLEX_IS_ZERO(C) (C.r==0 && C.i==0)
#define COMPLEXP_IS_ZERO(PC) (PC->r==0 && PC->i==0)
extern void complex_init(complex *c);
extern void complex_copy(const complex *src, complex *dest);
extern float complex_abs(const complex *c);
extern void complex_div(const complex *num, const complex *denom, complex *dest);
extern void complex_div_real(const complex *num, float denom, complex *dest);
extern void complex_mul(const complex *c1, const complex *c2, complex *dest);
extern void complex_mul_real(const complex *num, float mul, complex *dest);
extern void complex_add(const complex *c, const complex *add, complex *dest);
extern void complex_sub(const complex *c, const complex *sub, complex *dest);
extern void complex_sqrt(const complex *c, complex *dest);
extern void complex_conjugate(const complex *c, complex *dest);
extern void complex_invert(const complex *c, complex *dest);
extern void complex_self_add(complex *self, const complex *add);
extern void complex_self_sub(complex *self, const complex *sub);
extern void complex_self_div(complex *self, const complex *denom);
extern void complex_self_div_real(complex *self, float real);
extern void complex_self_mul_real(complex *self, float real);
extern void complex_self_sqrt(complex *c);

typedef struct matrix {
    unsigned int rows;
    unsigned int cols;
    complex *data;
} matrix;

#define MATRIX(M,ROW,COL) M.data[(ROW) * M.cols + (COL)]
#define MATRIXP(PM,ROW,COL) (PM)->data[(ROW) * (PM)->cols + (COL)]

typedef struct vector {
    unsigned int count;
    complex *data;
} vector;

#define VECTOR(V,IDX) V.data[IDX]
#define VECTORP(PV,IDX) (PV)->data[IDX]

extern void matrix_create(unsigned int rows, unsigned int cols, matrix* m);
extern void matrix_init(matrix* m);
extern void matrix_destroy(matrix* m);
extern void matrix_fill(matrix* m, ...);
extern void matrix_fill_int(matrix* m, ...);
extern void matrix_copy(const matrix* src, matrix *dest);
extern void matrix_dup(const matrix* src, matrix *dest);
//matrix_to_vector: dest->count dictates the row count to fetch from matrix
extern void matrix_to_colvector(const matrix* src, unsigned int startrow, unsigned int stoprow, unsigned int col, vector *dest);
//matrix_to_vector: dest->count dictates the row count to fetch from matrix
extern void matrix_to_rowvector(const matrix* src, unsigned int startcol, unsigned int stopcol, unsigned int row, vector *dest);
//matrix_apply_vector: src->count dictates the row count to apply to matrix
extern void matrix_apply_vector(const vector* src, unsigned int startrow, unsigned int col, matrix *dest);
extern void matrix_mul(const matrix* m1, const matrix* m2, matrix* dest, bool create_buffer);
extern void vector_product_to_matrix(const vector* vcol, const vector* vrow, matrix *dest, bool create_buffer);
extern void matrix_unit_vector(matrix *m, int colcount, ...);
extern void rowvector_mul_submatrix(const vector* vrow, const matrix* m, unsigned int matrix_startcol, unsigned int matrix_startrow, vector *vrow_dest);
extern void submatrix_mul_colvector(const matrix* m, const vector* vcol, unsigned int matrix_startrow, unsigned int matrix_startcol, vector *vcol_dest);
extern float matrix_norm2_col(const matrix *m, unsigned int startrow, unsigned int stoprow, unsigned int col);

extern void vector_create(unsigned int count, vector* v);
extern void vector_init(vector* v);
extern void vector_destroy(vector* v);
extern float vector_norm2(const vector *v);
extern void vector_copy(const vector* src, vector *dest);
extern void vector_dup(const vector* src, vector *dest);
extern void vector_scalar_mul(vector* v, const complex *c);
extern void vector_scalar_div(vector* v, const complex *denom);
extern void vector_conjugate(const vector* v, vector* dest);

#endif //_COMPLEX_MATRIX_
