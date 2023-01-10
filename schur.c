#include <math.h>
#include <values.h>
#include "schur.h"

static void backsearch(matrix *h, unsigned int l, unsigned int* i1, unsigned int* i2);
static void wilkshift(const complex *a, const complex *b, const complex *c, const complex *d, complex *k);
static void rotgen(complex *a, complex *b, complex *cos, complex *sin);
static void rotapp_col(const complex *cos, const complex *sin,
        matrix *m,
        unsigned int xrow, unsigned int yrow,
        unsigned int startcol, unsigned int stopcol);
static void rotapp_row(const complex *cos, const complex *sin,
        matrix *m,
        unsigned int xcol, unsigned int ycol,
        unsigned int startrow, unsigned int stoprow);

unsigned int schur(matrix *h, matrix *q, int maxiter)
{
    unsigned int i, i1, i2, oldi2, iter;
    complex k, sin_conjugate;
    vector cos, sin;

    vector_create(h->cols - 1, &cos);
    vector_create(h->cols - 1, &sin);

    i1 = 0;
    i2 = h->rows - 1;
    iter = 0;

    while(true)
    {
        iter++;
        if(iter>maxiter)
            return SCHUR_CODE_TIMEOUT;

        oldi2 = i2;
        backsearch(h, i2, &i1, &i2);
        if(i2==1)
            break;
        if(i2!=oldi2)
            iter=0;

        wilkshift(
            &MATRIXP(h, i2-1, i2-1),
            &MATRIXP(h, i2-1, i2),
            &MATRIXP(h, i2, i2-1),
            &MATRIXP(h, i2, i2),
            &k);

        complex_self_sub(&MATRIXP(h, i1, i1), &k);

        for(i=i1; i<=i2-1; i++)
        {
            rotgen(
                &MATRIXP(h, i, i),
                &MATRIXP(h, i+1, i),
                &VECTOR(cos, i),
                &VECTOR(sin, i)
            );
            complex_self_sub(&MATRIXP(h, i+1, i+1), &k);
            rotapp_col(&VECTOR(cos, i), &VECTOR(sin, i), h,
                i/*xrow*/, i+1/*yrow*/,
                i+1/*startcol*/, h->cols-1/*stopcol*/);
        }

        for(i=i1; i<=i2-1; i++)
        {
            complex_conjugate(&VECTOR(sin, i), &sin_conjugate);
            rotapp_row(&VECTOR(cos, i), &sin_conjugate, h,
                i/*xcol*/, i+1/*ycol*/,
                0/*startrow*/, i+1/*stoprow*/);
            rotapp_row(&VECTOR(cos, i), &sin_conjugate, q,
                i/*xcol*/, i+1/*ycol*/,
                0/*startrow*/, q->rows-1/*stoprow*/);
            complex_self_add(&MATRIXP(h, i, i), &k);
        }

        complex_self_add(&MATRIXP(h, i2, i2), &k);
    }//while(true)

    vector_destroy(&cos);
    vector_destroy(&sin);

    return SCHUR_CODE_OK;
}

const char* schur_err(unsigned int err)
{
    switch(err)
    {
        case SCHUR_CODE_OK:
            return "Successful";
        case SCHUR_CODE_TIMEOUT:
            return "Schur is timed out";
    }
    return "Unknown";
}

static void backsearch(matrix *h, unsigned int l, unsigned int* i1, unsigned int* i2)
{
////finding deflation rows in an upper Hessenberg matrix
    complex *c;
    float negligible;

    negligible = FLT_MIN;//do not change this!!!
    *i1 = *i2 = l;
    while(*i1 > 0)
    {
        c = &MATRIXP(h, *i1, *(i1)-1);
        if( complex_abs(c) <= negligible )
        {
            complex_init(c);
            if(*i1==*i2)
                *i2 = *i1 = (*i1)-1;
            else
                break;
        }
        else {
            (*i1)--;
        }
    }
}

static void wilkshift(const complex *a, const complex *b, const complex *c, const complex *d, complex *k)
{
////B = {  a b
////       c d
/////   }
////wilkshift returns the eigenvalue of B that is nearest d
    float s;
    complex b_over_s, c_over_s, a_over_s, d_over_s;
    complex q, p, p_squared, r, tmp;
    complex p_plus_r, q_over_p_plus_r;

    complex_copy(d, k);
    s = complex_abs(a) + complex_abs(b) + complex_abs(c) + complex_abs(d);
    if(s!=0)
    {
        complex_div_real(b, s, &b_over_s);
        complex_div_real(c, s, &c_over_s);
        complex_mul(&b_over_s, &c_over_s, &q);
        if(!COMPLEX_IS_ZERO(q))
        {
            complex_div_real(a, s, &a_over_s);
            complex_div_real(d, s, &d_over_s);
            complex_sub(&a_over_s, &d_over_s, &tmp);
            complex_mul_real(&tmp, 0.5, &p);

            complex_mul(&p, &p, &p_squared);
            complex_add(&p_squared, &q, &tmp);
            complex_sqrt(&tmp, &r);

            if( (p.r * r.r) + (p.i * r.i) < 0 )
                complex_self_mul_real(&r, -1);

            complex_add(&p, &r, &p_plus_r);
            complex_div(&q, &p_plus_r, &q_over_p_plus_r);
            complex_mul_real(&q_over_p_plus_r, s, &tmp);
            complex_self_sub(k, &tmp);
        }
    }
    float negligible = powf(2, -24);
    if( fabsf(k->r)<=negligible ){ k->r=0; }
    if( fabsf(k->i)<=negligible ){ k->i=0; }
}

void rotgen(complex *a, complex *b, complex *cos, complex *sin)
{
    float a_abs, b_abs, tau, v, a_over_tau_abs, b_over_tau_abs;
    complex mu, a_over_tau, b_over_tau, tmp;

    if(COMPLEXP_IS_ZERO(b))
    {
        cos->r = 1;
        cos->i = 0;
        complex_init(sin);
        return;
    }

    if(COMPLEXP_IS_ZERO(a))
    {
        complex_init(cos);
        sin->r = 1;
        sin->i = 0;
        complex_copy(b, a);
        complex_init(b);
        return;
    }

    a_abs = complex_abs(a);
    complex_div_real(a, a_abs, &mu);

    b_abs = complex_abs(b);
    tau = a_abs + b_abs;

    complex_div_real(a, tau, &a_over_tau);
    a_over_tau_abs = complex_abs(&a_over_tau);

    complex_div_real(b, tau, &b_over_tau);
    b_over_tau_abs = complex_abs(&b_over_tau);

    v = tau * sqrt( (a_over_tau_abs * a_over_tau_abs) + (b_over_tau_abs * b_over_tau_abs)) ;

    cos->r = a_abs / v;
    cos->i = 0;

    complex_conjugate(b, &tmp);
    complex_self_div_real(&tmp, v);
    complex_mul(&mu, &tmp, sin);

    complex_mul_real(&mu, v, a);
    complex_init(b);
}

static void rotapp_col(const complex *cos, const complex *sin,
        matrix *m,
        unsigned int xrow, unsigned int yrow,
        unsigned int startcol, unsigned int stopcol)
{
    unsigned int col;
    complex a, b, t, sin_conjugate;

    complex_conjugate(sin, &sin_conjugate);
    for(col=startcol; col<=stopcol; col++)
    {
        complex_mul(cos, &MATRIXP(m, xrow, col), &a);
        complex_mul(sin, &MATRIXP(m, yrow, col), &b);
        complex_add(&a, &b, &t);

        complex_mul(cos, &MATRIXP(m, yrow, col), &a);
        complex_mul(&sin_conjugate, &MATRIXP(m, xrow, col), &b);
        complex_sub(&a, &b, &MATRIXP(m, yrow, col));

        complex_copy(&t, &MATRIXP(m, xrow, col));
    }
}

static void rotapp_row(const complex *cos, const complex *sin,
        matrix *m,
        unsigned int xcol, unsigned int ycol,
        unsigned int startrow, unsigned int stoprow)
{
    unsigned int row;
    complex a, b, t, sin_conjugate;

    complex_conjugate(sin, &sin_conjugate);
    for(row=startrow; row<=stoprow; row++)
    {
        complex_mul(cos, &MATRIXP(m, row, xcol), &a);
        complex_mul(sin, &MATRIXP(m, row, ycol), &b);
        complex_add(&a, &b, &t);

        complex_mul(cos, &MATRIXP(m, row, ycol), &a);
        complex_mul(&sin_conjugate, &MATRIXP(m, row, xcol), &b);
        complex_sub(&a, &b, &MATRIXP(m, row, ycol));

        complex_copy(&t, &MATRIXP(m, row, xcol));
    }
}
