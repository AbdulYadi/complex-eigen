#include <stdio.h>
#include <math.h>
#include <values.h>
#include <sys/param.h>
#include "hessenberg.h"
#include "schur.h"
#include "eigen.h"

#define TAU 0.000001

////check values in https://www.real-statistics.com/linear-algebra-matrix-topics/hessenberg-decomposition/

static void matrix_display(const matrix *m);
static void vector_display(const vector *v);

int main(int argc, char* argv[])
{
    matrix a, eigenvec;
    vector eigenval;
    unsigned int code;
    const char *codestr;

    matrix_create(4, 4, &a);

    vector_init(&eigenval);
    matrix_init(&eigenvec);

    matrix_fill_int(&a,
        16, -1, 1, 2,
        2, 12, 1, -1,
        1, 3, -24, 2,
        4, -2, 1, 20
    );

    printf("Input Matrix\n");
    matrix_display(&a);

    code = eigen(&a, &eigenval, &eigenvec);
    codestr = eigen_err(code);
    printf("\nEigen status: %s\n", codestr);

    if(eigen_success(code))
    {
        printf("\nEigen values\n");
        vector_display(&eigenval);

        printf("\nEigen vectors (for respective eigen value column)\n");
        matrix_display(&eigenvec);
    }

    matrix_destroy(&a);
    matrix_destroy(&eigenvec);
    vector_destroy(&eigenval);

    return 0;
}

static void matrix_display(const matrix *m)
{
    int r, c;
    float f;
    for(r=0; r<m->rows; r++)
    {
        for(c=0; c<m->cols; c++) {
            f = MATRIXP(m, r, c).i;
            if(f!=0)
                printf("%s%.08f%s%.08fi", c>0 ? "\t" : "",
                    MATRIXP(m, r, c).r,
                    f<0 ? "-" : "+",
                    fabsf(f));
            else
                printf("%s%.08f", c>0 ? "\t" : "", MATRIXP(m, r, c).r);
        }
        printf("\n");
    }
}

static void vector_display(const vector *v)
{
    int i;
    float f;
    for(i=0; i<v->count; i++) {
        f = VECTORP(v,i).i;
        if(f!=0)
            printf("%s%.08f%s%.08fi",
                i>0 ? "\t" : "",
                VECTORP(v,i).r,
                f<0 ? "-" : "+",
                fabsf(f));
        else
            printf("%s%.08f",
                i>0 ? "\t" : "",
                VECTORP(v,i).r);
    }
    printf("\n");
}
