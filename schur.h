#ifndef _SCHUR_
#define _SCHUR_

#include "complex_matrix.h"

#define SCHUR_CODE_OK 100
#define SCHUR_CODE_TIMEOUT 101

#define SCHUR_CODE_MIN SCHUR_CODE_OK
#define SCHUR_CODE_MAX SCHUR_CODE_TIMEOUT

extern unsigned int schur(matrix *h, matrix *q, int maxiter);
extern const char* schur_err(unsigned int err);

#endif //_SCHUR_
