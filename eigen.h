#ifndef _EIGEN_
#define _EIGEN_

#include <stdbool.h>
#include "complex_matrix.h"

extern unsigned int eigen(const matrix *m, vector *eigenval, matrix *eigenvec);
extern const char* eigen_err(unsigned int code);
extern bool eigen_success(unsigned int code);

#endif //_EIGEN_
