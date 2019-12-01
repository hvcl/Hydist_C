#include "engine.h"


__device__ void tridiag(int sn, DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD, DOUBLE *x, DOUBLE *Ap, DOUBLE *Bp, DOUBLE *ep);

__global__ void  tridiagSolver(bool print, bool isU, int startidx, int endidx, int jumpstep, int tridiag_coeff_width, Argument_Pointers* arg, Array_Pointers * arr);
__device__ void bienrandau(int i, int first, int last,  DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD,
    DOUBLE *a1, DOUBLE *b1, DOUBLE *c1, DOUBLE *d1, DOUBLE *a2, DOUBLE *c2, DOUBLE *d2);

__device__ void bienlongdau(int i, int first, int last,  DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD,
    DOUBLE *a1, DOUBLE *b1, DOUBLE *c1, DOUBLE *d1, DOUBLE *a2, DOUBLE *c2, DOUBLE *d2);

__device__ void  _vzSolver_calculate_preindex(DOUBLE t, int i, int j, int width, int first, int last,  Argument_Pointers* arg, Array_Pointers *arr, Constant_Coeffs* coeffs);