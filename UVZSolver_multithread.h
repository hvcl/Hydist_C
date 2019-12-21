#include "engine.h"


__device__ void tridiag(int sn, DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD, DOUBLE *x, DOUBLE *Ap, DOUBLE *Bp, DOUBLE *ep);

__global__ void  tridiagSolver(bool print, bool isU, int startidx, int endidx, int jumpstep, int tridiag_coeff_width, Argument_Pointers* arg, Array_Pointers * arr);




__device__ void bienrandau(int i, int first, int last,  DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD,
    DOUBLE *a1, DOUBLE *b1, DOUBLE *c1, DOUBLE *d1, DOUBLE *a2, DOUBLE *c2, DOUBLE *d2);

__device__ void bienlongdau(int i, int first, int last,  DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD,
    DOUBLE *a1, DOUBLE *b1, DOUBLE *c1, DOUBLE *d1, DOUBLE *a2, DOUBLE *c2, DOUBLE *d2);

__device__ void  _vzSolver_calculate_preindex(DOUBLE t, int i, int j, int width, int first, int last,  Argument_Pointers* arg, Array_Pointers *arr, Constant_Coeffs* coeffs);

__device__ void  _uzSolver_calculate_preindex(int i, int j, int width, int first, int last, Argument_Pointers* arg, Array_Pointers *arr, Constant_Coeffs* coeffs);


__device__ void _calculate_abcd(int i, int j, int first, int last, DOUBLE f4, int support_array_width,  bool bienran1, bool bienran2, Array_Pointers* arr);

__device__ void _calculate_matrix_coeff(bool isU, int i, int j, int support_array_width, int tridiag_coeff_width, int first, int last, int seg_no, bool bienran1, bool bienran2, 
    DOUBLE ubt_or_vbd, DOUBLE ubp_or_vbt,  DOUBLE TZ_r, DOUBLE TZ_l, bool dkBienQ_1, bool dkBienQ_2, int dkfr, Argument_Pointers* arg, Array_Pointers *arr);




__device__ void _vzSolver_extract_solution( int i,int j, int sn, int width, int first, int last, bool bienran1, bool bienran2, Argument_Pointers *arg, Array_Pointers * arr);

__device__ void _uzSolver_extract_solution( int i, int j, int sn, int width, int first, int last, bool bienran1, bool bienran2, Argument_Pointers *arg, Array_Pointers * arr);

__device__ void vSolver(DOUBLE t, int offset, int first, int last, int row, int col, bool bienran1, bool bienran2, DOUBLE* VISCOIDX, DOUBLE* Tsyw, 
    DOUBLE *v, DOUBLE *t_v, DOUBLE *u, DOUBLE *t_u, DOUBLE *z, DOUBLE *t_z, DOUBLE *Ky1, DOUBLE *Htdv, DOUBLE *H_moi, Constant_Coeffs* coeffs);

__device__ void uSolver(DOUBLE t, int offset, int N, int first, int last, int row, int col, bool bienran1, bool bienran2, DOUBLE* VISCOIDX, DOUBLE* Tsxw,
    DOUBLE *v, DOUBLE *t_v, DOUBLE *u, DOUBLE *t_u, DOUBLE *z, DOUBLE *t_z, DOUBLE *Kx1, DOUBLE *Htdu, DOUBLE *H_moi, Constant_Coeffs* coeffs);



__global__ void VZSolver_calculate_preindex(int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);
__global__ void VZSolver_calculate_abcd(int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);
__global__ void VZSolver_calculate_matrix_coeff(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg, Array_Pointers* arr);
__global__ void VZSolver_extract_solution(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg, Array_Pointers* arr);
__global__ void solveV(DOUBLE t, int startidx, int endidx, Argument_Pointers* arg, Constant_Coeffs* coeffs);
__global__ void update_margin_elem_V(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg);

__global__ void UZSolver_extract_solution(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg, Array_Pointers* arr);
__global__ void UZSolver_calculate_preindex(int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);
__global__ void UZSolver_calculate_abcd(int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);
__global__ void UZSolver_calculate_matrix_coeff(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg, Array_Pointers* arr);
__global__ void solveU(DOUBLE t, int startidx, int endidx, Argument_Pointers* arg, Constant_Coeffs* coeffs);
__global__ void update_margin_elem_U(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg);
