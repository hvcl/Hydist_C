#include "cuda.h"
#include "engine.h"
#define powf pow


__global__ void Onetime_init( Argument_Pointers *arg, Constant_Coeffs* coeffs);

__global__ void update_h_moi(Argument_Pointers* arg);

__global__ void Reset_states_horizontal(Argument_Pointers* arg, Constant_Coeffs* coeffs);


__global__ void Reset_states_vertical(Argument_Pointers* arg, Constant_Coeffs* coeffs);


__device__ void Interpolate_FS_d(int location, int offset, int sign, Argument_Pointers* arg, Constant_Coeffs* coeffs);

__device__ void Interpolate_FS_ng(int location, int offset, int sign, Argument_Pointers* arg, Constant_Coeffs* coeffs);

__global__ void Find_Calculation_limits_Horizontal( Argument_Pointers *arg, Constant_Coeffs* coeffs);
__global__ void Find_Calculation_limits_Vertical(Argument_Pointers *arg, Constant_Coeffs* coeffs);



__global__ void Htuongdoi(Argument_Pointers* arg);
__global__ void boundary_up (DOUBLE t, Argument_Pointers* arg, Constant_Coeffs* coeffs);
__global__ void boundary_down(DOUBLE t, Argument_Pointers* arg, Constant_Coeffs* coeffs);
__global__ void boundary_left(DOUBLE t, Argument_Pointers* arg, Constant_Coeffs* coeffs);
__global__ void boundary_right(DOUBLE t, Argument_Pointers* arg, Constant_Coeffs* coeffs);