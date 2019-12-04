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

__global__ void preprocess_data(Argument_Pointers* arg, Constant_Coeffs* coeffs);

__device__ void Boundary_value(bool isU, DOUBLE t, int location, int location_extension, int width, int total_time,
	int* boundary_type, DOUBLE* hi, DOUBLE* boundary_array, DOUBLE* t_z, DOUBLE* boundary_condition, int* moc, int* dau, int* cuoi);


__device__ void FS_boundary(bool isU, DOUBLE t, int width, int total_time, int location, DOUBLE hmax,  
	int* boundary_type, DOUBLE* htaiz_bd, DOUBLE* FS, DOUBLE* CC, int* moc, int* dau, int*cuoi );

__global__ void Update_Boundary_Value(DOUBLE t, int total_time, Argument_Pointers* arg );

__global__ void update_uvz(Argument_Pointers* arg, Constant_Coeffs* coeffs);

__device__ void _normalize (DOUBLE coeff, int N, int M, int closenb_dist, int farnb_dist, DOUBLE* tmp_buffer, DOUBLE* val_buff, int* khouot);
__global__ void Normalize(DOUBLE isU, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);
__global__ void update_buffer(bool updateU, Argument_Pointers* arg, Array_Pointers* arr);

__device__ int locate_segment_v(int N, int M, bool* bienran1, bool* bienran2, int* first, int* last, int row, int col,  int* daui, int* cuoii, int* moci, DOUBLE* h, DOUBLE NANGDAY);

__device__ int locate_segment_u(int N, int M, bool* bienran1, bool* bienran2, int* first, int* last, int row, int col,  int* dauj, int* cuoij, int* mocj, DOUBLE* h, DOUBLE NANGDAY);
