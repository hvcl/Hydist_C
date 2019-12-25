#include "engine.h"





__global__ void hesoK(Constant_Coeffs* coeffs, Argument_Pointers* arg);


__global__ void Find_VTH(Constant_Coeffs* coeffs, Argument_Pointers* arg);

__global__ void Scan_FSi(DOUBLE t, DOUBLE s_start, bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);

__global__ void FSi_extract_solution( bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs*coeffs);

__global__ void Scan_FSj(DOUBLE t, DOUBLE s_start, bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);

__global__ void FSj_extract_solution(bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);

__global__ void Update_FS(Argument_Pointers* arg);

__global__ void BedLoad(DOUBLE t, bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs);