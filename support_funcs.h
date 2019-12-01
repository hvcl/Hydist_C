#include "cuda.h"
#include "engine.h"
#define powf pow


__global__ void Onetime_init( Argument_Pointers *arg, Constant_Coeffs* coeffs);