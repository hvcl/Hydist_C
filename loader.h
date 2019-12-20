#include "engine.h"
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <cmath>

void Load_coeffs(Constant_Coeffs& ce);


template <typename T>
pair<int, int> load_file(const char* filename, vector<T> &arr);


pair<int, int> load_inputs(string dir, vector<DOUBLE> &h,
				vector<DOUBLE> &hsnham,
				vector<DOUBLE> &VISCOIDX, 
				vector <DOUBLE> &Fw,
				vector<DOUBLE> &bc_up,
				vector<DOUBLE> &bc_right,
				vector<DOUBLE> &bc_left, 
				vector<DOUBLE> &bc_down, 
				vector<DOUBLE> &CC_u, 
				vector<DOUBLE> &CC_d,
				vector<DOUBLE> &CC_l,
				vector<DOUBLE> &CC_r,
				vector<int> &bienQ,
				int* total_time);

// cudaError_t initialize(int device);
void load_initial_condition(string dir, 
							vector<DOUBLE> &u, 
							vector<DOUBLE> &v,
							vector<DOUBLE> &z,
							vector<DOUBLE> &FS,
							vector<int> &khouot);


template <typename T>
void device_copy(vector<T> &source, T* des);

template <typename T>
T* device_alloc(int nBytes);

template <typename T>
T* device_alloc_and_copy(vector<T> &h_array);

Array_Pointers supporting_arrays_alloc(int M, int N, Array_Pointers** device_arr_ptr);

Argument_Pointers attribute_arrays_memory_alloc(int device, Host_arrays &ap, Argument_Pointers** device_arg_ptr);
// void update_boundary()