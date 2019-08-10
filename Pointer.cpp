#include "Pointer.h"
using namespace std;

/*
Constructor: - allocate memory on device according to arrays on host
			 - allocate supporting memory on device for tridiag solver and for calculating parameters
*/

Pointer::Pointer(map<string, pointer_struct> &ptr_list, int M, int N)
{
	// allocate memory for arrays inputs from host
	for (auto it = ptr_list.begin(); it!=ptr_list.end(); it++) {
		double* d_ptr;
		int size = it->second.size;
		string key = it->first;
		cudaError_t err = cudaMalloc((void**)d_ptr, size);
		// TODO:  make sure err is success here
		device_pointers[key] = pointer_struct(d_ptr, size);
	}

	// allocate memory for arrays that reside on device only
	string tridiag_ptr[] = {"x", "AA", "BB", "CC", "DD", "Ap", "Bp", "ep"};
	int name_list_size = *(&tridiag_ptr + 1) - tridiag_ptr;
	int array_size = 2 * M * N + 4 * max(M, N);
	double * d_ptr;
	for (int i = 0; i < name_list_size; i++) {
		// check error here
		cudaMalloc((void**)d_ptr, array_size);
		device_pointers[tridiag_ptr[i]] = pointer_struct(d_ptr, array_size);
	}

	string aux_ptr[] = {"a1", "b1", "c1", "d1", "a2", "c2", "d2", "f1", "f2", "f3", "f5"};
	name_list_size = *(&aux_ptr + 1) - aux_ptr;
	array_size = M * N + 2 * max(M, N);
	for (int i = 0; i < name_list_size; i++) {
		// check error here
		cudaMalloc((void**)d_ptr, array_size);
		device_pointers[aux_ptr[i]] = pointer_struct(d_ptr, array_size);
	}
}


Pointer::~Pointer()
{
}

/*  This function copy arrays on host to device
*/
void Pointer::to_device() {
	for (auto it = host_pointers.begin(); it != host_pointers.end(); it++) {
		double* h_ptr = it->second.ptr;
		double* d_ptr = device_pointers[it->first].ptr;
		int size = it->second.size;
		// check error here
		cudaMemcpy(d_ptr, h_ptr, size, cudaMemcpyHostToDevice);
	}
}


// This function extract data from device to host
void Pointer::extract(pointer_struct* arg_list[], int size)
{
	for (int i = 0; i < size; i++) {
		string key = find_value(arg_list[i]->ptr);
		cudaMemcpy(arg_list[i]->ptr, get_device_pointer(key), arg_list[i]->size, cudaMemcpyDeviceToHost);
	}
}

// this function get the pointer on device of the according key
double * Pointer::get_device_pointer(string key)
{
	return device_pointers[key].ptr;
}

