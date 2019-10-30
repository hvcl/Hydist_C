// #include "cuda.h"
#include "constant.h"
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <sstream>
#define DOUBLE double
// #include "cuda.h"
// #include "cuda_runtime.h"
using namespace std;


struct Options{
	int Tmax, interval, sediment_start, bed_change_start, kenhhepd, kenhhepng;
	bool ketdinh, channel, debug, plot;
	Options(int T, int itv, int sds, int bcs, int cohesive, int kn,
			bool kenhd, bool db, bool plt){
		Tmax = T;
		interval = itv;
		sediment_start = sds;
		bed_change_start = bcs;
		kenhhepng = kn;
		kenhhepd = kenhd;
		ketdinh = cohesive;
		channel = kenhd || kn;
		debug = db;
		plot = plt;
	}
	Options(){

	}
};



struct  Host_arrays
{
	int M, N, total_time;
	vector<DOUBLE> h, hsnham, VISCOINDX, Fw, FS, dH;
	vector<DOUBLE> bc_up, bc_down, bc_left, bc_right;
	vector<DOUBLE> CC_u, CC_d, CC_l, CC_r;
	vector<DOUBLE> u, v, z;
	vector<int> khouot, bienQ, boundary_type;
};


template <typename T>
pair<int, int> load_file(const char* filename, vector<T> &arr);


pair<int, int> load_inputs(string dir, vector<DOUBLE> &h,
				vector<DOUBLE> &hsnham,
				vector<DOUBLE> &VISCOINDX, 
				vector<DOUBLE> &bc_up,
				vector<DOUBLE> &bc_right,
				vector<DOUBLE> &bc_left, 
				vector<DOUBLE> &bc_down, 
				vector<DOUBLE> &CC_u, 
				vector<DOUBLE> &CC_d,
				vector<DOUBLE> &CC_l,
				vector<DOUBLE> &CC_r,
				vector<int> &bienQ);

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

Array_Pointers supporting_arrays_alloc(int M, int N, Array_Pointers** device_arr_ptr)

// Argument_Pointers attribute_arrays_memory_alloc(int device, Host_arrays &ap, Argument_Pointers** device_arg_ptr)
// void update_boundary()

// void simulate(Options ops);

