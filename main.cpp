// #include <cuda.h>
// #include <cuda_runtime.h>
#include <vector>
#include "engine.h"
#include "support_funcs.h"
#include "loader.h"
#include <algorithm>
#include <unistd.h>
#include <stdlib.h>
#include <utility>
#include <string>
#include <iostream>
#define DOUBLE double 

using namespace std;

// load input

// initialize environment


// launch kernel


int main (int argc, char ** argv){
	int M, N, total_time;
	vector<DOUBLE> h, hsnham, VISCOIDX, Fw, FS, dH;
	vector<DOUBLE> bc_up, bc_down, bc_left, bc_right;
	vector<DOUBLE> CC_u, CC_d, CC_l, CC_r;
	vector<DOUBLE> u, v, z;
	vector<int> khouot, bienQ, boundary_type;

	// 1. parse arguments 
	// a. input directory 
		/* -f: file 
	   -h: hour
	   -m: min
	   -s: second
	   -p: plot
	   -d: debug
	   -v: visualize
	   -b: sediment calculation begin
	   -B: bed change calculaition Begin
	   -e: sediment end
	   -E: bedchange end
	   -i : save interval - in minute
	*/
	string dir;
	int hour(1), minute(0), sec(0), sediment_start(20), bed_change_start(500), bed_change_end(1000);
	int save_interval(60), sediment_end(40);
	int t_start(0);
	bool kenhhepd(0), kenhhepng(0), cohesive(1);
	bool plot(true), visualize(true), debug(false), load_initial_condition(false);
	int c;


	while ((c = getopt(argc, argv, "f:h:m:s::p::d::v::b::B::e::E::i::c::T::k::n::")) != -1){
		switch (c)
		{
			case('f'):
				dir = optarg;
				continue;
			case ('h'):
				hour = atoi(optarg);
				continue;
			case ('m'):
				minute = atoi(optarg);
				continue;
			case ('s'): 
				sec = atoi(optarg);
				continue;
			case ('p'):
				plot = (bool) atoi(optarg);
				continue;
			case ('d'):
				debug = (bool) atoi(optarg);
				continue;
			case ('v'):
				visualize = (bool) atoi(optarg);
				continue;
			case ('b'):
				sediment_start = atoi(optarg);
				continue;
			case ('B'):
				bed_change_start = atoi(optarg);
				continue;
			case ('e'):
				sediment_end = atoi(optarg);
				continue;
			case ('E'):
				bed_change_end = atoi(optarg);
				continue;
			case ('i'):
				save_interval = atoi(optarg);
				continue;
			case ('c'):
				load_initial_condition = atoi(optarg);
				continue;
			case ('T'):
				t_start = atoi(optarg);
				continue;
			case('k'):
				kenhhepd = atoi(optarg);
				continue;
			case('n'):
				kenhhepng = atoi(optarg);
				continue;
		}
	}

	int Tmax = hour * 3600 + minute * 60 + sec;
	
	Options ops(Tmax, t_start, save_interval, sediment_start, bed_change_start, cohesive, kenhhepng, kenhhepd, debug, plot);

	// cout << "dir " << dir << endl
	// 	<< "hour " << hour << endl
	// 	<< "min " << minute << endl
	// 	<< "sec " << sec << endl
	// 	<< "plot " << plot << endl
	// 	<< "debug " << debug << endl
	// 	<< "visualize " << visualize << endl
	// 	<< "sediment start " << sediment_start << endl
	// 	<< "bed_change_start " << bed_change_start << endl
	// 	<< "sediment_end " << sediment_end << endl
	// 	<< "save_interval " << save_interval << endl 
	// 	<< "load initial condition " << load_initial_condition << endl; 

	// // 3. Call initialize functions
	// // a. load inputs

	// // test this first
	// pair<int, int> mesh_size = load_file<DOUBLE>("dosausongLUY.txt", h);
	// cout << mesh_size.first << " " << mesh_size.second << endl;

	
	pair <int, int> mesh_size =  load_inputs(dir, h, hsnham, VISCOIDX, Fw, bc_up, bc_down, bc_left, bc_right,
				CC_u, CC_d, CC_l, CC_r,  bienQ, boundary_type, &total_time);
	cout << mesh_size.first << " " << mesh_size.second << endl;


	// ATENDTION: input files that store 2d meshes must be preprocess so that the coordinates align with python indexing
	// i.e bottom left origin must be tranpose, flip and pad to make the padded array in top left origin 
	N = mesh_size.first - 3;
	M = mesh_size.second - 3;
	ops.M = M;
	ops.N = N;
	ops.total_time = total_time;


	// b. alocate memory on GPU and transfer data to GPU

	// argument pointers struct's pointer on device, this is a pointer to device's memory
	Argument_Pointers* d_argument_pointers;

	// a copy of argument pointer struct that is stored on the host
	// from this struct we can access pointers to arrays on device
	Argument_Pointers h_argument_pointers;

	Array_Pointers* d_arr_pointers, h_arr_pointers; 
	// host pointers here 
	Host_arrays host_ap; 

	host_ap.M = M;
	host_ap.N = N;
	host_ap.h = h;
	host_ap.hsnham = hsnham;
	host_ap.VISCOIDX = VISCOIDX;
	host_ap.Fw = Fw;
	host_ap.FS = FS;
	host_ap.dH = dH;
	host_ap.bc_up = bc_up;
	host_ap.bc_down = bc_down;
	host_ap.bc_left = bc_left;
	host_ap.bc_right = bc_right;
	host_ap.CC_u = CC_u;
	host_ap.CC_d = CC_d;
	host_ap.CC_r = CC_r;
	host_ap.CC_l = CC_l;
	host_ap.u = u;
	host_ap.v = v;
	host_ap.z = z;
	host_ap.khouot = khouot;
	host_ap.bienQ = bienQ;
	host_ap.boundary_type = boundary_type;



	// cout << "checking vector assignment: " << endl
	// 	 << "size in GPU : " << host_ap.bienQ.size() 
	// 	 << "  size in CPU: " << bienQ.size() << endl;
	// 	 // << "first 2 elem:  " << boundary_type[ M * 2 + 100] << " " << boundary_type[M * 2 + 100] << endl
	// 	 // << " " << host_ap.boundary_type[M * 2 + 100] << " " << host_ap.boundary_type[M * 2 + 100] << endl;

	// bool fail = false;

	// for (int i = 0; i < bienQ.size(); i++){
	// 	if (bienQ[i] != 0){
	// 		cout << bienQ[i] << endl; 
	// 		if (host_ap.bienQ[i] != bienQ[i]){
	// 			fail = true;
	// 			cout << "test failed"
	// 				<< bienQ[i] << " != " << host_ap.bienQ[i] << endl;
	// 			break;
	// 		}
	// 	}

	// }
	// if (!fail) 
	// 		cout << "sucess! " << endl;
	DOUBLE* gpu_h; 
	gpu_h = (DOUBLE*) malloc(sizeof(DOUBLE) * host_ap.h.size());

	// test function : attribute_arrays_memory_alloc
	h_argument_pointers = attribute_arrays_memory_alloc(0, host_ap, &d_argument_pointers);

	// gpu_h: the grid depth map on gou
	
	// cout << host_ap.h.size() << endl;
	cudaError_t status = cudaMemcpy((void*) gpu_h, h_argument_pointers.h, sizeof(DOUBLE) * host_ap.h.size(), cudaMemcpyDeviceToHost);
	assert(status == cudaSuccess);

	h_arr_pointers = supporting_arrays_alloc(M, N, &d_arr_pointers);


	// load coefficients used in Device code
	Constant_Coeffs  h_const_coeffs, *d_const_coeffs;
	Load_coeffs (h_const_coeffs);

	cudaMalloc((void**) &d_const_coeffs, sizeof(Constant_Coeffs));
	status = cudaMemcpy(d_const_coeffs, &h_const_coeffs, sizeof(Constant_Coeffs), cudaMemcpyHostToDevice);
	assert(status == cudaSuccess);

	// check struct values on GPU:
	Constant_Coeffs* coeffs;
	coeffs = (Constant_Coeffs*) malloc(sizeof(Constant_Coeffs));
	status = cudaMemcpy((void*) coeffs, d_const_coeffs, sizeof(Constant_Coeffs), cudaMemcpyDeviceToHost);
	assert(status == cudaSuccess);

	// cout << coeffs->dX << "  " << coeffs->dY << " " << coeffs->Ks << endl;
	// cout << "data transfer Done" << endl;

	// check if values on device are the same with values on host, and if we has stored the right pointers
	// done
	
	for (int i = 0; i < host_ap.h.size(); i++){
		if (gpu_h[i] != host_ap.h[i]){
			cout << "failed" << endl
				<< "i = " << i << endl
				<< "gpu: " << gpu_h[i] << endl
				<< "cpu: " << host_ap.h[i] << endl;

		}
	}

	// grid size and block size 
	dim3 block_2d(min(32, M + 3), 1, 1);
	dim3 grid_2d((M + 3) / min(32, M + 3) + 1,N + 3, 1);
	// dim3 grid_2d(4, 30, 1);
	cout << grid_2d.x << " " << grid_2d.y <<" " << grid_2d.z <<  " " << endl;
	Onetime_init <<<grid_2d, block_2d >>>(d_argument_pointers,d_const_coeffs);

 	
	synch_and_check();


	int* arr;
	arr = (int*) malloc(sizeof(int) * host_ap.h.size());
	status = cudaMemcpy((void*) arr, h_argument_pointers.khouot, sizeof(int) * host_ap.h.size(), cudaMemcpyDeviceToHost);
	assert(status == cudaSuccess);
	ofstream ofs;
	ofs.open("Outputs/khouot.txt");
	for (int i= 0; i < N + 3; i ++){
			for (int j = 0; j < M + 3; j ++){
				ofs << arr[i * (M + 3) + j] << " ";
	
			}
		ofs << endl;
	}



	// load initial condition
	if (load_initial_condition)
		printf("need to write load initial condition here\n");

	Find_Calculation_limits_Horizontal <<<(1, N, 1), (32, 1, 1)>>> (d_argument_pointers, d_const_coeffs);

	Find_Calculation_limits_Vertical <<<(1, M, 1), (32, 1, 1)>>> (d_argument_pointers, d_const_coeffs);


	Htuongdoi<<<grid_2d, block_2d>>> (d_argument_pointers);

	synch_and_check();

	preprocess_data <<<(1, 1, 1), (32, 1,1 )>>> (d_argument_pointers, d_const_coeffs);
	synch_and_check();


	// enter main loop here

	Hydraulic_Calculation(d_argument_pointers, d_arr_pointers, d_const_coeffs, ops);



	return 0;
}
