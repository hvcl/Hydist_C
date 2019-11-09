// #include <cuda.h>
// #include <cuda_runtime.h>
#include <vector>
#include "Engine.h"
#include "constant.h"
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
	int hour(1), min(0), sec(0), sediment_start(20), bed_change_start(500), bed_change_end(1000);
	int save_interval(60), sediment_end(40);
	bool plot(true), visualize(true), debug(false), load_initial_condition(false);
	int c;


	while ((c = getopt(argc, argv, "f:h:m:s::p::d::v::b::B::e::E::i::c::")) != -1){
		switch (c)
		{
			case('f'):
				dir = optarg;
				continue;
			case ('h'):
				hour = atoi(optarg);
				continue;
			case ('m'):
				min = atoi(optarg);
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
		}
	}

	cout << "dir " << dir << endl
		<< "hour " << hour << endl
		<< "min " << min << endl
		<< "sec " << sec << endl
		<< "plot " << plot << endl
		<< "debug " << debug << endl
		<< "visualize " << visualize << endl
		<< "sediment start " << sediment_start << endl
		<< "bed_change_start " << bed_change_start << endl
		<< "sediment_end " << sediment_end << endl
		<< "save_interval " << save_interval << endl 
		<< "load initial condition " << load_initial_condition << endl; 

	// // 3. Call initialize functions
	// // a. load inputs

	// // test this first
	// pair<int, int> mesh_size = load_file<DOUBLE>("dosausongLUY.txt", h);
	// cout << mesh_size.first << " " << mesh_size.second << endl;

	
	pair <int, int> mesh_size =  load_inputs(dir, h, hsnham, VISCOIDX, Fw, bc_up, bc_down, bc_left, bc_right,
				CC_u, CC_d, CC_l, CC_r,  bienQ);
	cout << mesh_size.first << " " << mesh_size.second << endl;
	M = mesh_size.first;
	N = mesh_size.second;

	// b. alocate memory on GPU and transfer data to GPU

	// argument pointers struct's pointer on device, this is a pointer to device's memory
	Argument_Pointers* argument_pointer_struct_on_device;

	// a copy of argument pointer struct that is stored on the host
	// from this struct we can access pointers to arrays on device
	Argument_Pointers argument_pointer;
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
	// 	 << "size : " << host_ap.h.size() 
	// 	 << " " << h.size() << endl
	// 	 << "first 2 elem:  " << h[ M * 2 + 100] << " " << h[M * 2 + 100] << endl
	// 	 << " " << host_ap.h[M * 2 + 100] << " " << host_ap.h[M * 2 + 100] << endl;

	// bool fail = false;

	// for (int i = 0; i < h.size(); i++){
	// 	if (h[i] != 0){
	// 		cout << h[i] << endl; 
	// 		if (host_ap.h[i] != h[i]){
	// 			fail = true;
	// 			cout << "test failed"
	// 				<< h[i] << " != " << host_ap.h[i] << endl;
	// 			break;
	// 		}
	// 	}

	// }
	// if (!fail) 
	// 		cout << "sucess! " << endl;

	// what do I need: 
	// a struct that contain addresses of pointers in host: done: host_ap
	// a struct that contain addresses of pointer in device 

	// test function : attribute_arrays_memory_alloc
	argument_pointer = attribute_arrays_memory_alloc(0, host_ap, &argument_pointer_struct_on_device);

	DOUBLE* gpu_h; 
	gpu_h = (DOUBLE*) malloc(host_ap.h.size());
	cudaError_t status = cudaMemcpy((void*) gpu_h, argument_pointer.h, sizeof(DOUBLE) * host_ap.h.size(), cudaMemcpyDeviceToHost);
	assert(status == cudaSuccess);

	cout << "sucess!" << endl;

	for (int i = 0; i < host_ap.h.size(); i++){
		if (gpu_h[i] != host_ap.h[i]){
			cout << "failed" << endl
				<< "i = " << i << endl
				<< "gpu: " << gpu_h[i] << endl
				<< "cpu: " << host_ap.h[i] << endl;

		}
	}

	// check if values on device are the same with values on host, and if we has stored the right pointers






	// // c. calculate presequiste coeffs on GPU 
	// if (load_initial_condition)
	// 	// call load initial condition here


	// 
	



	// call initial function here 

	// 4. enter main loop
	// todo : write kernels on engine.h 











	return 0;
}
