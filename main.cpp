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

	
	pair <int, int> mesh_size =  load_inputs(dir, h, hsnham, VISCOIDX, bc_up, bc_down, bc_left, bc_right,
				CC_u, CC_d, CC_l, CC_r,  bienQ);
	cout << mesh_size.first << " " << mesh_size.second << endl;
	M = mesh_size.first;
	N = mesh_size.second;

	// // b. alocate memory on GPU and transfer data to GPU
	// // argument pointers struct's pointer on device
	// Argument_Pointers* argument_pointer_struct_on_device;
	// Argument_Pointers argument_pointer;
	// // host pointers here 
	


	// // allocate memory and transfer here
	// // todo: test the function 




	// // c. calculate presequiste coeffs on GPU 
	// if (load_initial_condition)
	// 	// call load initial condition here





	// call initial function here 

	// 4. enter main loop
	// todo : write kernels on engine.h 











	return 0;
}
