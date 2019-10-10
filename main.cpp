// #include <cuda.h>
// #include <cuda_runtime.h>
#include <iostream>
#include <vector>
// #include "Engine.h"
#include <unistd.h>


using namespace std;

// load input

// initialize environment


// launch kernel


int main (int* argc, char** argv[]){
	int bienQ[4];
	int M, N, total_time;
	vector<DOUBLE> h, hsnham, VISCOINDX, Fw, FS, dH;
	vector<DOUBLE> bc_up, bc_down, bc_left, bc_right;
	vector<DOUBLE> CC_u, CC_d, CC_l, CC_r;
	vector<DOUBLE> u, v, z;
	vector<int> khouot;

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
	string f;
	int hour(1), min(0), sec(0), sediment_start(20), bed_change_start(500), bed_change_end(1000);
	int save_interval(60), sediment_end(40);
	bool plot(true), visualize(true), debug(false), load_initial_condition(false);

	switch (getopt(argc, argv, "f:h:m:s::p::d::v::b::B::e::E::i::c"))
	{
		case(f):
			dir = optarg;
			continue;
		case (h):
			hour = atoi(optarg);
			continue;
		case (m):
			min = atoi(optarg);
			continue;
		case (s): 
			sec = atoi(optarg);
			continue;
		case (p):
			plot = (bool) atoi(optarg);
			continue;
		case (d):
			debug = (bool) atoi(optarg);
			continue;
		case (v):
			visualize = (bool) atoi(optarg);
			continue;
		case (b):
			sediment_start = atoi(optarg);
			continue;
		case (B):
			bed_change_start = atoi(optarg);
			continue;
		case (e):
			sediment_end = atoi(optarg);
			continue;
		case (E):
			bed_change_end = atoi(optarg);
			continue;
		case (i):
			save_interval = atoi(optarg);
			continue;
		case (c):
			load_initial_condition = atoi(optarg);
			continue;


	}
	// 2. initialize cuda

	// // 3. Call initialize functions
	// // a. load inputs
	// // load_inputs(h, hsnham, VISCOINDX, bc_up, bc_down, bc_left, bc_right,
	// // 			CC_u, CC_d, CC_l, CC_r, , bienQ);
	// pair<int, int> mesh_size;
	// // test this first
	// mesh_size = load_file("dosausongLUY.txt");
	// cout << mesh_size.first << " " << mesh_size.second << endl;
	// // if okay then 
	
	// mesh_size = // load_inputs(dir, h, hsnham, VISCOINDX, bc_up, bc_down, bc_left, bc_right,
	// // 			CC_u, CC_d, CC_l, CC_r,  bienQ);
	// cout << mesh_size.first << " " << mesh_size.second << endl;
	



	// M = mesh_size.first;
	// N = mesh_size.second;

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












	retun 0;
}
