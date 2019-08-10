// #include <cuda.h>
// #include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include "Engine.h"


using namespace std;

// load input

// initialize environment


// launch kernel


int main (int* argc, char** argv[]){
	int bienQ[4];
	int M, N, total_time
	vector<DOUBLE> h, hsnham, VISCOINDX, Fw, FS, dH;
	vector<DOUBLE> bc_up, bc_down, bc_left, bc_right;
	vector<DOUBLE> CC_u, CC_d, CC_l, CC_r;
	vector<DOUBLE> u, v, z;
	vector<int> khouot;

	// extract parameter

	// call load input accordingly
	// load_inputs(h, hsnham, VISCOINDX, bc_up, bc_down, bc_left, bc_right,
	// 			CC_u, CC_d, CC_l, CC_r, , bienQ);
	pair<int, int> mesh_size;
	mesh_size = load_file("dosausongLUY.txt");
	cout << mesh_size.first << " " << mesh_size.second << endl;

	// call intialize environment




	retun 0;
}