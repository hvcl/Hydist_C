
#include "engine.h"
#include "UVZSolver_multithread.h"
#include "support_funcs.h"
#include <iostream>
using namespace std;


void synch_and_check(){
	cudaDeviceSynchronize();
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess)
		cout << "error: " <<  cudaGetErrorString(err) << endl;
}


void update_boundary_at_t(int M, int N, float t, bool channel, int total_time, Argument_Pointers* d_arg_ptr, Constant_Coeffs* coeffs)
{
	if (channel){
		// int offset = M + 3;
		// boundary_up(t, d_arg_ptr, coeffs);
		// boundary_down(t, d_arg_ptr, coeffs);
		// boundary_left(t, d_arg_ptr, coeffs);
		// boundary_right(t, d_arg_ptr, coeffs);

	} else{
		Update_Boundary_Value<<<(1, 1024, 1), (1, max(M, N) / 1024 + 1)>>>(t, total_time, d_arg_ptr);
	}
}


void Hydraulic_Calculation(Argument_Pointers* d_arg_ptr, Array_Pointers* d_arr_ptr, Constant_Coeffs* coeffs, Options ops){
	// note: blocksize in this case is fixed to be 1024 threads, can change later
	int blocksize = 1024;
	int M = ops.M; 
	int N = ops.N;
	int M1 = M + 3;
	int N1 = N + 3;
	dim3 block_u = (min(M1, blocksize), 1, 1);
	dim3 grid_u = ( M1 / block_u.x + 1, 1, 1);
	dim3 block_v = (min(N1, blocksize), 1, 1);
	dim3 grid_v = (N1 / block_v.x, 1, 1);
	dim3 block_2d = (min(blocksize, M1), 1, 1);
	dim3 grid_2d = ((int) ceil((DOUBLE)(M1) / min(blocksize, M1)), N1, 1) ;
	int start_idx, end_idx, jump_step;
	bool channel = ops.kenhhepng xor ops.kenhhepd;
	int t = ops.t_start;
	// int Tmax = ops.Tmax;
	int Tmax = 0.5;
	cout << "t = " << t << endl;
	while (t < Tmax){
		t += 0.5 * t;
		cout << "Tmax = " << ops.Tmax << endl;

		// update boundary conditionmake
		update_boundary_at_t(M, N, t, channel, ops.total_time, d_arg_ptr, coeffs);
		synch_and_check();

		// set start/ end index for kernels
		// start_idx = 2;
		// end_idx = M;
		// jump_step = 2;
		// if ((channel) && (ops.kenhhepd)){
		// 	start_idx = 3;
		// 	end_idx = M - 1;
		// }
		// is U 
		// jump_step
		

		// set block size

		// call kernels



	}
}




