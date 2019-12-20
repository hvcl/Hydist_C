
#include "engine.h"
#include "UVZSolver_multithread.h"
#include "support_funcs.h"
using namespace std;




void update_boundary_at_t(float t, Argument_Pointers* d_arg_ptr, Constant_Coeffs* coeffs, bool channel);
{
	if (channel){
		// int offset = M + 3;
		// boundary_up(t, d_arg_ptr, coeffs);
		// boundary_down(t, d_arg_ptr, coeffs);
		// boundary_left(t, d_arg_ptr, coeffs);
		// boundary_right(t, d_arg_ptr, coeffs);

	} else{
		int total_time = coeffs->total_time;
		Update_Boundary_Value(t, total_time,  )
	}
}


void Hydraulic_Calculation(int Tmax, Argument_Pointers* d_arg_ptr, Array_Pointers* d_arr_ptr, Constant_Coeffs* coeffs, Options ops){
	// note: blocksize in this case is fixed to be 1024 threads, can change later
	int blocksize = 1024;
	int M1 = ops.M + 3;
	int N1 = ops.N + 3;
	dim3 block_u = (min(M1, blocksize), 1, 1);
	dim3 grid_u = ( M1 / block_u.x + 1, 1, 1);
	dim3 block_v = (min(N1, blocksize), 1, 1);
	dim3 grid_v = (N1 / block_v.x, 1, 1);
	dim3 block_2d = (min(blocksize, M1), 1, 1);
	dim3 grid_2d = ((int) ceil((DOUBLE)(M1) / min(blocksize, M1)), N1, 1) ;
	int start_idx, end_idx, jump_step;
	bool channel = ops.kenhhepng xor ops.kenhhepd;
	int t = ops.t_start;
	while (t < Tmax){
		t += 0.5 * t;

		// update boundary conditionmake
		update_boundary_at_t(t, d_arg_ptr, coeffs, channel);

		// set start/ end index for kernels
		start_idx = 2;
		end_idx = M;
		jump_step = 2;
		if ((channel) && (ops.kenhhepd)){
			start_idx = 3;
			end_idx = M - 1;
		}
		// is U 
		// jump_step
		



		// set block size

		// call kernels



	}
}