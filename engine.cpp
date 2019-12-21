
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


void Hydraulic_Calculation(DOUBLE dT, DOUBLE NANGDAY, Argument_Pointers* d_arg_ptr, Array_Pointers* d_arr_ptr, Constant_Coeffs* coeffs, Options ops){
	// note: blocksize in this case is fixed to be 1024 threads, can change later
	int blocksize = 1024;
	int M = ops.M; 
	int N = ops.N;
	int M1 = M + 3;
	int N1 = N + 3;

	dim3 block_2d = (min(blocksize, M1), 1, 1);
	dim3 grid_2d = ((int) ceil((DOUBLE)(M1) / min(blocksize, M1)), N1, 1) ;

	dim3 block_shape;
	dim3 grid_shape;

	int start_idx, end_idx, jump_step;
	bool isU;

	bool channel = ops.channel;

	DOUBLE t = ops.t_start;
	// int Tmax = ops.Tmax;
	DOUBLE Tmax = 0.5;
	cout << "t = " << t << endl;
	while (t < Tmax){
		t += 0.5 * dT;
		cout << "Tmax = " << ops.Tmax << endl;

		// update boundary conditionmake
		update_boundary_at_t(M, N, t, channel, ops.total_time, d_arg_ptr, coeffs);
		synch_and_check();

		// set start/ end index for kernels
		start_idx = 2;
		end_idx = M;
		jump_step = 2;
		if ((channel) && (ops.kenhhepd)){
			start_idx = 3;
			end_idx = M - 1;
		}

		isU = true;
		jump_step = 2;

		// set block size
		block_shape = dim3(1, 1024, 1);
		grid_shape = dim3(M1, 1, 1);

		// block_shape = (1, 1024, 1);
		// grid_shape = (M, (int) ceil( N / 1024.0), 1);

		cout << "block_shape: " << block_shape.x << " " << block_shape.y <<  " " << block_shape.z << endl;
		cout << grid_shape.x << " " << grid_shape.y << " " << grid_shape.z << endl;
		cout << "M : " << M1 << " N: " << N1 << endl;

		UZSolver_calculate_preindex <<<grid_shape, block_shape>>> (start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
		synch_and_check();
		UZSolver_calculate_abcd<<<grid_shape, block_shape>>>(start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
		synch_and_check();
		// UZSolver_calculate_matrix_coeff<<<grid_shape, block_shape>>>(start_idx, end_idx, NANGDAY, d_arg_ptr, d_arr_ptr);
		// synch_and_check();

		// tridiagSolver<<<(1, M-1, 1), (32, 1, 1)>>> (false,isU, start_idx, end_idx, jump_step, 2 * N + 1, d_arg_ptr, d_arr_ptr);
		// synch_and_check();


		// UZSolver_extract_solution <<<grid_shape, block_shape>>>(start_idx, end_idx, NANGDAY, d_arg_ptr, d_arr_ptr);
		// synch_and_check();

		// Normalize<<<grid_2d, block_2d>>> (isU,d_arg_ptr, d_arr_ptr, coeffs);
		// synch_and_check();
		// update_buffer <<<grid_2d, block_2d>>>(isU, d_arg_ptr, d_arr_ptr);
		// synch_and_check();

		// solveV <<<grid_shape, block_shape>>>(t, 2, N, d_arg_ptr, coeffs);
		// synch_and_check();
		// update_margin_elem_V<<<(1, N, 1), (32, 1, 1)>>> (2, N, NANGDAY, d_arg_ptr);
		// synch_and_check();

		// // note that isU here is false since it normalize value of v after solving for v
		// Normalize<<<grid_2d, block_2d>>> (false,d_arg_ptr, d_arr_ptr, coeffs);
		// synch_and_check();
		// update_buffer <<<grid_2d, block_2d>>>(false, d_arg_ptr, d_arr_ptr);
		// synch_and_check();


		// update_h_moi <<<grid_2d, block_2d>>> (d_arg_ptr);
  //       synch_and_check();;


  //       update_uvz <<<grid_2d, block_2d>>> (d_arg_ptr, coeffs);
  //       synch_and_check();
   

  //       Find_Calculation_limits_Horizontal <<<(1, N, 1), (32, 1, 1)>>> (d_arg_ptr, coeffs);
  //       Find_Calculation_limits_Vertical <<<(1, M, 1), (32, 1, 1)>>>(d_arg_ptr, coeffs);
  //       Htuongdoi <<<grid_2d, block_2d>>> (d_arg_ptr);
  //       synch_and_check();

        // get result from device here and check


        // sediment transport simulation condition start here


        // second half of the simulation


  //       t += dT * 0.5;


		// update_boundary_at_t(M, N, t, channel, ops.total_time, d_arg_ptr, coeffs);        
  //       synch_and_check();

        

  //       block_shape = dim3(1024, 1, 1) ;
  //       grid_shape = dim3((int) (ceil(M / 1024.0)), N, 1);
  //       start_idx = 2;
  //       end_idx = N;
  //       jump_step = 2;
  //       isU = false;
  //       if ((ops.channel) && (ops.kenhhepng)){
  //           start_idx = 3;
  //           end_idx = N;
  //       }


  //       VZSolver_calculate_preindex <<<grid_shape,block_shape>>> (start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
  //       synch_and_check();
            
  //       VZSolver_calculate_abcd<<<grid_shape,block_shape>>> (start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
  //       synch_and_check();
  //       VZSolver_calculate_matrix_coeff<<<grid_shape,block_shape>>> (start_idx, end_idx, NANGDAY, d_arg_ptr, d_arr_ptr);
  //       synch_and_check();


  //       tridiagSolver<<<(1, N - 1, 1), (32, 1, 1)>>> (false, isU, start_idx, end_idx, jump_step, 2 * M + 1, d_arg_ptr, d_arr_ptr);
  //       synch_and_check();

  //       VZSolver_extract_solution<<<grid_shape,block_shape>>> (start_idx, end_idx, NANGDAY, d_arg_ptr, d_arr_ptr);
  //       synch_and_check();      


  //       Normalize<<<grid_2d, block_2d>>> (isU, d_arg_ptr, d_arr_ptr, coeffs);
  //       synch_and_check();
  //       update_buffer<<<grid_2d, block_2d>>> (isU, d_arg_ptr, d_arr_ptr);
  //       synch_and_check();

  //       solveU<<<grid_2d, block_2d>>> (t, 2, M, d_arg_ptr, coeffs);
  //       synch_and_check();
  //       update_margin_elem_U<<<(1, M, 1),(32, 1, 1)>>> (2, M, NANGDAY, d_arg_ptr);
  //       synch_and_check();

  //       // similar to first haft, isU here is true, since it normalize u value after solving for u
  //       Normalize<<<grid_2d, block_2d>>> (true, d_arg_ptr, d_arr_ptr, coeffs);
		// synch_and_check();
		// update_buffer <<<grid_2d, block_2d>>>(true, d_arg_ptr, d_arr_ptr );
		// synch_and_check();


		// update_h_moi <<<grid_2d, block_2d>>> (d_arg_ptr);
  //       synch_and_check();
  //       Reset_states_vertical <<<(M, 1, 1), (1, 32, 1)>>> (d_arg_ptr, coeffs);
  //       synch_and_check();;


  //       update_uvz <<<grid_2d, block_2d>>> (d_arg_ptr, coeffs);
  //       synch_and_check();
   

  //       Find_Calculation_limits_Horizontal <<<(1, N, 1), (32, 1, 1)>>> (d_arg_ptr, coeffs);
  //       Find_Calculation_limits_Vertical <<<(1, M, 1), (32, 1, 1)>>>(d_arg_ptr, coeffs);
  //       Htuongdoi <<<grid_2d, block_2d>>> (d_arg_ptr);
  //       synch_and_check();

       // sediment transport simulation here


	}
}




