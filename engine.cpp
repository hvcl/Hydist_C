
#include "engine.h"
#include "UVZSolver_multithread.h"
#include "support_funcs.h"
#include "sediment_transport.h"
#include "loader.h"
#include <iostream>
#include <fstream>
#include <string>
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
		dim3 grid(1, max(M, N) / 1024 + 1);
		dim3 block(1, 1024);
		Update_Boundary_Value<<<grid, block>>>(t, total_time, d_arg_ptr);
	}
}




void Hydraulic_Calculation(DOUBLE dT, DOUBLE NANGDAY, Argument_Pointers* d_arg_ptr, Array_Pointers* d_arr_ptr, 
						Constant_Coeffs* coeffs, Argument_Pointers h_arg_pointer, Options ops){
	// note: blocksize in this case is fixed to be 1024 threads, can change later
	int blocksize = 1024;
	int M = ops.M; 
	int N = ops.N;
	int M1 = M + 3;
	int N1 = N + 3;

	dim3 block_2d(min(blocksize, M1), 1, 1);
	dim3 grid_2d((int) ceil((DOUBLE)(M1) / min(blocksize, M1)), N1, 1) ;

	dim3 block_shape;
	dim3 grid_shape;

	int start_idx, end_idx, jump_step;
	bool isU;

	bool channel = ops.channel;

	DOUBLE t = ops.t_start;
	int Tmax = ops.Tmax;
	// DOUBLE Tmax = 0.5;
	cout << "t = " << t << " sediment start = " << ops.sediment_start << endl;
	while (t < Tmax){
		t += 0.5 * dT;
// 		cout << t << endl;

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
		block_shape = dim3(1, 512, 1);
		grid_shape = dim3(M1, 1, 1);



		UZSolver_calculate_preindex <<<grid_shape, block_shape>>> (start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
		synch_and_check();
		UZSolver_calculate_abcd<<<grid_shape, block_shape>>>(start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
		synch_and_check();
		UZSolver_calculate_matrix_coeff<<<grid_shape, block_shape>>>(start_idx, end_idx, NANGDAY, d_arg_ptr, d_arr_ptr);
		synch_and_check();

		dim3 grid(1, M - 1, 1);
        tridiagSolver_v2<<<1, M - 1>>> (false,isU, start_idx, end_idx, jump_step, 2 * N + 1, d_arg_ptr, d_arr_ptr);
		// tridiagSolver<<<grid, 32>>> (false,isU, start_idx, end_idx, jump_step, 2 * N + 1, d_arg_ptr, d_arr_ptr);
		synch_and_check();


		UZSolver_extract_solution <<<grid_shape, block_shape>>>(start_idx, end_idx, NANGDAY, d_arg_ptr, d_arr_ptr);
		synch_and_check();



		Normalize<<<grid_2d, block_2d>>> (isU,d_arg_ptr, d_arr_ptr, coeffs);
		synch_and_check();
		update_buffer <<<grid_2d, block_2d>>>(isU, d_arg_ptr, d_arr_ptr);
		synch_and_check();

		solveV <<<grid_shape, block_shape>>>(t, 2, N, d_arg_ptr, coeffs);
		synch_and_check();

		grid = dim3(1, N, 1);
		update_margin_elem_V<<<grid, 32>>> (2, N, NANGDAY, d_arg_ptr);
		synch_and_check();

		// note that isU here is false since it normalize value of v after solving for v
		Normalize<<<grid_2d, block_2d>>> (false,d_arg_ptr, d_arr_ptr, coeffs);
		synch_and_check();
		update_buffer <<<grid_2d, block_2d>>>(false, d_arg_ptr, d_arr_ptr);
		synch_and_check();


		update_h_moi <<<grid_2d, block_2d>>> (d_arg_ptr);
        synch_and_check();

        
        update_uvz <<<grid_2d, block_2d>>> (d_arg_ptr, coeffs);
        synch_and_check();
   
        grid = dim3(1, N, 1);
        Find_Calculation_limits_Horizontal <<<grid, 32>>> (d_arg_ptr, coeffs);
        grid = dim3(1, M, 1);
        Find_Calculation_limits_Vertical <<<grid, 32>>>(d_arg_ptr, coeffs);
        Htuongdoi <<<grid_2d, block_2d>>> (d_arg_ptr);
        synch_and_check();

        // get result from device here and check


        if ( ((int) t % ops.interval == 0) && (t - (int) t == 0))
        	save_result(h_arg_pointer, (int) t);
        // sediment transport simulation condition start here

        if (t >= ops.sediment_start){
            Find_VTH<<<grid_2d, block_2d>>>(coeffs, d_arg_ptr);
            hesoK<<<grid_2d, block_2d>>> (coeffs, d_arg_ptr);
            synch_and_check();
         
            start_idx = 3;
            end_idx = M - 1;
            Scan_FSj<<<grid_shape, block_shape>>>(t, ops.bed_change_start, ops.ketdinh, start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
            synch_and_check();
        //      Tridiag
            jump_step = 1;
            grid= dim3(1, M - 1 , 1);
            tridiagSolver<<<grid, 32>>>(false, isU, start_idx, 
                            end_idx, jump_step, N + 3, 
                           d_arg_ptr, d_arr_ptr);
            synch_and_check();
 
        //     # Extract Solution
            FSj_extract_solution<<<grid_shape, block_shape>>>(ops.ketdinh, start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
            synch_and_check();
            Update_FS<<<grid_2d, block_2d>>>(d_arg_ptr);
        }
        


        // second half of the simulation

        t += dT * 0.5;


		update_boundary_at_t(M, N, t, channel, ops.total_time, d_arg_ptr, coeffs);        
        synch_and_check();

        

        block_shape = dim3(512, 1, 1) ;
        grid_shape = dim3((int) (ceil(M / 1024.0)), N, 1);
        start_idx = 2;
        end_idx = N;
        jump_step = 2;
        isU = false;
        if ((ops.channel) && (ops.kenhhepng)){
            start_idx = 3;
            end_idx = N;
        }


        VZSolver_calculate_preindex <<<grid_shape,block_shape>>> (start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
        synch_and_check();
            
        VZSolver_calculate_abcd<<<grid_shape,block_shape>>> (start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
        synch_and_check();
        VZSolver_calculate_matrix_coeff<<<grid_shape,block_shape>>> (start_idx, end_idx, NANGDAY, d_arg_ptr, d_arr_ptr);
        synch_and_check();

        grid = dim3(1, N-1, 1);
        tridiagSolver<<<grid, 32>>> (false, isU, start_idx, end_idx, jump_step, 2 * M + 1, d_arg_ptr, d_arr_ptr);
        synch_and_check();

        VZSolver_extract_solution<<<grid_shape,block_shape>>> (start_idx, end_idx, NANGDAY, d_arg_ptr, d_arr_ptr);
        synch_and_check();    



        Normalize<<<grid_2d, block_2d>>> (isU, d_arg_ptr, d_arr_ptr, coeffs);
        synch_and_check();
        update_buffer<<<grid_2d, block_2d>>> (isU, d_arg_ptr, d_arr_ptr);
        synch_and_check();



        solveU<<<grid_2d, block_2d>>> (t, 2, M, d_arg_ptr, coeffs);
        synch_and_check();
        grid = dim3(1, M, 1);
        update_margin_elem_U<<<grid,32>>> (2, M, NANGDAY, d_arg_ptr);
        synch_and_check();



        // similar to first haft, isU here is true, since it normalize u value after solving for u
        Normalize<<<grid_2d, block_2d>>> (true, d_arg_ptr, d_arr_ptr, coeffs);
		synch_and_check();
		update_buffer <<<grid_2d, block_2d>>>(true, d_arg_ptr, d_arr_ptr );
		synch_and_check();

		
		update_h_moi <<<grid_2d, block_2d>>> (d_arg_ptr);
        synch_and_check();

        grid = dim3(M, 1, 1);
        dim3 block(1, 32, 1);
        Reset_states_vertical <<<grid, block>>> (d_arg_ptr, coeffs);
        synch_and_check();

        if ( ((int) t % ops.interval == 0) && (t - (int) t == 0))
        	save_result(h_arg_pointer, (int) t);


        update_uvz <<<grid_2d, block_2d>>> (d_arg_ptr, coeffs);
        synch_and_check();
   

        grid = dim3(1, N, 1);
        Find_Calculation_limits_Horizontal <<<grid, 32>>> (d_arg_ptr, coeffs);
        grid = dim3(1, M, 1);
        Find_Calculation_limits_Vertical <<<grid, 32>>>(d_arg_ptr, coeffs);
        Htuongdoi <<<grid_2d, block_2d>>> (d_arg_ptr);
        synch_and_check();

        if (t >= ops.sediment_start){            
            Find_VTH<<<grid_shape, block_shape>>>(coeffs, d_arg_ptr);
            hesoK<<<grid_shape, block_shape>>>(coeffs, d_arg_ptr);
            synch_and_check();

        //     # sediment kernels come here
        //     # Scan FSj
            start_idx = 3;
            end_idx = N - 1;
            Scan_FSi<<<grid_shape, block_shape>>>(t, ops.bed_change_start, ops.ketdinh, start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
            synch_and_check();

            jump_step = 1;
        //     # Tridiag
            grid = dim3(1, N-3, 1);
            tridiagSolver<<<grid, 32>>>(false, isU, start_idx,
                            end_idx, jump_step, M + 3, 
                            d_arg_ptr, d_arr_ptr);
            synch_and_check();

        //     # Extract Solution
            FSi_extract_solution<<<grid_shape, block_shape>>>(ops.ketdinh, start_idx, end_idx, d_arg_ptr, d_arr_ptr, coeffs);
            synch_and_check();
            Update_FS<<<grid_2d, block_2d>>>(d_arg_ptr);
            synch_and_check();
            if ( ((int) t % 360 == 0) && (t - (int) t == 0)){
                // note: ask prof Bay about ros coefficient, ros here should be the same with the one in constant coeffs
                DOUBLE ros = 2000;
                save_FS(h_arg_pointer, (int) t, ros);
            }
           }

        //  only here is new. Verify if this work 
        // if int(t) % bed_change_start == 0 and t - int(t) == 0:
        //     Calculate_Qb(ketdinh, arg_struct_ptr, arr_struct_ptr, block=block_2d, grid=grid_2d)
        //     BedLoad(floattype(t), ketdinh, start_idx, end_idx, arg_struct_ptr, arr_struct_ptr, block=block_2d, grid=grid_2d)

	}
}




