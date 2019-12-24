
#include "engine.h"
#include "UVZSolver_multithread.h"
#include "support_funcs.h"
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
		dim3 block(1, 1024, 1);
		Update_Boundary_Value<<<grid, block>>>(t, total_time, d_arg_ptr);
	}
}

template <typename T>
void save_file(T* array, int width, int height, const char* filename){
	ofstream ofs;
	ofs.open(filename);
	if (!ofs)
		cout << "cannot open file " << filename << endl;
	for (int i = 1; i <= height; i++){
		for (int j = 1; j <= width; j++)
			ofs << array[i * (width + 3) + j] << " ";
		ofs << endl;
	}
}



template void save_file<double>(double*, int, int, const char*);
template void save_file<float>(float*, int, int, const char*);
template void save_file<int>(int*, int, int, const char*);



void save_result(Argument_Pointers h_arg_pointer, int t, bool save_FS=false){
	DOUBLE *u, *v, *z;

	int size = sizeof(DOUBLE) * (h_arg_pointer.M + 3) * (h_arg_pointer.N + 3);
	u = (DOUBLE*) malloc(size);
	cudaError_t status = cudaMemcpy((void*) u, h_arg_pointer.t_u, size, cudaMemcpyDeviceToHost);
	assert(status == cudaSuccess);
	v = (DOUBLE*) malloc(size);
	status = cudaMemcpy((void*) v, h_arg_pointer.t_v, size, cudaMemcpyDeviceToHost);
	assert(status == cudaSuccess);
	z = (DOUBLE*) malloc(size);
	status = cudaMemcpy((void*) z, h_arg_pointer.t_z, size, cudaMemcpyDeviceToHost);
	assert(status == cudaSuccess);
	save_file <double> (u, h_arg_pointer.M, h_arg_pointer.N, ("Outputs/u_" + to_string(t) + ".txt").c_str());
	save_file <double> (v, h_arg_pointer.M, h_arg_pointer.N, ("Outputs/v_" + to_string(t) + ".txt").c_str());
	save_file <double> (z, h_arg_pointer.M, h_arg_pointer.N, ("Outputs/z_" + to_string(t) + ".txt").c_str());


}


void Hydraulic_Calculation(DOUBLE dT, DOUBLE NANGDAY, Argument_Pointers* d_arg_ptr, Array_Pointers* d_arr_ptr, 
						Constant_Coeffs* coeffs, Argument_Pointers h_arg_pointer, Options ops){
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
	int Tmax = ops.Tmax;
	// DOUBLE Tmax = 0.5;
	cout << "t = " << t << endl;
	while (t < Tmax){
		t += 0.5 * dT;
		cout << t << endl;

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
		tridiagSolver<<<grid, 32>>> (false,isU, start_idx, end_idx, jump_step, 2 * N + 1, d_arg_ptr, d_arr_ptr);
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
        if ( ((int) t % ops.interval == 0) && (t - (int) t == 0))
        	save_result(h_arg_pointer, (int) t);  


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
        synch_and_check();;


        update_uvz <<<grid_2d, block_2d>>> (d_arg_ptr, coeffs);
        synch_and_check();
   

        grid = dim3(1, N, 1);
        Find_Calculation_limits_Horizontal <<<grid, 32>>> (d_arg_ptr, coeffs);
        grid = dim3(1, M, 1);
        Find_Calculation_limits_Vertical <<<grid, 32>>>(d_arg_ptr, coeffs);
        Htuongdoi <<<grid_2d, block_2d>>> (d_arg_ptr);
        synch_and_check();

        

       // sediment transport simulation here


	}
}




