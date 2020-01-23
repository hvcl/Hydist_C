/*******************************************************************************************************
                              University of Illinois/NCSA Open Source License
                                 Copyright (c) 2012 University of Illinois
                                          All rights reserved.

                                        Developed by: IMPACT Group
                                          University of Illinois
                                      http://impact.crhc.illinois.edu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
to deal with the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

  Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimers in the documentation and/or other materials provided with the distribution.
  Neither the names of IMPACT Group, University of Illinois, nor the names of its contributors may be used to endorse or promote products derived from this Software without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.

*******************************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "spike_kernel.hxx"
#include <assert.h>

__device__
void findBestGrid( int m, int tile_marshal, int *p_m_pad, int *p_b_dim, int *p_s, int *p_stride)
{
    int b_dim, m_pad, s, stride;
    int B_DIM_MAX, S_MAX;
    
    if ( sizeof(DOUBLE) == 4) {
        B_DIM_MAX = 256;
        S_MAX     = 512;    
    }
    else if (sizeof(DOUBLE) == 8){ /* double and complex */
        B_DIM_MAX = 128;
        S_MAX     = 256;     
    }
    else { /* doubleComplex */
        B_DIM_MAX = 64;
        S_MAX      = 128;    
    }
    
    /* b_dim must be multiple of 32 */
    if ( m < B_DIM_MAX * tile_marshal ) {
        b_dim = max( 32, (m/(32*tile_marshal))*32);
        s = 1;
        m_pad = ((m + b_dim * tile_marshal -1)/(b_dim * tile_marshal)) * (b_dim * tile_marshal);
        stride = m_pad/(s*b_dim);    
    }
    else {
        b_dim = B_DIM_MAX;
        
        s = 1;
        do {       
            int s_tmp = s * 2;
            int m_pad_tmp = ((m + s_tmp*b_dim*tile_marshal -1)/(s_tmp*b_dim*tile_marshal)) * (s_tmp*b_dim*tile_marshal);           
            float diff = (float)(m_pad_tmp - m)/float(m);
            /* We do not want to have more than 20% oversize */
            if ( diff < .2 ) {
                s = s_tmp;      
            }
            else {
                break;
            }
        } while (s < S_MAX);
                       
        m_pad = ((m + s*b_dim*tile_marshal -1)/(s*b_dim*tile_marshal)) * (s*b_dim*tile_marshal);        
        stride = m_pad/(s*b_dim);
    }
      
    *p_stride = stride;
    *p_m_pad  = m_pad;
    *p_s      = s;
    *p_b_dim  = b_dim;        
}


__device__
void gtsv_spike_partial_diag_pivot_v1(const DOUBLE* dl, const DOUBLE* d, const DOUBLE* du, DOUBLE* b, DOUBLE* x, const int m)
{


    // cudaFuncSetCacheConfig(tiled_diag_pivot_x1<DOUBLE,DOUBLE>,cudaFuncCachePreferL1);
    // cudaFuncSetCacheConfig(spike_GPU_back_sub_x1<DOUBLE>,cudaFuncCachePreferL1);
    //parameter declaration
    int s; //griddim.x
    int stride;
    int b_dim, m_pad;
    int tile = 128;
    int tile_marshal = 16;
    int T_size = sizeof(DOUBLE);
    
    findBestGrid( m, tile_marshal, &m_pad, &b_dim, &s, &stride);
   
    printf("m=%d m_pad=%d s=%d b_dim=%d stride=%d\n", m, m_pad, s, b_dim, stride);    
            
    int local_reduction_share_size = 2*b_dim*3*T_size;
    int global_share_size = 2*s*3*T_size;
    int local_solving_share_size = (2*b_dim*2+2*b_dim+2)*T_size;
    int marshaling_share_size = tile_marshal*(tile_marshal+1)*T_size;
    
    
    dim3 g_data(b_dim/tile_marshal,s);
    dim3 b_data(tile_marshal,tile_marshal);
    
    bool* flag; // tag for pivoting
    DOUBLE* dl_buffer;   //dl buffer
    DOUBLE* d_buffer;    //b
    DOUBLE* du_buffer; 
    DOUBLE* b_buffer;
    DOUBLE* w_buffer;
    DOUBLE* v_buffer;
    DOUBLE* c2_buffer;
    
    DOUBLE* x_level_2;
    DOUBLE* w_level_2;
    DOUBLE* v_level_2;
    
    
    //buffer allocation
    cudaError_t status =  cudaMalloc((void **)&flag, sizeof(bool)*m_pad); 
    assert(status == cudaSuccess);
    status = cudaMalloc((void **)&dl_buffer, T_size*m_pad); 
    assert(status == cudaSuccess);
    status = cudaMalloc((void **)&d_buffer, T_size*m_pad); 
    assert(status == cudaSuccess);
    status = cudaMalloc((void **)&du_buffer, T_size*m_pad); 
    assert(status == cudaSuccess);
    status = cudaMalloc((void **)&b_buffer, T_size*m_pad); 
    assert(status == cudaSuccess);
    status = cudaMalloc((void **)&w_buffer, T_size*m_pad); 
    assert(status == cudaSuccess);
    status = cudaMalloc((void **)&v_buffer, T_size*m_pad); 
    assert(status == cudaSuccess);
    status = cudaMalloc((void **)&c2_buffer, T_size*m_pad); 
    assert(status == cudaSuccess);
    
    cudaMalloc((void **)&x_level_2, T_size*s*2);
    assert(status == cudaSuccess);
    cudaMalloc((void **)&w_level_2, T_size*s*2);
    assert(status == cudaSuccess); 
    cudaMalloc((void **)&v_level_2, T_size*s*2);
    assert(status == cudaSuccess); 
    
    
    
    //kernels 
    
    //data layout transformation
    foward_marshaling_bxb<DOUBLE><<<g_data ,b_data, marshaling_share_size >>>(dl_buffer, dl, stride, b_dim,m, cuGet(0));
    foward_marshaling_bxb<DOUBLE><<<g_data ,b_data, marshaling_share_size >>>(d_buffer,  d,  stride, b_dim,m, cuGet(1));
    foward_marshaling_bxb<DOUBLE><<<g_data ,b_data, marshaling_share_size >>>(du_buffer, du, stride, b_dim,m, cuGet(0));
    foward_marshaling_bxb<DOUBLE><<<g_data ,b_data, marshaling_share_size >>>(b_buffer,  b,  stride, b_dim,m, cuGet(0));
     
    //partitioned solver
    //tiled_diagonal_pivoting<<<s,b_dim>>>( x,w,v,c2_buffer,flag, dl,d,du,b, stride,tile);
    tiled_diag_pivot_x1<DOUBLE,DOUBLE><<<s,b_dim>>>(b_buffer, w_buffer, v_buffer, c2_buffer, flag, dl_buffer, d_buffer, du_buffer, stride, tile);
    
    
    //SPIKE solver
    spike_local_reduction_x1<DOUBLE><<<s,b_dim,local_reduction_share_size>>>(b_buffer,w_buffer,v_buffer,x_level_2, w_level_2, v_level_2,stride);
    spike_GPU_global_solving_x1<<<1,32,global_share_size>>>(x_level_2,w_level_2,v_level_2,s);
    spike_GPU_local_solving_x1<DOUBLE><<<s,b_dim,local_solving_share_size>>>(b_buffer,w_buffer,v_buffer,x_level_2,stride);
    spike_GPU_back_sub_x1<DOUBLE><<<s,b_dim>>>(b_buffer,w_buffer,v_buffer, x_level_2,stride);

    back_marshaling_bxb<DOUBLE><<<g_data ,b_data, marshaling_share_size >>>(x,b_buffer,stride,b_dim,m);
    // back_marshaling_bxb<DOUBLE><<<g_data ,b_data, marshaling_share_size >>>(b,b_buffer,stride,b_dim,m);
    
    //free
    
    cudaFree(flag);
    cudaFree(dl_buffer);
    cudaFree(d_buffer);
    cudaFree(du_buffer);
    cudaFree(b_buffer);
    cudaFree(w_buffer);
    cudaFree(v_buffer);
    cudaFree(c2_buffer);
    cudaFree(x_level_2);
    cudaFree(w_level_2);
    cudaFree(v_level_2);    
    printf("tId %d finished successfully \n",threadIdx.x);          
    
}
