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

#include "ops_device.hxx"


//Data layout transformation for inputs
template <typename T> 
__global__ void foward_marshaling_bxb ( T* x,
                                        const T* y,
                                        const int h_stride,
                                        const int l_stride,
                                        int m,
                                        T pad
                                        )
{	
	int b_dim;


	int global_in;
	int global_out;
	int shared_in;
	int shared_out;
	
	
	b_dim = blockDim.x; 	//16

	global_in = blockIdx.y*l_stride*h_stride +  (blockIdx.x*b_dim+threadIdx.y)*h_stride+threadIdx.x;
	global_out = blockIdx.y*l_stride*h_stride + threadIdx.y*l_stride + blockIdx.x*b_dim+threadIdx.x;
	shared_in = threadIdx.y*(b_dim+1)+threadIdx.x;
	shared_out = threadIdx.x*(b_dim+1)+threadIdx.y;

	int k;

    struct __dynamic_shmem__<T> shmem; 
    T *share = shmem.getPtr();

	for(k=0;k<h_stride;k+=b_dim)
	{	
		share[shared_in]= global_in >= m ? pad : y[global_in];
		global_in += b_dim;
		
		__syncthreads();
		
		x[global_out] = share[shared_out];
		global_out+=b_dim*l_stride;
		__syncthreads();
	}		
}

//Data layout transformation for results
template <typename T> 
__global__ void  back_marshaling_bxb (
                                      T* x,
                                      const T* y,
                                      const int h_stride,
                                      const int l_stride,
                                      int m
                                      )
{	
	int b_dim;
	
	int global_in;
	int global_out;
	int shared_in;
	int shared_out;
	
	b_dim = blockDim.x; 	//16

	global_out = blockIdx.y*l_stride*h_stride +  (blockIdx.x*b_dim+threadIdx.y)*h_stride+threadIdx.x;
	global_in = blockIdx.y*l_stride*h_stride + threadIdx.y*l_stride + blockIdx.x*b_dim+threadIdx.x;
	shared_in = threadIdx.y*(b_dim+1)+threadIdx.x;
	shared_out = threadIdx.x*(b_dim+1)+threadIdx.y;
	
	int k;

    struct __dynamic_shmem__<T> shmem; 
    T *share = shmem.getPtr();

	for(k=0;k<h_stride;k+=b_dim)
	{
	
		share[shared_in]=y[global_in];
		global_in += b_dim*l_stride;
		
		__syncthreads();
		
        if (global_out < m) {
		    x[global_out] = share[shared_out];
        }
		global_out+=b_dim;
		__syncthreads();
	}		
}

//Partitioned solver with tiled diagonal pivoting
// T_ELEM_REAL to hold the type of sgema (it is necesary for complex variants
template <typename T_ELEM , typename T_ELEM_REAL> 
__global__ void tiled_diag_pivot_x1(
                                      T_ELEM* x,
                                      T_ELEM* w,  //left halo
                                      T_ELEM* v,  //right halo
                                      T_ELEM* b_buffer,  //modified main diag
                                      bool *flag,  //buffer to tag pivot
                                      const T_ELEM* a,    //lower diag
                                      const T_ELEM* b,    //main diag
                                      const T_ELEM* c,    //upper diag
                                      const int stride,
                                      const int tile
                                      )                                    
{
	
	int b_dim;
	int ix;
	int bx;
	
	bx = blockIdx.x;
	b_dim = blockDim.x;
	ix = bx*stride*b_dim+threadIdx.x;

	
	int k=0;
	T_ELEM b_k,b_k_1,a_k_1,c_k,c_k_1,a_k_2;
	T_ELEM x_k,x_k_1;
	T_ELEM w_k,w_k_1;
	T_ELEM v_k_1;
	
	T_ELEM_REAL kia = (sqrt(5.0)-1.0)/2.0;
	b_k = b[ix];
	c_k = c[ix];
	//x_k = d[ix];
    x_k = x[ix];
	w_k = a[ix];
	
	a_k_1 = a[ix+b_dim];
	b_k_1 = b[ix+b_dim];
	c_k_1 = c[ix+b_dim];
	//x_k_1 = d[ix+b_dim];
    x_k_1 = x[ix+b_dim];
	
	a_k_2 = a[ix+2*b_dim];
		
	int i;
		
	//forward
	for(i=1;i<=tile;i++)
	{
		while(k<(stride*i)/tile)
		{        
			T_ELEM_REAL sgema;
			
			// math.h has an intrinsics for float, double 
            sgema = max(cuAbs(c_k), cuAbs(a_k_1));
			sgema = max(sgema, cuAbs(b_k_1));
			sgema = max( sgema, cuAbs(c_k_1));
			sgema = max( sgema, cuAbs(a_k_2));			
			
			if( cuMul(cuAbs(b_k),sgema) >= cuMul( kia, cuMul(cuAbs(c_k), cuAbs(a_k_1)) ))
			{    
                T_ELEM b_inv = cuDiv( cuGet(1), b_k);
				//write back
				flag[ix]=true;
				
				x_k = cuMul( x_k, b_inv );
				w_k = cuMul( w_k, b_inv );
				
				x[ix] = x_k;		//k
				w[ix] = w_k;		//k
				b_buffer[ix]=b_k;
				//

				if( k < stride-1)
				{
					ix+=b_dim;
					//update					        
                    x_k = cuFma( cuNeg(a_k_1) , x_k, x_k_1); //k+1
					w_k = cuMul( cuNeg(a_k_1), w_k);        //k+1                    					
                    b_k = cuFma( cuNeg(a_k_1), cuMul(c_k,b_inv) , b_k_1);         //k+1
					
					if( k < stride-2)				
					{
						//load new data
						b_k_1 = b[ix+b_dim];  //k+2
						a_k_1 = a_k_2;		  //k+2
						//x_k_1 = d[ix+b_dim];
                        x_k_1 = x[ix+b_dim];
						c_k   = c_k_1;		  //k+1
						c_k_1 = c[ix+b_dim];  //k+2
						
						a_k_2 = k< (stride-3) ? a[ix+2*b_dim] : cuGet(0); //k+3
						
					}
					else			//k =stride -2
					{
						b_k_1 = cuGet(0);
						a_k_1 = cuGet(0);
						x_k_1 = cuGet(0);
						c_k   = cuGet(0);
						c_k_1 = cuGet(0);
						a_k_2 = cuGet(0);
					}
				}
				else		//k=stride -1
				{
					v[ix] = cuMul( c[ix], b_inv );
					ix   += b_dim;
				}
				
				k+=1;             							
			}
			else
			{		
				T_ELEM delta;
								
                delta = cuFma( b_k, b_k_1, cuNeg(cuMul(c_k,a_k_1)) );
				delta = cuDiv( cuGet(1) , delta );				                
                x[ix] = cuFma( x_k, b_k_1, cuNeg( cuMul(c_k,x_k_1))); //k
                x[ix] = cuMul( x[ix], delta);
                
				w[ix]  = cuMul(w_k, cuMul(b_k_1,delta)); //k
				b_buffer[ix]=b_k;				
				flag[ix] = false;
				
                x_k_1 = cuFma( b_k,x_k_1, cuNeg(cuMul(a_k_1,x_k))); //k+1
                x_k_1 = cuMul(x_k_1,delta);
				w_k_1 = cuMul(cuMul(cuNeg(a_k_1),w_k), delta);	  //k+1
				
				x[ix+b_dim]        = x_k_1;	  //k+1
				w[ix+b_dim]        = w_k_1;	  //k+1
				b_buffer[ix+b_dim] = b_k_1;				
				flag[ix+b_dim]=false;	
				
				if(k<stride-2)
				{
					ix+=2*b_dim;		
					//update					
                    //x_k = cuFma(cuNeg(a_k_2), x_k_1, d[ix]);     //k+2      
                    x_k = cuFma(cuNeg(a_k_2), x_k_1, x[ix]);       //k+2              
					w_k = cuMul(cuNeg(a_k_2),w_k_1);     //k+2					
                    b_k = cuMul(cuMul(a_k_2,b_k), cuMul(c_k_1,delta));		//k+2
                    b_k = cuSub( b[ix],b_k);
					
					if(k<stride-3)
					{
						//load new data
						c_k = c[ix];     //k+2
						b_k_1 = b[ix+b_dim];  //k+3
						a_k_1 = a[ix+b_dim];  //k+3
						c_k_1 = c[ix+b_dim];  //k_3
						//x_k_1 = d[ix+b_dim];  //k+3
                        x_k_1 = x[ix+b_dim];  //k+3
						a_k_2 = k<stride-4? a[ix+2*b_dim] : cuGet(0);
					}
					else		//k=stride-3
					{
						b_k_1 = cuGet(0);
						a_k_1 = cuGet(0);
						c_k   = cuGet(0);
						c_k_1 = cuGet(0);
						x_k_1 = cuGet(0);
						a_k_2 = cuGet(0);
					}
				}
				else		////k=stride -2
				{
					
					T_ELEM v_temp;
					v_temp = cuMul( c[ix+b_dim], delta);					
                    v[ix]  = cuMul(v_temp, cuNeg(c_k));
					v[ix+b_dim] = cuMul(v_temp,b_k);
					ix+= 2*b_dim;
				}				
				k+=2;	                			
			}
			
		}
	}


	//k=stride
	
	//backward
	//last one
	k-=1;
	ix-=b_dim;
	if(flag[ix])
	{
		x_k_1=x[ix];
		w_k_1=w[ix];
		v_k_1=v[ix];
		k-=1;
		//x[ix]
	}
	else		//2-by-2
	{
		ix-=b_dim;
		x_k_1=x[ix];
		w_k_1=w[ix];
		v_k_1=v[ix];
		k-=2;
	}
	ix-=b_dim;
	
	for(i=tile-1;i>=0;i--)	{
		while(k>=(i*stride)/tile){
			if(flag[ix]){		//1-by-1			
				c_k = c[ix];
				b_k = b_buffer[ix];				
                
                T_ELEM tempDiv = cuDiv(cuNeg(c_k),b_k);
                x_k_1 = cuFma( x_k_1, tempDiv, x[ix]);                				
                w_k_1 = cuFma( w_k_1, tempDiv , w[ix]);                
				v_k_1 = cuMul(v_k_1, tempDiv);
                
				x[ix] = x_k_1;
				w[ix] = w_k_1;
				v[ix] = v_k_1;
				k-=1;
			}
			else {
			
				T_ELEM delta;
				b_k   = b_buffer[ix-b_dim];
				c_k   = c[ix-b_dim];
				a_k_1 = a[ix];
				b_k_1 = b_buffer[ix];
				c_k_1 = c[ix];
                delta = cuFma( b_k, b_k_1, cuNeg(cuMul(c_k,a_k_1)) );
				delta = cuDiv( cuGet(1) , delta );	                
                
                T_ELEM prod = cuMul(c_k_1 , cuMul(b_k , delta));
                
				x[ix] =  cuFma(cuNeg(x_k_1), prod, x[ix]);
				w[ix] =  cuFma(cuNeg(w_k_1), prod, w[ix]);
				v[ix] =  cuMul(cuNeg(v_k_1),prod);
				
                ix  -= b_dim;
                prod = cuMul(c_k_1 , cuMul(c_k , delta));
                
				x_k_1 = cuFma( x_k_1,prod, x[ix]);
				w_k_1 = cuFma( w_k_1,prod, w[ix]);
				v_k_1 = cuMul(v_k_1,prod);
				x[ix] = x_k_1;
				w[ix] = w_k_1;
				v[ix] = v_k_1;
				k-=2;
			}
			ix-=b_dim;
		}		
	}	
}

//
// first version of x32 and long x32
//
template <typename T_ELEM , typename T_ELEM_REAL> 
__global__ void tiled_diag_pivot_x32(
                                      T_ELEM* x,
                                      bool *flag,  //buffer to tag pivot
                                      const T_ELEM* a,    //lower diag
                                      const T_ELEM* b,    //modified one main diag****
                                      const T_ELEM* c,    //upper diag
                                      const int stride,
                                      const int tile
                                      )                                    
{
	
	int b_dim;
	int ix;
	b_dim = blockDim.x;
	ix = blockIdx.x*stride*b_dim+threadIdx.x;
	int iy;
	iy = blockIdx.y*stride;

	
	int k=0;
	T_ELEM b_k,b_k_1,a_k_1,c_k,c_k_1;
	T_ELEM x_k,x_k_1;
	
	T_ELEM_REAL kia = (sqrt(5.0)-1.0)/2.0;
	x_k = x[ix];
	b_k = b[iy];
	c_k = c[iy];
	
	
	a_k_1 = a[iy+1];
	b_k_1 = b[iy+1];
	c_k_1 = c[iy+1];
    x_k_1 = x[ix+b_dim];
	
			
	int i;
		
	//forward
	for(i=1;i<=tile;i++)
	{
		while(k<(stride*i)/tile)
		{        

			if( flag[iy])
			{    
                T_ELEM b_inv = cuDiv( cuGet(1), b_k);
				x_k = cuMul( x_k, b_inv );
				x[ix] = x_k;		//k
				
				ix+=b_dim;
				iy+=1;
				if( k < stride-1)
				{
					//update					        
                    x_k = cuFma( cuNeg(a_k_1) , x_k, x_k_1); //k+1
					
					if( k < stride-2)				
					{
						//load new data
						a_k_1 = a[iy+1];
	                    x_k_1 = x[ix+b_dim];
						
						
					}
					else			//k =stride -2
					{
						a_k_1 = cuGet(0);
						x_k_1 = cuGet(0);
						

					}
				}			
				k+=1;             							
			}
			else
			{		
				T_ELEM delta;
                delta = cuFma( b_k, b_k_1, cuNeg(cuMul(c_k,a_k_1)) );
				delta = cuDiv( cuGet(1) , delta );				                
                x[ix] = cuFma( x_k, b_k_1, cuNeg( cuMul(c_k,x_k_1))); //k
                x[ix] = cuMul( x[ix], delta);
                x_k_1 = cuFma( b_k,x_k_1, cuNeg(cuMul(a_k_1,x_k))); //k+1
                x_k_1 = cuMul(x_k_1,delta);
	
				x[ix+b_dim]        = x_k_1;	  //k+1
				ix+=2*b_dim;	
				
				
				if(k<stride-2)
				{
					//update					
					x_k = cuFma(cuNeg(a[iy+2]), x_k_1, x[ix]);       //k+2   
					iy+=2;
					if(k<stride-3)
					{
						//load new data
						a_k_1 = a[iy+1];  //k+3
                        x_k_1 = x[ix+b_dim];  //k+3
	
					}
					else		//k=stride-3
					{
						a_k_1 = cuGet(0);
						x_k_1 = cuGet(0);
					}
					     
				}
		
				k+=2;	                			
			}
			
		}
	}
	//k=stride
	//backward
	//last one
	k-=1;
	ix-=b_dim;
	iy-=1;
	if(flag[ix])
	{
		x_k_1=x[ix];
		k-=1;
	}
	else		//2-by-2
	{
		ix-=b_dim;
		iy-=1;
		x_k_1=x[ix];
		k-=2;
	}
	ix-=b_dim;
	iy-=1;
	
	for(i=tile-1;i>=0;i--)	{
		while(k>=(i*stride)/tile){
			if(flag[iy]){		//1-by-1			
				c_k = c[iy];
				b_k = b[iy];				
               
                T_ELEM tempDiv = cuDiv(cuNeg(c_k),b_k);
                x_k_1 = cuFma( x_k_1, tempDiv, x[ix]);                				
                
				x[ix] = x_k_1;
				k-=1;
			}
			else {
			
				T_ELEM delta;
				b_k   = b[iy-1];
				c_k   = c[iy-1];
				a_k_1 = a[iy];
				b_k_1 = b[iy];
				c_k_1 = c[iy];
                delta = cuFma( b_k, b_k_1, cuNeg(cuMul(c_k,a_k_1)) );
				delta = cuDiv( cuGet(1) , delta );	                
                
                T_ELEM prod = cuMul(c_k_1 , cuMul(b_k , delta));
                
				x[ix] =  cuFma(cuNeg(x_k_1), prod, x[ix]);
				
                ix  -= b_dim;
				iy  -=1;
                prod = cuMul(c_k_1 , cuMul(c_k , delta));
                
				x_k_1 = cuFma( x_k_1,prod, x[ix]);
				x[ix] = x_k_1;
				k-=2;
			}
			ix-=b_dim;
			iy-=1;
		}		
	}	
}

//SPIKE solver within a thread block for 1x rhs
template <typename T> 
__global__ void 
spike_local_reduction_x1
(
T* x,
T* w,
T* v,
T* x_mirror,
T* w_mirror,
T* v_mirror,
const int stride  //stride per thread
)
{
	int tx;
	int b_dim;
	int bx;
	
	tx = threadIdx.x;
	b_dim = blockDim.x;
	bx = blockIdx.x;
	//
	//extern __shared__ T shared[];
    struct __dynamic_shmem__<T> shmem; 
    T *shared = shmem.getPtr();
        
	T* sh_w=shared;				
	T* sh_v=sh_w+2*b_dim;				
	T* sh_x=sh_v+2*b_dim;			
	
	//a ~~ w
	//b ~~ I
	//c ~~ v
	//d ~~ x
	
	int base = bx*stride*b_dim;
	
	//load halo to scratchpad
	sh_w[tx] = w[base+tx];
	sh_w[tx+b_dim] = w[base+tx+(stride-1)*b_dim];
	sh_v[tx] = v[base+tx];
	sh_v[tx+b_dim] = v[base+tx+(stride-1)*b_dim];
	sh_x[tx] = x[base+tx];
	sh_x[tx+b_dim] = x[base+tx+(stride-1)*b_dim];
	
	__syncthreads();



	
	int scaler = 2;
	
	while(scaler<=b_dim)
	{
		if(tx < b_dim/scaler)
		{
			int index;
			int up_index;
			int down_index;
			index = scaler*tx+scaler/2-1;
			up_index= scaler*tx;
			down_index = scaler*tx + scaler-1;
			T det = cuGet(1);
			det = cuFma( cuNeg(sh_v[index+b_dim]), sh_w[index+1], det);
			det = cuDiv( cuGet(1) , det);
			
			T d1,d2;
			d1 = sh_x[index+b_dim];
			d2 = sh_x[index+1];
			
            sh_x[index+b_dim] = cuMul( cuFma( sh_v[index+b_dim], cuNeg(d2), d1), det);
            sh_x[index+1]     = cuMul( cuFma(sh_w[index+1],cuNeg(d1), d2), det);			
            sh_w[index+1] = cuMul( sh_w[index+b_dim], cuMul(sh_w[index+1], cuNeg(det)));	            
			sh_w[index+b_dim] = cuMul(sh_w[index+b_dim],det);
									
			sh_v[index+b_dim] = cuMul(sh_v[index+b_dim], cuMul(sh_v[index+1], cuNeg(det)));            
			sh_v[index+1] = cuMul(sh_v[index+1],det);
			
			
			//boundary
            sh_x[up_index] 		= cuFma( sh_x[index+1], cuNeg(sh_v[up_index]), sh_x[up_index]);            
			sh_x[down_index+b_dim] = cuFma(sh_x[index+b_dim], cuNeg(sh_w[down_index+b_dim]), sh_x[down_index+b_dim]);
            
            sh_w[up_index] = cuFma( sh_w[index+1], cuNeg(sh_v[up_index]), sh_w[up_index]);
			
			sh_v[up_index] 		= cuMul( cuNeg(sh_v[index+1]), sh_v[up_index]);

            sh_v[down_index+b_dim] = cuFma(sh_v[index+b_dim], cuNeg(sh_w[down_index+b_dim]), sh_v[down_index+b_dim]);
			
            sh_w[down_index+b_dim] 	= cuMul( cuNeg(sh_w[index+b_dim]), sh_w[down_index+b_dim]);
			
		}
		scaler*=2;
		__syncthreads();
	}
	
	//write out
	
/*
	scaler = b_dim/2;
	
	while(scaler>=2)
	{
		if(tx < b_dim/scaler)
		{
			int index;
			int up_index;
			int down_index;
			index = scaler*tx+scaler/2-1;
			up_index= scaler*tx-1;
			down_index = scaler*tx + scaler;
			//up_index=up_index<0?0:up_index;
		//	down_index=down_index<len?down_index:len-1;
		
			double up_value,down_value;
			up_value = up_index<0? 0.0: sh_x[up_index+b_dim];
			down_value = down_index<b_dim?sh_x[down_index]:0.0;
			sh_x[index+b_dim] -= sh_w[index+b_dim]*up_value+sh_v[index+b_dim]*down_value;
			sh_x[index+1]   -= sh_w[index+1]*up_value+sh_v[index+1]*down_value;
			
			double temp_1,temp_2;
			temp_1 = sh_w[index+b_dim];
			temp_2 = sh_w[index+1];		
			
			
			up_value = up_index<0? 0.0: sh_w[up_index+b_dim];
			down_value = down_index<b_dim?sh_w[down_index]:0.0;
			sh_w[index+b_dim] = -sh_w[index+b_dim] *up_value-sh_v[index+b_dim]*down_value;
			sh_w[index+1] = -sh_w[index+1] *up_value-sh_v[index+1]*down_value;
			
			up_value = up_index<0? 0.0: sh_v[up_index+b_dim];
			down_value = down_index<b_dim?sh_v[down_index]:0.0;
			sh_v[index+b_dim] = -temp_1 *up_value-sh_v[index+b_dim]*down_value;
			sh_v[index+1] = -temp_2 *up_value-sh_v[index+1]*down_value;
			
		}
		scaler/=2;
		__syncthreads();
	}
	*/
	
	w[base+tx] =sh_w[tx];
	w[base+tx+(stride-1)*b_dim] = sh_w[tx+b_dim];
	
	v[base+tx] =sh_v[tx];
	v[base+tx+(stride-1)*b_dim] = sh_v[tx+b_dim];
	
	x[base+tx] =sh_x[tx];
	x[base+tx+(stride-1)*b_dim] = sh_x[tx+b_dim];
	
	//write mirror
	if(tx<1)
	{
		int g_dim=gridDim.x;
		w_mirror[bx] = sh_w[0];
		w_mirror[g_dim+bx] = sh_w[2*b_dim-1];
		
		v_mirror[bx] = sh_v[0];
		v_mirror[g_dim+bx] = sh_v[2*b_dim-1];
		
		x_mirror[bx] = sh_x[0];
		x_mirror[g_dim+bx] = sh_x[2*b_dim-1];
	}

}

///////////////////////////
/// a global level SPIKE solver for oneGPU
/// One block version
///
////////////////////
template <typename T> 
__global__ void 
spike_GPU_global_solving_x1
(
T* x,
T* w,
T* v,
const int len
)
{

	int ix;
	int b_dim;
	
	b_dim = blockDim.x;
	//
	//extern __shared__ T shared[];
    struct __dynamic_shmem__<T> shmem; 
    T *shared = shmem.getPtr();    
    
    
	T* sh_w=shared;				
	T* sh_v=sh_w+2*len;				
	T* sh_x=sh_v+2*len;	

	
	//a ~~ w
	//b ~~ I
	//c ~~ v
	//d ~~ x
	
	//read data
	ix = threadIdx.x;
	while(ix<len)
	{
		sh_w[ix] = w[ix];
		sh_w[ix+len] = w[ix+len];
		
		sh_v[ix] = v[ix];
		sh_v[ix+len] = v[ix+len];
		
		sh_x[ix] = x[ix];
		sh_x[ix+len] = x[ix+len];
		
		ix+=b_dim;
	}
	__syncthreads();
	
	
	
	int scaler = 2;
	while(scaler<=len)
	{
		ix = threadIdx.x;
		while(ix < len/scaler)
		{
			int index;
			int up_index;
			int down_index;
			index = scaler*ix+scaler/2-1;
			up_index= scaler*ix;
			down_index = scaler*ix + scaler-1;
			T det = cuGet(1);
			det = cuFma( cuNeg(sh_v[index+len]), sh_w[index+1], det);
			det = cuDiv( cuGet(1) , det);			
            
            
			T d1,d2;
			d1 = sh_x[index+len];
			d2 = sh_x[index+1];
			
            sh_x[index+len] = cuMul( cuFma( sh_v[index+len], cuNeg(d2), d1), det);
            
            sh_x[index+1]     = cuMul( cuFma(sh_w[index+1],cuNeg(d1), d2), det);
						
            sh_w[index+1] = cuMul( sh_w[index+len], cuMul(sh_w[index+1], cuNeg(det)));		
            
			sh_w[index+len] = cuMul(sh_w[index+len],det);
									
            sh_v[index+len] = cuMul(sh_v[index+len], cuMul(sh_v[index+1], cuNeg(det)));    
            
			sh_v[index+1] = cuMul(sh_v[index+1],det);
			
			
			//boundary
            sh_x[up_index] 		= cuFma( sh_x[index+1], cuNeg(sh_v[up_index]), sh_x[up_index]); 
            
            sh_x[down_index+len] = cuFma(sh_x[index+len], cuNeg(sh_w[down_index+len]), sh_x[down_index+len]);
						
            sh_w[up_index] = cuFma( sh_w[index+1], cuNeg(sh_v[up_index]), sh_w[up_index]);

            sh_v[up_index] 		= cuMul( cuNeg(sh_v[index+1]), sh_v[up_index]);	
            
            sh_v[down_index+len] = cuFma(sh_v[index+len], cuNeg(sh_w[down_index+len]), sh_v[down_index+len]);
            
            sh_w[down_index+len] 	= cuMul( cuNeg(sh_w[index+len]), sh_w[down_index+len]);
			
			ix+=b_dim;
			
		}
		scaler*=2;
		__syncthreads();
	}
	
	//backward reduction
	
	scaler = len/2;
	while(scaler>=2)
	{
		ix = threadIdx.x;
		while(ix < len/scaler)
		{
			int index;
			int up_index;
			int down_index;
			index = scaler*ix+scaler/2-1;
			up_index= scaler*ix-1;
			down_index = scaler*ix + scaler;
			up_index=up_index<0?0:up_index;
			down_index=down_index<len?down_index:len-1;
			
            sh_x[index+len] = cuFma( cuNeg(sh_w[index+len]), sh_x[up_index+len], sh_x[index+len]);
            sh_x[index+len] = cuFma( cuNeg(sh_v[index+len]), sh_x[down_index], sh_x[index+len]);
            
			sh_x[index+1]   = cuFma( cuNeg(sh_w[index+1]), sh_x[up_index+len], sh_x[index+1]);
            sh_x[index+1]   = cuFma( cuNeg(sh_v[index+1]), sh_x[down_index], sh_x[index+1]);	
            
			ix+=b_dim;
		}
		scaler/=2;
		__syncthreads();
	}
	
	//write out
	
	ix = threadIdx.x;
	while(ix<len)
	{
	
		x[ix] = sh_x[ix];
		x[ix+len] = sh_x[ix+len];
		ix+=b_dim;
	}
	
}



/// a thread-block level SPIKE solver 
template <typename T> 
__global__ void spike_GPU_local_solving_x1(
                                            T* x,
                                            const T* w,
                                            const T* v,
                                            const T* x_mirror,
                                            const int stride
                                          )
{
	int tx;
	int b_dim;
	int bx;
	
	tx = threadIdx.x;
	b_dim = blockDim.x;
	bx = blockIdx.x;
	//	
    struct __dynamic_shmem__<T> shmem; 
    T *shared = shmem.getPtr();
        
	T *sh_w=shared;				
	T *sh_v=sh_w+2*b_dim;				
	T *sh_x=sh_v+2*b_dim;		//sh_x is 2*b_dim + 2
	
	//a ~~ w
	//b ~~ I
	//c ~~ v
	//d ~~ x
	
	int base = bx*stride*b_dim;
	
	//load halo to scratchpad
	sh_w[tx] = w[base+tx];
	sh_w[tx+b_dim] = w[base+tx+(stride-1)*b_dim];
	sh_v[tx] = v[base+tx];
	sh_v[tx+b_dim] = v[base+tx+(stride-1)*b_dim];
	
	//swap the order of x
	//why
	sh_x[tx+1] = x[base+tx+(stride-1)*b_dim];
	sh_x[tx+b_dim+1] = x[base+tx];
	
	__syncthreads();
	
	if(tx<1)
	{
		int g_dim=gridDim.x;
		sh_x[0]=bx>0?x_mirror[bx-1+g_dim]: cuGet(0);
		
		sh_x[2*b_dim+1]=bx<g_dim-1?x_mirror[bx+1]: cuGet(0);
		
		sh_x[b_dim+1]=x_mirror[bx];
		sh_x[b_dim]=x_mirror[bx+g_dim];
		//sh_x[1]=x_mirror[bx];
		//sh_x[2*b_dim] = x_mirror[bx+g_dim];
		/*
		sh_w[0]=0.0;
		sh_w[2*b_dim-1]=0.0;
		sh_v[0]=0.0;
		sh_v[2*b_dim-1]=0.0;
		/**/
		
	}
	__syncthreads();
	
	int scaler = b_dim;
	while(scaler>=2)
	{
		if(tx < b_dim/scaler)
		{
			int index;
			int up_index;
			int down_index;
			index = scaler*tx+scaler/2-1;
			up_index= scaler*tx;
			down_index = scaler*tx + scaler+1;
	
			
		//	sh_x[index+b_dim+1] -= sh_w[index+b_dim]*sh_x[up_index]+sh_v[index+b_dim]*sh_x[down_index+b_dim];
		//	sh_x[index+2]   -= sh_w[index+1]*sh_x[up_index]+sh_v[index+1]*sh_x[down_index+b_dim];
			
            sh_x[index+1] = cuFma( sh_w[index+b_dim], cuNeg(sh_x[up_index]) , sh_x[index+1]);
            sh_x[index+1] = cuFma( sh_v[index+b_dim], cuNeg(sh_x[down_index+b_dim]) , sh_x[index+1]);
            
            sh_x[index+b_dim+2]   = cuFma( sh_w[index+1], cuNeg(sh_x[up_index]), sh_x[index+b_dim+2]);
            sh_x[index+b_dim+2]  = cuFma( sh_v[index+1],cuNeg(sh_x[down_index+b_dim]), sh_x[index+b_dim+2]);
						
		//	sh_x[index+1] = sh_w[index+b_dim];
		//	sh_x[index+b_dim+2]   = sh_w[index+1];
		}
		scaler/=2;
		__syncthreads();
	}
	
	
	
	//write out
	x[base+tx] =sh_x[tx+b_dim+1];
	
/*	
	int k;

	//x[base+tx] = sh_x[tx+b_dim+1];
	
	for(k=1;k<stride-1;k++)
	{
		x[base+tx+k*b_dim] -= w[base+tx+k*b_dim]*sh_x[tx] + v[base+tx+k*b_dim]*sh_x[tx+b_dim+2];
		//x[base+tx+k*b_dim] =sh_x[tx];
	}
	*/
	
	x[base+tx+(stride-1)*b_dim] = sh_x[tx+1];


}


//backward substitution for SPIKE solver
template <typename T> 
__global__ void 
spike_GPU_back_sub_x1
(
T* x,
const T* w,
const T* v,
const T* x_mirror,
const int stride
)
{
	int tx;
	int b_dim;
	int bx;

	
	tx = threadIdx.x;
	b_dim = blockDim.x;
	bx = blockIdx.x;
	int base = bx*stride*b_dim;
	T x_up,x_down;
	

	if(tx>0 && tx<b_dim-1)
	{
		x_up =x[base+tx-1+(stride-1)*b_dim];
		x_down = x[base+tx+1];
	}
	else
	{
		int g_dim=gridDim.x;
		if(tx==0)
		{
			x_up   = bx>0 ? x_mirror[bx-1+g_dim]:cuGet(0);
			x_down = x[base+tx+1];
		}
		else
		{
			x_up   = x[base+tx-1+(stride-1)*b_dim];
			x_down = bx<g_dim-1?x_mirror[bx+1]:cuGet(0);
		}
	}
	
	//x[base+tx] = sh_x[tx+b_dim+1];
	int k;
	for(k=1;k<stride-1;k++)
	{        
        x[base+tx+k*b_dim] = cuFma( w[base+tx+k*b_dim], cuNeg(x_up), x[base+tx+k*b_dim]);
        x[base+tx+k*b_dim] = cuFma( v[base+tx+k*b_dim], cuNeg(x_down), x[base+tx+k*b_dim]);
        
		//x[base+tx+k*b_dim] =sh_x[tx];
	}
}






