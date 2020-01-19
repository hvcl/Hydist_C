/* Copyright (c) 2009-2013,  NVIDIA CORPORATION

   All rights reserved.
   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
   Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
   Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the distribution.
   Neither the name of the NVIDIA CORPORATION nor the names of its contributors may be used to endorse
   or promote products derived from this software without specific prior written permission.
   
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
   OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/   
 
#ifndef _OPS_DEVICE_HXX_
#define _OPS_DEVICE_HXX_
#define DOUBLE double


/* Multiplication */


static __inline__ __device__ __host__  DOUBLE cuMul( DOUBLE x , DOUBLE y )
{
    return(x * y);
}



/* Negation */


static __inline__ __device__ __host__  DOUBLE cuNeg( DOUBLE x )
{
    return(-x);
}

/* Addition */


static __inline__ __device__ __host__  DOUBLE cuAdd( DOUBLE x , DOUBLE y )
{
    return(x + y);
}

 
/* Subtraction */


static __inline__ __device__ __host__  DOUBLE cuSub( DOUBLE x , DOUBLE y )
{
    return(x - y);
}

  
/* Division */


static __inline__ __device__ __host__  DOUBLE cuDiv( DOUBLE x , DOUBLE y )
{
    return (x / y);
}


/* Fma */


static __inline__ __device__ __host__  DOUBLE cuFma( DOUBLE x , DOUBLE y, DOUBLE d )
{
    return ((x * y) + d);
}

  

/* absolute value */


static __inline__ __device__ __host__  DOUBLE cuAbs( DOUBLE x )
{
    return (fabs(x));
}

 __inline__ __device__ __host__  DOUBLE cuGet(int x)
{
    return DOUBLE(x);

}

template <typename T_ELEM> struct __dynamic_shmem__{
    __device__ T_ELEM * getPtr() { 
        extern __device__ void error(void);
        error();
        return NULL;
    }
}; 
/* specialization of the above structure for the desired types */
template <> struct __dynamic_shmem__<float>{
    __device__ float * getPtr() { 
        extern __shared__ float Sptr[];
        return Sptr;
    }
};
template <> struct __dynamic_shmem__<DOUBLE>{
    __device__ DOUBLE * getPtr() { 
        extern __shared__ DOUBLE Dptr[];
        return Dptr;
    }
};




#endif
