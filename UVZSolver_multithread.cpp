/**
ULSAN NATIONAL INSTIUTE OF SCIENCE AND TECHNOLOGY
Copyright (c) 2019 HVCL lab
Created by Huong Nguyen
**/


#include "UVZSolver_multithread.h"
#include "support_funcs.h"
#include "gtsv.h"
// #include "spike_kernel.hxx"
// #define DOUBLE double
#define epsilon 1e-16

// m is size of the matrix, which is SN



__device__ void tridiag(int sn, DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD, DOUBLE *x, 
    DOUBLE *Ap, DOUBLE *Bp, DOUBLE *ep){
    if (sn == 0){
        x[0] = DD[0] / BB[0];
        // printf("sn is not set\n");
        return;
    }

    Ap[0] = - CC[0] / BB[0];
    Bp[0] = DD[0] / BB[0];
  
    for (int i = 1; i < sn; i++){
        
        ep[i] = AA[i] * Ap[i - 1] + BB[i]; 
        Ap[i] = -CC[i] / ep[i];
        Bp[i] = (DD[i] - (AA[i] * Bp[i - 1]) ) / ep[i];
    }
    
    x[sn] = (DD[sn] - (AA[sn] * Bp[sn - 1])) / (BB[sn] + (AA[sn] * Ap[sn - 1]));

    for (int i = sn - 1; i >= 0; i--){
        x[i] = Bp[i] + (Ap[i] * x[i + 1]);     

    } 
}



__global__ void  tridiagSolver_v2(bool print, bool isU, int startidx, int endidx, int jumpstep, int tridiag_coeff_width, Argument_Pointers* arg, Array_Pointers * arr){

    int i = threadIdx.x +  startidx;
    if (i > endidx) return;
    int number_of_segments;
    int *dau, *cuoi;
    if (isU){
        number_of_segments = arg->mocj[i];
        dau = arg->dauj;
    }
    else{ 
        number_of_segments = arg->moci[i];
        dau = arg->daui;

    }
    for (int j = 0; j < number_of_segments; j++){
        int first = dau[i * segment_limit + j];
        int pos = i * tridiag_coeff_width + first * jumpstep + jumpstep % 2; 

        DOUBLE* Dl = &(arr->AA[pos]);
        DOUBLE* D = &(arr->BB[pos]);
        DOUBLE* Du = &(arr->CC[pos]);
        DOUBLE* B = &(arr->DD[pos]);
        DOUBLE* x = &(arr->x[pos]);
        gtsv_spike_partial_diag_pivot_v1(Dl, D, Du, B, x, arr->SN[i * segment_limit + j]);
        
    }

}


__global__ void  tridiagSolver(bool print, bool isU, int startidx, int endidx, int jumpstep, int tridiag_coeff_width, Argument_Pointers* arg, Array_Pointers * arr){
    int i = blockIdx.y +  startidx;
    if (i > endidx) return;
    int number_of_segments;
    int *dau, *cuoi;
    if (isU){
        number_of_segments = arg->mocj[i];
        dau = arg->dauj;
    }
    else{ 
        number_of_segments = arg->moci[i];
        dau = arg->daui;

    }
    for (int j = 0; j < number_of_segments; j++){
        int first = dau[i * segment_limit + j];
        int pos = i * tridiag_coeff_width + first * jumpstep + jumpstep % 2; 

        DOUBLE* Dl = &(arr->AA[pos]);
        DOUBLE* D = &(arr->BB[pos]);
        DOUBLE* Du = &(arr->CC[pos]);
        DOUBLE* B = &(arr->DD[pos]);
        DOUBLE* x = &(arr->x[pos]);
        DOUBLE* Ap = &(arr->Ap[pos]);
        DOUBLE* Bp = &(arr->Bp[pos]);
        DOUBLE* ep = &(arr->ep[pos]);
        // if ((threadIdx.x ==0))
        //     printf("i= %d, j = %d, SN_offset = %d, sn = %d\n", i, j, i*segment_limit + j, arr->SN[i * segment_limit + j] );
        tridiag(arr->SN[i * segment_limit + j], Dl, D, Du, B, x, Ap, Bp, ep);

    }
}
__device__ void bienrandau(int i, int first, int last,  DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD,
    DOUBLE *a1, DOUBLE *b1, DOUBLE *c1, DOUBLE *d1, DOUBLE *a2, DOUBLE *c2, DOUBLE *d2){
    // printf("bienrandau is calledsn");
    if (first > last) return;
    AA[i * 2] = a2[i];
    AA[i * 2 + 1] = a1[i];
    BB[i * 2] = 1;
    BB[i * 2 + 1] = b1[i];
    CC[i * 2] = c2[i];
    CC[i * 2 + 1] = c1[i];
    DD[i * 2] = d2[i];
    DD[i * 2 + 1] = d1[i];
}


__device__ void bienlongdau(int i, int first, int last,  DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE*DD,
    DOUBLE *a1, DOUBLE *b1, DOUBLE *c1, DOUBLE *d1, DOUBLE *a2, DOUBLE *c2, DOUBLE *d2){

    if (first > last) return;
    AA[i * 2] = a1[i];
    AA[i * 2 + 1] = a2[i + 1];
    BB[i * 2] = b1[i]; 
    BB[i * 2 + 1] = 1;
    CC[i * 2] = c1[i];
    CC[i * 2 + 1] = c2[i + 1];
    DD[i * 2] = d1[i];
    DD[i * 2 + 1] = d2[i + 1];
    printf("AA %llx, BB %llx, CC %llx, DD %llx\n",AA, BB, CC, DD);

}


__device__ void  _vzSolver_calculate_preindex(int i, int j, int width, int first, int last,  Argument_Pointers* arg, Array_Pointers *arr, Constant_Coeffs* coeffs){
    if ((first > last) || (j < first) || ( j > last)) return;
    DOUBLE p = 0.0;
    DOUBLE q = 0.0;    
    __shared__ DOUBLE *u, *v, *z, *Htdv, *Htdu, *H_moi, *VISCOIDX, *Tsyw, *Ky1, *a2, *c2, *d2, *f1, *f2, *f3, *f5;
    __shared__ int support_array_width;
    __shared__ DOUBLE H_TINH, dYbp, dT, dTchia2dY,dTchia2dX, dXbp, dX2, CORIOLIS_FORCE, Windy ;
    support_array_width = arg->M + 2;
    u = arg->u;
    v = arg->v;
    z = arg->z;
    Htdv = arg->Htdv;
    Htdu = arg->Htdu;
    VISCOIDX = arg->VISCOIDX;
    Tsyw = arg->Tsyw;
    H_moi = arg->H_moi;
    Ky1 = arg-> Ky1;
    a2 = &(arr->a2[i * support_array_width]);
    c2 = &(arr->c2[i * support_array_width]);
    d2 = &(arr->d2[i * support_array_width]);
    f1 = &(arr->f1[i * support_array_width]);
    f2 = &(arr->f2[i * support_array_width]);
    f3 = &(arr->f3[i * support_array_width]);
    f5 = &(arr->f5[i * support_array_width]);

    H_TINH = coeffs->H_TINH;
    dT = coeffs->dT;
    dXbp = coeffs->dXbp;
    dYbp = coeffs->dYbp;
    dX2 = coeffs->dX2;
    dTchia2dY = coeffs->dTchia2dY;
    dTchia2dX = coeffs->dTchia2dX;
    Windy = coeffs->Windy;
    CORIOLIS_FORCE = coeffs->CORIOLIS_FORCE;
    // can phai control dk bien o ngoai function
    // ham nay ko coarseless
    // can actually get rid of fs, and only need to use 3 variables instead of 3 m * n arrays
    // reduce a significant number of read and write operation to global memory
    // make more spaces to accomodate larger data set
    DOUBLE utb = (u[(i - 1) * width + j] + u[i * width + j] + u[(i - 1) * width + j + 1] + u[i * width + j + 1]) * 0.25;
    f1[j] = dTchia2dY * v[i * width + j] + VISCOIDX[i * width + j] * dT / dYbp;
    f2[j] = -(2 + Ky1[i * width + j] * dT * sqrt(v[i * width + j] * v[i * width + j] + utb * utb) / Htdv[i * width + j] + (2 * dT * VISCOIDX[i * width + j]) / dYbp);
    f3[j] = dT * VISCOIDX[i * width + j] / dYbp - dTchia2dY * v[i * width + j];
    // printf("%d f1[%d] = %.10f\n", i, j, f1[j]);
    if (H_moi[(i - 1) * width + j] <= H_TINH){
        if (utb < 0){
            q = utb * (-3 * v[i * width + j] + 4 * v[(i + 1) * width + j] - v[(i + 2) * width + j]) / dX2;
            p = (v[i * width + j] - 2 * v[(i + 1) * width + j] + v[(i + 2) * width + j] ) / dXbp;
        }
    }else{
        if (H_moi[(i + 1) * width + j] <= H_TINH){
            if ((H_moi[(i - 2) * width + j] > H_TINH) && (utb > 0)){
                q = utb * (3 * v[i * width + j] - 4 * v[(i - 1) * width + j] + v[(i - 2) * width + j]) / dX2;
                p = (v[i * width + j] - 2 * v[(i - 1) * width + j] + v[(i - 2) * width + j] ) / dXbp;
            }
        }else{
            q = utb * (v[(i + 1) * width + j] - v[(i - 1) * width + j]) / dX2;
            p = (v[(i + 1) * width + j] - 2 * v[i * width + j] + v[(i - 1) * width + j]) / dXbp;
        }
    }
    f5[j] = -2 * v[i * width + j] + dT * q + dT * CORIOLIS_FORCE * utb - dT * VISCOIDX[i * width + j] * p - dT * (Windy - Tsyw[i * width + j]) / Htdv[i * width + j];
    c2[j] = dTchia2dY * Htdv[i * width + j];             
    a2[j] = - dTchia2dY * Htdv[i * width + j - 1];
    d2[j] = z[i * width + j] - dTchia2dX * (Htdu[i * width + j] * u[i * width + j] - Htdu[(i - 1) * width + j] * u[(i - 1) * width + j]);


}

__device__ void  _uzSolver_calculate_preindex(int i, int j, int width, int first, int last, Argument_Pointers* arg, Array_Pointers *arr, Constant_Coeffs* coeffs){
    if ((first > last) || (i < first) || ( i > last)) return;
    DOUBLE p = 0.0;
    DOUBLE q = 0.0;
    __shared__ DOUBLE *u, *v, *z, *Htdv, *Htdu, *H_moi, *VISCOIDX, *Tsxw, *Kx1, *a2, *c2, *d2, *f1, *f2, *f3, *f5;
    __shared__ int support_array_width;
    __shared__ DOUBLE H_TINH, dYbp, dXbp, dT, dTchia2dY,dTchia2dX, dY2, CORIOLIS_FORCE, Windx ;

    support_array_width = arg->N + 2;
    u = arg->u;
    v = arg->v;
    z = arg->z;
    Htdv = arg->Htdv;
    Htdu = arg->Htdu;
    VISCOIDX = arg->VISCOIDX;
    Tsxw = arg->Tsxw;
    H_moi = arg->H_moi;
    Kx1 = arg-> Kx1;
    a2 = &(arr->a2[j * support_array_width]);
    c2 = &(arr->c2[j * support_array_width]);
    d2 = &(arr->d2[j * support_array_width]);
    f1 = &(arr->f1[j * support_array_width]);
    f2 = &(arr->f2[j * support_array_width]);
    f3 = &(arr->f3[j * support_array_width]);
    f5 = &(arr->f5[j * support_array_width]);

    H_TINH = coeffs->H_TINH;
    dT = coeffs->dT;
    dXbp = coeffs->dXbp;
    dYbp = coeffs->dYbp;
    dY2 = coeffs->dY2;
    dTchia2dY = coeffs->dTchia2dY;
    dTchia2dX = coeffs->dTchia2dX;
    Windx = coeffs->Windx;
    CORIOLIS_FORCE = coeffs->CORIOLIS_FORCE;

    DOUBLE vtb = (v[i * width +  j - 1] + v[i * width + j] + v[(i + 1) * width + j - 1] + v[(i + 1) * width + j]) * 0.25;
    f1[i] = dTchia2dX * u[i * width + j] + VISCOIDX[i * width + j] * dT / dXbp;
    f2[i] = -(2 + Kx1[i * width + j] * dT * sqrt(u[i * width + j] * u[i * width + j] + vtb * vtb) / Htdu[i * width + j] + (2 * dT * VISCOIDX[i * width + j]) / dXbp); // chua tinh muc nuoc trung binh

    f3[i] = dT * VISCOIDX[i * width + j] / dXbp - dTchia2dX * u[i * width + j];
    
    if (H_moi[i * width + j - 1] <= H_TINH){
        if (vtb < 0){
            p = vtb * (-3 * u[i * width + j] + 4 * u[i * width + j + 1] - u[i * width + j + 2]) / dY2;
            q = (u[i * width + j] - 2 * u[i * width + j + 1] + u[i * width + j + 2]) / dYbp;
        }
    }
    else{
        if (H_moi[i * width + j + 1] <= H_TINH){
            if ((H_moi[i * width + j - 2] > H_TINH) && (vtb > 0)){
                p = vtb * (3 * u[i * width + j] - 4 * u[i * width + j - 1] + u[i * width + j - 2]) / dY2;
                q = (u[i * width + j] - 2 * u[i * width + j - 1] + u[i * width + j - 2] ) / dYbp;
            }
        }else{
            p = vtb * (u[i * width + j + 1] - u[i * width + j - 1]) / dY2;
            q = (u[i * width + j + 1] - 2 * u[i * width + j] + u[i * width + j - 1]) / dYbp;
        }
    }

    f5[i] = -2 * u[i * width + j] + dT * p  - dT * CORIOLIS_FORCE * vtb - dT * VISCOIDX[i * width + j] * q - dT * (Windx - Tsxw[i * width + j]) / Htdu[i * width + j];

    
    c2[i] = dTchia2dX * Htdu[i * width + j];
    a2[i] = - dTchia2dX * Htdu[(i - 1) * width + j];

    d2[i] = z[i * width + j] - dTchia2dY * (Htdv[i * width + j] * v[i * width + j] - Htdv[i * width + j - 1] * v[i * width + j - 1]);
    
}

__device__ void _calculate_abcd(int i, int j, int first, int last, DOUBLE f4, int support_array_width,  bool bienran1, bool bienran2, Array_Pointers* arr){
    

    if ((first > last) || (j < first) || ( j > last)) return;

    __shared__ DOUBLE *a1, *b1, *c1, *d1, *a2, *c2, *d2, *f1, *f2, *f3, *f5;
    a1 = &(arr->a1[i * support_array_width]);
    b1 = &(arr->b1[i * support_array_width]);
    c1 = &(arr->c1[i * support_array_width]);
    d1 = &(arr->d1[i * support_array_width]);
    a2 = &(arr->a2[i * support_array_width]);
    c2 = &(arr->c2[i * support_array_width]);
    d2 = &(arr->d2[i * support_array_width]);
    f1 = &(arr->f1[i * support_array_width]);
    f2 = &(arr->f2[i * support_array_width]);
    // f4 = 2 * g * dTchia2dY;
    f3 = &(arr->f3[i * support_array_width]);
    f5 = &(arr->f5[i * support_array_width]);



    if (last - 1 > first){
            {
                a1[j] = f4 - f1[j] / a2[j];
                c1[j] = -f4 - (f3[j] / c2[j + 1]);

                b1[j] = f2[j] - (f1[j] * c2[j] / a2[j]) - (f3[j] * a2[j + 1] / c2[j + 1]);
                
                d1[j] = f5[j] - (f1[j] * d2[j] / a2[j]) - (f3[j] * d2[j + 1] / c2[j + 1]);

            }
        //limit this for the warp with tid range contains first only
        __syncthreads();
        if (bienran1){
            a1[first] = f4;
            b1[first] = f2[first] - (f3[first] * a2[first+ 1] / c2[first+ 1]);
            c1[first] = - f4 - (f3[first] / c2[first+ 1]);
            d1[first] = f5[first] - (f3[first] * d2[first+ 1] / c2[first+ 1]);
            // if (i == 245){
            //     printf("d1[first] %d %d %d %llx %llx %llx\n",j, first, last, f5[first], f3[first], d2[first+ 1]);

            // }
        }
        // limit this for the warp with tid range contains last only
        if (bienran2){
            a1[last - 1] = f4 - f1[last -1] / a2[last -1];
            b1[last - 1] = f2[last - 1] - (f1[last - 1] * c2[last -1] / a2[last - 1]);
            c1[last - 1] = - f4;
            d1[last - 1] = f5[last - 1] - (f1[last - 1] * d2[last - 1] / a2[last - 1]);

        }
    }
    else{
        //printf("Overflow, comparision btw DOUBLE and int\n");
        a1[first] = f4;
        b1[first] = f2[first];
        d1[first] = f5[first];
        c1[first] = -f4;
    }

    // if (i == 245) printf("_calculate_abcd %d %llx %llx %llx %llx %llx \n",j, d1[j], f1[j], d2[j], f3[j], a2[j]);

}


__device__ void _calculate_matrix_coeff(bool isU, int i, int j, int support_array_width, int tridiag_coeff_width, int first, int last, int seg_no, bool bienran1, bool bienran2, 
    DOUBLE ubt_or_vbd, DOUBLE ubp_or_vbt,  DOUBLE TZ_r, DOUBLE TZ_l, bool dkBienQ_1, bool dkBienQ_2, int dkfr, Argument_Pointers* arg, Array_Pointers *arr){
    //DOUBLE *a1, DOUBLE *b1, DOUBLE *c1, DOUBLE *d1, DOUBLE *a2, DOUBLE *b2, DOUBLE *c2, DOUBLE *d2, DOUBLE* AA, DOUBLE* BB, DOUBLE* CC, DOUBLE* DD){
    


    if ((first > last) || (j < first) || ( j > last)) {
        return;
    }

    // int array_width = arg->M + 2;
    // int tridiag_coeff_width = 2 * arg->M  + 1;
    __shared__ DOUBLE *a1, *b1, *c1, *d1, *a2, *c2, *d2, *AA, *BB, *CC, *DD;
    int offset = first * 2;
    a1 = &(arr->a1[i * support_array_width]);
    b1 = &(arr->b1[i * support_array_width]);
    c1 = &(arr->c1[i * support_array_width]);
    d1 = &(arr->d1[i * support_array_width]);
    a2 = &(arr->a2[i * support_array_width]);
    c2 = &(arr->c2[i * support_array_width]);
    d2 = &(arr->d2[i * support_array_width]);
    AA = &(arr->AA[i * tridiag_coeff_width]);
    BB = &(arr->BB[i * tridiag_coeff_width]);
    CC = &(arr->CC[i * tridiag_coeff_width]);
    DD = &(arr->DD[i * tridiag_coeff_width]);
    int sn = 2 * (last - first);
    bool isBienran;

    // if (j == first)
    //     printf("debug %d %x\n", i, BB);
    // TODO: can re-organize this function so that only warp with right index range execute right snippet, avoid unnescessary memory access and execution cost.

    if (bienran1){
        
        bienrandau(j, first, last, AA, BB, CC, DD, a1, b1, c1, d1, a2, c2, d2);

        DD[offset] = d2[first];
        // ran - long
        if (bienran2 == false){
            if ((dkBienQ_2) && (last == dkfr)){         //r == dkfr:   // Kiem tra lai phan nay
                // if (j == last)
                // printf("herrrrreeeeee giai vz ran long q %d \n", last);
                // attention 
                if (j == last){
                    AA[sn + offset] = a2[last];
                    BB[sn + offset] = 1;
                    DD[sn + offset] = d2[last] - c2[last] * ubp_or_vbt;
                }

            }
            else{
                //printf("ran - long \n");
                sn --;
                if (j == last){
                    AA[sn + offset] = a1[last - 1];
                    BB[sn + offset] = b1[last - 1];
                    DD[sn + offset] = d1[last - 1] - c1[last - 1] * TZ_r;
                }

            }
        }
        // ran - ran
        else{
            if (j == last){
                AA[sn + offset] = a2[last];
                BB[sn + offset] = 1;
                DD[sn + offset] = d2[last];
            }
        }
    }
    // long
    else{
        if ((dkBienQ_1) && (first == 2)){
            // printf("debug\n");
            bienrandau(j,first, last,  AA, BB, CC, DD, a1, b1, c1, d1, a2, c2, d2);

            DD[offset] = d2[first] - a2[first] * ubt_or_vbd;
            // thieu bb[1] va cc[1] cho truong hop vz, hoi lai co
            isBienran = true;
        }
        else{
            bienlongdau(j,first, last, AA, BB, CC, DD, a1, b1, c1, d1, a2, c2, d2);
            if (j == first){
                BB[offset] = b1[first];
                CC[offset] = c1[first];
                DD[offset] = d1[first] - a1[first] * TZ_l;
            }
            isBienran = false;
        }
        // long - long
        // TODO: this part should be done by warp with index range contain sn.

        if (bienran2 == false){ // variable isbienran is equivalent with variable text in original code
            if ((dkBienQ_1) && (last == dkfr)){     //r == dkfr: // BienQ[0] && r == M trong truong hop giaianv
                if (!isBienran) sn --;
                if(j == last){
                    AA[sn + offset] = a2[last];
                    BB[sn + offset] = 1;
                    DD[sn + offset] = d2[last] - c2[last] * ubp_or_vbt;
                }
            }
            else{
                sn --;
                if (!isBienran)
                    sn -= 1;
                if (j == last){
                    AA[sn + offset] = a1[last - 1];
                    BB[sn + offset] = b1[last - 1];
                    DD[sn + offset] = d1[last - 1] - c1[last - 1] * TZ_r;
                }
            }
        }
        else{
            if (!isBienran) sn --;
            if (j == last){
                AA[sn + offset] = a2[last];
                BB[sn + offset] = 1;
                DD[sn + offset] = d2[last];
            }
        }
    }
    if (j == last){
        // printf("first = %d, last = %d, sn = %d, i  = %d, segment_limit = %d seg_no = %d\n",first, last, sn, i, segment_limit, seg_no );
        arr->SN[i * segment_limit + seg_no] = sn;
    }
}

__device__ void _vzSolver_extract_solution( int i,int j, int sn, int width, int first, int last, bool bienran1, bool bienran2, Argument_Pointers *arg, Array_Pointers * arr)
{   
    int M = arg->M;
    DOUBLE* t_z = arg->t_z;
    DOUBLE* t_v = arg->t_v;
    int* bienQ = arg->bienQ;
    DOUBLE* x = &(arr->x[i * (2 * M + 1)]);
    DOUBLE* a2 = &(arr->a2[i * (M + 2)]);
    DOUBLE* c2 = &(arr->c2[i * (M + 2)]);
    DOUBLE* d2 = &(arr->d2[i * (M + 2)]);
    // if ( j == 10){
    //     printf("%d %f \n", i, arg->t_u[i * width + j] );
    // }
     
    if (j > last || j < first) return;

    if (j < last){
        if (bienran1){
            t_z[i * width + j] = x[2 * j];
            t_v[i * width + j] = x[2 * j + 1];
            t_v[i * width + first - 1] = 0;

        }
        else{
            if( (bienQ[1]) && (first == 2)){
                t_z[i * width + j] = x[2 * j];
                t_v[i * width + j] = x[2 * j + 1];
                t_v[i * width + first - 1] = arg->vbd[i];
            }
            else{
                // warp divergent here. TODO
                // TODO: first warp can use one instruction other than other warps
                if (j == first)
                {       
                    t_v[i * width + first] = x[2 * first];
                    t_v[i * width + first - 1] = (d2[first] - t_z[i * width + first] - c2[first] * t_v[i * width + first]) / a2[first];
                }
                else{
                    t_z[i * width + j] = x[2 * j - 1];
                    t_v[i * width + j] = x[2 * j];
                }
                
            }
        }
    }  else {

        if (bienran2){
            
            t_v[i * width + last] = 0;
            t_z[i * width + last] = x[first * 2 + sn];
            // if (i == arg->N)
            //     printf("t_v[%d, %d] initially is %.15lf\n", i, j, t_v[i * width + last]);
        }
        else{

            if ((bienQ[0]) && (last == M)){
                t_v[i * width + last] = arg->vbt[i];
                // printf("t_v[%d, %d] = %.15lf\n",i, last, t_v[i * width  + last] );
                t_z[i * width + last] = x[2 * first + sn];

            }
            else{
                t_v[i * width + last] = (d2[last] - a2[last] * t_v[i * width + last - 1] - t_z[i * width + last]) / c2[last];
            }
        }
    }
    // if ((t_v[i * width + j] > 0) || (t_z[i * width + j] > 0)){
    //     printf("t_v[%d, %d] = %.9f\n", i, j, t_v[i * width + j]);
    //     printf("t_z[%d, %d] = %.9f\n", i, j, t_z[i * width + j]);
    // }


}


__device__ void _uzSolver_extract_solution( int i, int j, int sn, int width, int first, int last, bool bienran1, bool bienran2, Argument_Pointers *arg, Array_Pointers * arr)
{   
    int N = arg->N;
    DOUBLE* t_z = arg->t_z;
    DOUBLE* t_u = arg->t_u;
    int* bienQ = arg->bienQ;
    DOUBLE* x = &(arr->x[j * (2 * N + 1)]);
    DOUBLE* a2 = &(arr->a2[j * (N + 2)]);
    DOUBLE* c2 = &(arr->c2[j * (N + 2)]);
    DOUBLE* d2 = &(arr->d2[j * (N + 2)]);
    if (i > last || i < first) return;
    if (i < last){
        if (bienran1){

            // t_z[i * width + j] = x[2 * (i - first)];
            // t_u[i * width + j] = x[2 * (i - first) + 1];
            t_z[i * width + j] = x[2 * i];
            t_u[i * width + j] = x[2 * i + 1];
            t_u[(first - 1) * width + j] = 0;
        }else{
            if ((bienQ[2]) && (first == 2)){
                t_z[i * width + j] = x[2 * i];
                t_u[i * width + j] = x[2 * i + 1];
                t_u[(first - 1) * width + j] = arg->ubt[j];
            }
            else{
                if (i == first){
                    t_u[first * width + j] = x[first * 2];
                    t_u[(first - 1) * width + j] = (d2[first] - t_z[first * width + j] - c2[first] * t_u[first * width + j]) / a2[first];

                } else if (i < last){
                    t_z[i * width + j] = x[2 * i - 1];
                    t_u[i * width + j] = x[2 * i];    
                }
                
            }
        }       
    }
    else { //if (i == last)
        if (bienran2 ){
            t_u[last * width + j] = 0;
            t_z[last * width + j] = x[first * 2 + sn];

        }
        else{
            if ((bienQ[3]) && (last == N)){
                t_u[last * width + j] = arg->ubp[j];
                t_z[last * width + j] = x[first * 2 + sn];
            }

            else{
                t_u[last * width + j] = (d2[last] - a2[last] * t_u[(last - 1) * width + j] - t_z[last * width + j]) / c2[last];
            }
        }
    }

}



__device__ void vSolver(DOUBLE t, int offset, int first, int last, int row, int col, bool bienran1, bool bienran2, DOUBLE* VISCOIDX, DOUBLE* Tsyw, 
    DOUBLE *v, DOUBLE *t_v, DOUBLE *u, DOUBLE *t_u, DOUBLE *z, DOUBLE *t_z, DOUBLE *Ky1, DOUBLE *Htdv, DOUBLE *H_moi, Constant_Coeffs* coeffs)
{

    DOUBLE p, q, tmp;
    __shared__ DOUBLE H_TINH, dY2, dX2, dY, CORIOLIS_FORCE, Windy, dXbp, dYbp, g, HaiChiadT;

    H_TINH = coeffs->H_TINH; 
    dY = coeffs->dY;
    dY2 = coeffs->dY2;
    dX2 = coeffs->dX2;
    Windy = coeffs->Windy;
    dXbp = coeffs->dXbp;
    dYbp = coeffs->dYbp;
    CORIOLIS_FORCE = coeffs->CORIOLIS_FORCE;
    g = coeffs->g;
    HaiChiadT = coeffs->HaiChiadT;


    q = 0;
    p = 0;
    tmp = 0;
    DOUBLE utb = (u[(row - 1) * offset +  col] + u[row * offset +  col] + u[(row - 1) * offset +  col + 1] + u[row * offset +  col + 1]) * 0.25;
    DOUBLE t_utb = (t_u[(row - 1) * offset +  col] + t_u[row * offset +  col] + t_u[(row - 1) * offset +  col + 1] + t_u[row * offset +  col + 1]) * 0.25;

    p = (v[row * offset +  col + 1] - v[row * offset +  col - 1]) / dY2;
    p = HaiChiadT + p + Ky1[row * offset +  col] * sqrt(utb * utb + v[row * offset +  col] * v[row * offset +  col]) / Htdv[row * offset +  col];
    
    
    if (H_moi[(row - 1) * offset +  col] <= H_TINH){
        if (utb < 0){
            q = t_utb * (-3 * v[row * offset +  col] + 4 * v[(row + 1) * offset +  col] - v[(row + 2) * offset +  col]) / dX2;
            // q = tutb * (-3 * v(i, j) + 4 * v(i + 1, j) - v(i + 2, j)) / dX2
            tmp = (v[row * offset +  col] - 2 * v[(row + 1) * offset +  col] + v[(row + 2) * offset +  col] ) / dXbp;
        }
    }
    else{
        if (H_moi[(row + 1) * offset +  col] <= H_TINH){
            if ((H_moi[(row - 2) * offset +  col] > H_TINH) && (utb > 0)){
                q = t_utb * (3 * v[row * offset +  col] - 4 * v[(row - 1) * offset +  col] + v[(row - 2) * offset +  col]) /dX2;
                tmp = (v[row * offset +  col] - 2 * v[(row - 1) * offset +  col] + v[(row - 2) * offset +  col] ) / dXbp;
            }
        }
        else{
            q = t_utb * (v[(row + 1) * offset +  col] - v[(row - 1) * offset +  col]) / dX2;
            tmp = (v[(row + 1) * offset +  col] - 2 * v[row * offset +  col] + v[(row - 1) * offset +  col]) / dXbp;
        }
    }
    q = HaiChiadT * v[row * offset +  col] - q - CORIOLIS_FORCE * t_utb;

    q = (q - g * (z[row * offset +  col + 1] - z[row * offset +  col]) / dY + 
        VISCOIDX[row * offset +  col] * (tmp + (v[row * offset +  col + 1] - 2 * v[row * offset +  col] + v[row * offset +  col - 1]) / dYbp)) + 
        (Windy - Tsyw[row * offset +  col]) / Htdv[row * offset +  col];
    

    t_v[row * offset +  col] = q / p ;
}


// __device__ void uSolver(DOUBLE t, int offset, int N, int first, int last, int row, int col, bool bienran1, bool bienran2, DOUBLE* VISCOIDX, DOUBLE* Tsxw,
//     DOUBLE *v, DOUBLE *t_v, DOUBLE *u, DOUBLE *t_u, DOUBLE *z, DOUBLE *t_z, DOUBLE *Kx1, DOUBLE *Htdu, DOUBLE *H_moi, Constant_Coeffs* coeffs)
// {

//     __shared__ DOUBLE H_TINH, dY2, dX2, dX, CORIOLIS_FORCE, Windx, dXbp, dYbp, g, HaiChiadT;

//     H_TINH = coeffs->H_TINH; 
//     dX = coeffs->dY;
//     dY2 = coeffs->dY2;
//     dX2 = coeffs->dX2;
//     Windx = coeffs->Windy;
//     dXbp = coeffs->dXbp;
//     dYbp = coeffs->dYbp;
//     CORIOLIS_FORCE = coeffs->CORIOLIS_FORCE;
//     g = coeffs->g;
//     HaiChiadT = coeffs->HaiChiadT;



//     DOUBLE p, q, tmp;
//     p = 0; q = 0; tmp = 0;
//     // this can be optimised

//     DOUBLE vtb = (v[row * offset + col - 1] + v[row * offset + col] + v[(row + 1) * offset +  col - 1] + v[(row + 1) * offset +  col]) * 0.25;
//     DOUBLE t_vtb = (t_v[row * offset + col - 1] + t_v[row * offset + col] + t_v[(row + 1) * offset +  col - 1] + t_v[(row + 1) * offset +  col]) * 0.25;

//     p = (u[(row + 1) * offset +  col] - u[(row - 1) * offset +  col]) / dX2;

//     p = (HaiChiadT + p + Kx1[row * offset + col] * sqrt(vtb * vtb + u[row * offset + col] * u[row * offset + col]) / Htdu[row * offset + col]);
    
//     if (H_moi[row * offset + col - 1] <= H_TINH){
//         if (vtb < 0){
//             q = t_vtb * (-3 * u[row * offset + col] + 4 * u[row * offset + col + 1] - u[row * offset +  col + 2]) / dY2;
//             tmp = (u[row * offset +  col] - 2 * u[row * offset +  col + 1] + u[row * offset +  col + 2] ) / dYbp;
//         }
//     }
//     else {
//         if (H_moi[row * offset +  col + 1] <= H_TINH) {
//             if ((H_moi[row * offset +  col - 2] > H_TINH) && (vtb > 0)){ 
//                 q = t_vtb * (3 * u[row * offset +  col] - 4 * u[row * offset +  col - 1] + u[row * offset +  col - 2]) /dY2;
//                 tmp = (u[row * offset +  col] - 2 * u[row * offset +  col - 1] + u[row * offset +  col - 2] ) / dYbp;
//             }
//         }
//         else{
            
//             q = t_vtb * (u[row * offset +  col + 1] - u[row * offset +  col - 1]) / dY2;
//             tmp = (u[row * offset +  col + 1] - 2 * u[row * offset +  col] + u[row * offset +  col - 1]) / dYbp;
//         }
//     }


//     q = HaiChiadT * u[row * offset +  col] - q + CORIOLIS_FORCE * t_vtb;
//     q = (q - g * (z[(row + 1) * offset +  col] - z[row * offset +  col]) / dX + VISCOIDX[row * offset +  col] * (tmp + (u[(row + 1) * offset +  col] - 
//                         2 * u[row * offset +  col] + u[(row - 1) * offset +  col]) / dXbp )) + (Windx - Tsxw[row * offset +  col]) / Htdu[row * offset +  col];

//     t_u[row * offset +  col] = q / p;
// }

__device__ void uSolver(DOUBLE t, int offset, int N, int first, int last, int row, int col, bool bienran1, bool bienran2, Argument_Pointers* arg, Constant_Coeffs* coeffs)
{

    __shared__ DOUBLE H_TINH, dY, dX, CORIOLIS_FORCE, Windx, g, HaiChiadT;
    __shared__ DOUBLE *VISCOIDX, *Tsxw, *v, *t_v, *u, *t_u, *z, *Kx1, *Htdu, *H_moi;

    H_TINH = coeffs->H_TINH; 
    dX = coeffs->dY;
    dY = coeffs->dY;
    //dX2 = coeffs->dX2;
    Windx = coeffs->Windy;
    //dXbp = coeffs->dXbp;
    //dYbp = coeffs->dYbp;
    CORIOLIS_FORCE = coeffs->CORIOLIS_FORCE;
    g = coeffs->g;
    HaiChiadT = coeffs->HaiChiadT;

    VISCOIDX = arg->VISCOIDX;
    Tsxw = arg->Tsxw;
    v = arg->v;
    u = arg->u;
    z = arg->z;
    t_u = arg->t_u;
    t_v = arg->t_v;
    //t_z = arg->t_z;
    Kx1 = arg->Kx1;
    Htdu = arg->Htdu;
    H_moi = arg->H_moi;



    DOUBLE p, q, tmp;
    p = 0; q = 0; tmp = 0;
    // this can be optimised
    int i_j = row * offset + col;
    // [i+1, j]
    int i_1_j = i_j + offset;
    DOUBLE vtb = (v[i_j - 1] + v[i_j] + v[i_1_j - 1] + v[i_1_j]) * 0.25;

    DOUBLE t_vtb = (t_v[i_j - 1] + t_v[i_j] + t_v[i_1_j - 1] + t_v[i_1_j]) * 0.25;

    p = (u[i_1_j] - u[(row - 1) * offset +  col]) / (dX * 2);

    p = (HaiChiadT + p + Kx1[i_j] * sqrt(vtb * vtb + u[i_j] * u[i_j]) / Htdu[i_j]);
    
    if (H_moi[i_j - 1] <= H_TINH){
        if (vtb < 0){
            q = t_vtb * (-3 * u[i_j] + 4 * u[i_j + 1] - u[i_j + 2]) / (dY * 2);
            tmp = (u[i_j] - 2 * u[i_j + 1] + u[i_j + 2] ) / (dY * dY);
        }
    }
    else {
        if (H_moi[i_j + 1] <= H_TINH) {
            if ((H_moi[i_j - 2] > H_TINH) && (vtb > 0)){ 
                q = t_vtb * (3 * u[i_j] - 4 * u[i_j - 1] + u[i_j - 2]) / (dY * 2);
                tmp = (u[i_j] - 2 * u[i_j - 1] + u[i_j - 2] ) / (dY * dY);
            }
        }
        else{
            
            q = t_vtb * (u[i_j + 1] - u[i_j - 1]) / (dY * 2);
            tmp = (u[i_j + 1] - 2 * u[i_j] + u[i_j - 1]) / (dY * dY);
        }
    }


    q = HaiChiadT * u[i_j] - q + CORIOLIS_FORCE * t_vtb;
    q = (q - g * (z[i_1_j] - z[i_j]) / dX + VISCOIDX[i_j] * (tmp + (u[i_1_j] - 
                        2 * u[i_j] + u[(row - 1) * offset +  col]) / (dX * dX) )) + (Windx - Tsxw[i_j]) / Htdu[i_j];

    t_u[i_j] = q / p;
}


__global__ void 
VZSolver_calculate_preindex( int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
    // locate segment 
    int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if ((i > endidx) || (j > arg->M)) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    locate_segment_v(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->daui, arg->cuoii, arg->moci, arg->h, coeffs->NANGDAY);
    _vzSolver_calculate_preindex(i, j, arg->M + 3, first, last, arg, arr, coeffs);
}



__global__ void 
VZSolver_calculate_abcd(int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
    int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if ((i > endidx) || (j > arg->M)) return;

    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    DOUBLE dTchia2dY, g, NANGDAY;
    dTchia2dY = coeffs->dTchia2dY;
    g = coeffs->g;
    NANGDAY = coeffs->NANGDAY;
    locate_segment_v(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->daui, arg->cuoii, arg->moci, arg->h, NANGDAY);
    _calculate_abcd(i, j, first, last, 2 * g * dTchia2dY, arg->M + 2, bienran1 , bienran2, arr);

}

__global__ void 
VZSolver_calculate_matrix_coeff(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg, Array_Pointers* arr){
    int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if ((i > endidx) || (j > arg->M)) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    int seg_no = locate_segment_v(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->daui, arg->cuoii, arg->moci, arg->h, NANGDAY);

    _calculate_matrix_coeff(false,i, j, arg->M + 2, 2 * (arg->M) + 1, first, last, seg_no, bienran1, bienran2,
                          arg->vbd[i], arg->vbt[i], arg->t_z[i * (arg->M + 3) + last], arg->t_z[i * (arg->M + 3) + first], arg->bienQ[1], arg->bienQ[0], arg->M,arg,arr);
}

__global__ void 
VZSolver_extract_solution(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg, Array_Pointers* arr){

    int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;

    if ((j > endidx) || (i > arg->N)) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    int seg_no = locate_segment_v(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->daui, arg->cuoii, arg->moci, arg->h, NANGDAY);
    _vzSolver_extract_solution(i, j, arr->SN[i * segment_limit + seg_no], arg->M + 3, first, last, bienran1, bienran2, arg, arr);


}

__global__ void 
UZSolver_extract_solution(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg, Array_Pointers* arr){
    int i = blockIdx.y * blockDim.y + threadIdx.y + 2;
    int j = blockIdx.x * blockDim.x + threadIdx.x + startidx;
    if ((j > endidx) || (i > arg->N)) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    int seg_no = locate_segment_u(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->dauj, arg->cuoij, arg->mocj, arg->h, NANGDAY);
    
    _uzSolver_extract_solution(i, j, arr->SN[j * segment_limit + seg_no], arg->M + 3, first, last, bienran1, bienran2, arg, arr);
   

}

__global__ void 
UZSolver_calculate_preindex(int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
    // locate segment 

    int i = blockIdx.y * blockDim.y + threadIdx.y + 2;
    int j = blockIdx.x * blockDim.x + threadIdx.x + startidx;

    if ((j > endidx) || (i > arg->N)) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    locate_segment_u(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->dauj, arg->cuoij, arg->mocj, arg->h, coeffs->NANGDAY);

    _uzSolver_calculate_preindex( i, j, arg->M + 3, first, last, arg, arr, coeffs);
}

// same i are in same block

__global__ void 
UZSolver_calculate_abcd(int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
    // i runs from start index to M 
    int i = blockIdx.y * blockDim.y + threadIdx.y + 2;
    int j = blockIdx.x * blockDim.x + threadIdx.x + startidx;
    if ((j > endidx) || (i > arg->N)) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    DOUBLE dTchia2dX, g, NANGDAY;
    dTchia2dX = coeffs->dTchia2dX;
    g = coeffs->g;
    NANGDAY= coeffs->NANGDAY;

    locate_segment_u(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->dauj, arg->cuoij, arg->mocj, arg->h, NANGDAY);
    // if ((first) && (last))
    //     printf("i = %d, j = %d, first = %d, last = %d\n",i, j, first, last );
    _calculate_abcd(j, i, first, last, 2 * g * dTchia2dX, arg->N + 2, bienran1 , bienran2, arr);

}

// same i are in same block
__global__ void 
UZSolver_calculate_matrix_coeff(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg, Array_Pointers* arr){
    int i = blockIdx.y * blockDim.y + threadIdx.y + 2;
    int j = blockIdx.x * blockDim.x + threadIdx.x + startidx;
    if ((j > endidx) || (i > arg->N)) return;
    
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;

    int seg_no = locate_segment_u(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->dauj, arg->cuoij, arg->mocj, arg->h, NANGDAY);
    
    // printf("ubt = %llf, ubp = %llf\n", arg->ubt[j], arg->ubp[j] );

    _calculate_matrix_coeff(true, j, i, arg->N + 2, 2 * (arg->N) + 1, first, last, seg_no, bienran1, bienran2,
                          arg->ubt[j], arg->ubp[j], arg->t_z[last * (arg->M + 3) + j], arg->t_z[first * (arg->M + 3) + j], arg->bienQ[2], arg->bienQ[3], arg->N ,arg,arr);

}




__global__ void solveU(DOUBLE t, int startidx, int endidx, Argument_Pointers* arg, Constant_Coeffs* coeffs){

    // calculate i, j from thread index. Should use function here for general thread assigning patttern.
    // each thread has its own bienran1, bienran2, dau, cuoi that are passed to uSolver, and executes its own version of uSolver

    int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int j = blockIdx.x * blockDim.x + threadIdx.x + startidx;
    int N = arg->N;
    int M = arg->M;


    if ((i > N) || (j > endidx)) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    locate_segment_u(N, M, &bienran1, &bienran2, &first, &last, i, j, arg->dauj, arg->cuoij, arg->mocj, arg->h, coeffs->NANGDAY);
    
    if ( first > last || i < first - 1 || i >= last) return;
    

    uSolver(t, M + 3, N, first, last, i, j, bienran1, bienran2, arg, coeffs);    
}

__global__ void solveV(DOUBLE t, int startidx, int endidx, Argument_Pointers* arg, Constant_Coeffs* coeffs){

        int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
        int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
        int N = arg->N;
        int M = arg->M;
        if ((i > endidx) || (j > M)) return;
        //printf("thread no %d say hello from forth kernel\n", blockIdx.x*blockDim.x + threadIdx.x);
        // leverage memory ultilization by differrent grid shape. inner function remains unchanged
        bool bienran1 = false;
        bool bienran2 = false;
        int first = 0; int last = 0;
        locate_segment_v(N, M, &bienran1, &bienran2, &first, &last, i, j, arg->daui, arg->cuoii, arg->moci, arg->h, coeffs->NANGDAY);
        if ((first >= last) || (j >= last) || (j < first) ) return;
        //printf("dauu, cuoi: %d %d\n", first, last);
        vSolver(t, M +3, first, last, i, j, bienran1, bienran2, arg->VISCOIDX, arg->Tsyw, arg->v, arg->t_v, arg->u, 
            arg->t_u, arg->z, arg->t_z, arg->Ky1, arg->Htdv, arg->H_moi, coeffs);
}

__global__ void update_margin_elem_U(int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg){
    int j = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int N = arg->N;
    int M = arg->M;
    if (j > endidx) return;
    DOUBLE* t_u = arg->t_u;
    DOUBLE* h = arg->h;
    int* dauj = arg->dauj;
    int* cuoij = arg->cuoij;
    int* mocj = arg->mocj;
    
    int first = 0; int last = 0;
    for (int k = 0; k < mocj[j]; k++){
        int width = segment_limit;
        first = dauj[j * width + k];
        last = cuoij[j * width + k];
        bool bienran1 = false;
        bool bienran2 = false;
        width = M + 3;
        if ((first > 2) || ( (first == 2) && (( abs( h[1 * width + j] + h[1 * width + j - 1]) * 0.5 - NANGDAY)  < epsilon) ))
            bienran1 = true;
    
        if ((last < N) || ((last == N) && ( abs((h[N * width + j] + h[N * width + j - 1]) * 0.5 - NANGDAY ) < epsilon ) ))
            bienran2 = true;

        if (bienran1){
            t_u[(first - 1) * width + j] = 0;
        }
        else 
            t_u[(first - 1) * width + j] = 2 * t_u[first * width + j] - t_u[(first + 1) * width + j];

        if(bienran2)
            t_u[last * width + j] = 0;
        else 
            t_u[last * width + j] = 2 * t_u[(last - 1) * width +  j] - t_u[(last - 2) * width +  j];           
    }
}

__global__ void update_margin_elem_V( int startidx, int endidx, DOUBLE NANGDAY, Argument_Pointers* arg){
    int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int M = arg->M;
    if (i > endidx) return;
    DOUBLE* t_v = arg->t_v;
    DOUBLE* h = arg->h;
    int* daui = arg->daui;
    int* cuoii = arg->cuoii;
    int* moci = arg->moci;
    
    int first = 0; int last = 0;
    for (int k = 0; k < moci[i]; k++){
        bool bienran1 = false;
        bool bienran2 = false;
        int width = segment_limit;
        first = daui[i * width + k];
        last = cuoii[i * width + k];
        width = M + 3;
        if ((first > 2) || ((first == 2) && ( abs((h[i * width + first - 1] + h[(i - 1) * width + first - 1]) * 0.5 - NANGDAY) < epsilon))){
           bienran1 = true;
        }


        if ((last < M) || ((last == M) && ( abs((h[i * width +  last] + h[(i - 1) * width + last]) * 0.5 - NANGDAY) < epsilon))) {
            bienran2 = true;
        }

        if (bienran1)
            t_v[i * width +  first - 1] = 0;
        else{
            t_v[i * width +  first - 1] = 2 * t_v[i * width +  first] - t_v[i * width +  first + 1];
        }
        if (bienran2){

            t_v[i * width +  last] = 0;
        }
        else{
            
            t_v[i * width +  last] = 2 * t_v[i * width +  last - 1] - t_v[i * width +  last - 2];
        }         
    }
        
}





