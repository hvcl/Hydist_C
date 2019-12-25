/**
ULSAN NATIONAL INSTIUTE OF SCIENCE AND TECHNOLOGY
Copyright (c) 2019 HVCL lab
Created by Huong Nguyen
**/

#include "sediment_transport.h"
#include "support_funcs.h"

// Source function for cohesive case
__device__ DOUBLE source_chs(Constant_Coeffs* coeffs, int i, int j, DOUBLE wsm , DOUBLE Zf, DOUBLE fs, DOUBLE fw, DOUBLE vth, DOUBLE h_moi){
	DOUBLE toee = coeffs->Toe;
	DOUBLE todd = coeffs->Tod;
	DOUBLE S = 0;
	DOUBLE Pd, beta, Cb;
	// for Tan Chau
	DOUBLE tob = coeffs->ro * fw * (vth * vth);
	// if (vth != 0)
	// 	printf("i= %d, j = %d, fw %.15lf vth %.15lf %.15lf\n", i, j, fw, vth, ro);

	if (h_moi > coeffs->ghtoe)
	    toee = toee * (1 + coeffs->hstoe * ((h_moi - coeffs->ghtoe)));
	if (tob > toee)
		S = coeffs->Mbochat / coeffs->ros * (tob - toee) / toee;
	else{
		if (tob < todd){
			// this can be optimized
			Pd = 1 - (tob / todd);
	        beta = 1 + (Zf / (1.25 + 4.75 * pow(Pd , 2.5)));
	        Cb = beta * fs;
	        S = - wsm * Cb * Pd;
	   //      if (abs(S) != 0)
				// printf("saiiiiii, %d %d %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n",i, j, S, fs, Pd, beta, Cb, Zf); 
	         
		}
	}
	return S;
}

// Source function for non cohesive case
__device__ DOUBLE source_nchs(Constant_Coeffs* coeffs, DOUBLE wsm, DOUBLE Zf, DOUBLE Dxr, DOUBLE Ufr, DOUBLE Uf, DOUBLE h_moi, DOUBLE fs){
	DOUBLE S = 0;
	DOUBLE gamac = 0.434 + 5.975 * Zf + 2.888 * Zf * Zf;
    DOUBLE Tx = max(0.0, (Uf * Uf - Ufr * Ufr) / (Ufr * Ufr));

	DOUBLE Cac = 0.015 * coeffs->dm * powf(Tx , 1.5) * powf(Dxr , (-0.3)) / (0.05 * h_moi);
    S = wsm * (Cac - gamac * fs);
	return S;
}

__device__ void _FSi_calculate__mactrix_coeff(Constant_Coeffs* coeffs, DOUBLE t, DOUBLE s_start, bool ketdinh, int i, int j, int first, int last, int seg_no, bool bienran1, bool bienran2, Argument_Pointers* arg, Array_Pointers* arr){

	__shared__ DOUBLE* FS, *H_moi, *t_u, *t_v, *VTH, *Htdu, *Htdv, *Fw, *Kx, *Ky;
	__shared__ int width, tridiag_coeff_width;
	__shared__ DOUBLE Ks, g, wss, Dxr, Ufr, dY, dX, dT, H_TINH; 
	DOUBLE *AA, *BB, *CC, *DD;
	// int sn = last - first - 2;
	width = arg-> M + 3;
	tridiag_coeff_width  = arg->M + 3;
	FS = arg->FS;
	H_moi = arg->H_moi;
	t_v = arg->t_v;
	t_u = arg->t_u;
	Htdv = arg->Htdv;
	Htdu = arg->Htdu;
	VTH = arg->VTH;
	Kx = arg->Kx;
	Ky = arg->Ky;
	Fw = arg->Fw;
	AA = &(arr->AA[i * tridiag_coeff_width]);
    BB = &(arr->BB[i * tridiag_coeff_width]);
    CC = &(arr->CC[i * tridiag_coeff_width]);
    DD = &(arr->DD[i * tridiag_coeff_width]);

    // coeffs:
    Ks = coeffs->Ks;
    g = coeffs->g;
    wss = coeffs->wss;
    Dxr = coeffs->Dxr;
    Ufr = coeffs->Ufr;
    dY = coeffs->dY;
    dX = coeffs->dX;
    dT = coeffs->dT;
    H_TINH = coeffs->H_TINH;

	if (last < first + 1 || j > last - 1 || j < first + 1) return;
	int pos = i * width + j;
	// int offset = first + 1;

	int offset = first + 1;

	DOUBLE c = 18 * log(12 * H_moi[pos] / Ks);
	DOUBLE Uf = sqrt(g) * abs(VTH[pos]) / c;
	DOUBLE wsm = wss *  powf((1 - FS[pos]), 4);
	DOUBLE Zf = 0;
	
	if (FS[pos] > 0.00001)
		Zf = wsm / (0.4 * (Uf + 2 * wsm));

	DOUBLE S;
	if (ketdinh)
		S = source_chs(coeffs, i, j, wsm, Zf, FS[pos], Fw[pos], VTH[pos], H_moi[pos]);
	else
		S = source_nchs(coeffs, wsm, Zf, Dxr, Ufr, Uf, H_moi[pos], FS[pos]);
	// if (S != 0)
	// 	printf("saiiiiii, %.15lf\n", S);
	if ( t < s_start) 
		S = 0;
	DOUBLE gamav = 0.98 - 0.198 * Zf + 0.032 * Zf * Zf;

	AA[j] = -gamav * 0.5 * (t_v[pos - 1] + t_v[pos]) / (dY * 2) - Htdv[pos - 1] * Ky[pos - 1] / (H_moi[pos] * (dY * dY));
	CC[j] = gamav * 0.5 * (t_v[pos - 1] + t_v[pos]) / (dY * 2) - Htdv[pos] * Ky[pos] / (H_moi[pos] * (dY * dY));
	BB[j] = 2.0 / dT + (Htdv[pos] * Ky[pos] + Htdv[pos - 1] * Ky[pos - 1]) / (H_moi[pos] * (dY * dY));

	DOUBLE p, q;
	p = 0;
	q = 0;
	if ((H_moi[pos + 2 * width] > H_TINH) &&  (H_moi[pos - width] > H_TINH))
	    p = (FS[pos + width] - FS[pos]);

	if  (H_moi[pos- 2 * width] > H_TINH && H_moi[pos + width] > H_TINH)
	    q = (FS[pos] - FS[pos - width]);

	// if (i == 463 && j == last - 1)
	// 		printf("trc %.15lf %.15lf %d %d\n", p,q, H_moi[pos - 2 * width] > H_TINH, H_moi[pos + width] > H_TINH);

	p = (1 / (H_moi[pos] * (dX * dX))) * (Htdu[pos] * Kx[pos] * p - Htdu[pos - width] * Kx[pos - width] * q);
	


	if (H_moi[pos + 2 * width] > H_TINH && H_moi[pos - 2 * width] > H_TINH){
	    q = (FS[pos + width] - FS[pos - width]) / (dX * 2);
		q = (t_u[pos] + t_u[pos - width]) * 0.5  * q * gamav;
	}
	else q = 0;


	DD[j] = FS[pos] / (dT * 0.5) - q + p + (S / H_moi[pos]);


	if (j == first + 1){
		if ((bienran1) || (t_v[i * width + first] == 0) ){
			BB[j] = BB[j] + AA[j];
				// printf("last row %d %d %.15lf\n",463, first + 1, BB[j]);

		} else{
			DD[j] = DD[j] - AA[j] * FS[i * width + first];
		}
	}
		// this is likely to be changed
	if (j == last - 1){

		if ((bienran2) || (t_v[i * width + last] == 0) ){
			// For Tan Chau
			BB[j] = BB[j] + CC[j];
		} else{
			DD[j] = DD[j] - CC[j] * FS[i * width + last];
		}
		arr->SN[i * segment_limit + seg_no] = last - first - 2;
		int sn = last - first - 2;
		// printf("i = %d sn = %d\n", i , sn );
	
	}
}


__device__ void _FSi_extract_solution(int i, int j, int first, int last, bool bienran1, bool bienran2, DOUBLE NDnen, Argument_Pointers* arg, Array_Pointers* arr){

	__shared__ DOUBLE* FS, *x, *t_v, *tFS;
	__shared__ int width;
	width = arg->M + 3;
	tFS = arg->tFS;
	FS = arg->FS;
	t_v = arg->t_v;
	int pos = i * width + j;
	x = arr->x + i * (arg->M + 3);
	if (j > last - 1 || j < first + 1 || last < first + 1)
		return;

	if (x[j] < 0) 
		x[j] = NDnen;
	tFS[pos] = x[j];

	if (j == first + 1){

		if ((bienran1) || (t_v[i * width + first] == 0))
			tFS[i * width + first] = tFS[pos];
		else 
			tFS[i * width + first] = FS[i * width + first];
	}
	if (j == last - 1){

		if ((bienran2) || (t_v[i * width + last] == 0))
			tFS[i * width + last] = tFS[pos];
		else {
			tFS[i * width + last] = FS[i * width + last];
			// printf("this should be here %d %d \n",i, j);
		}
	}

}


__device__ void _FSj_calculate__mactrix_coeff(Constant_Coeffs* coeffs, DOUBLE t, DOUBLE s_start, bool ketdinh, int i, int j, int first, int last, int seg_no, bool bienran1, bool bienran2, Argument_Pointers* arg, Array_Pointers* arr){
	if (i > last - 1 || i < first + 1 || last < first + 1)
		return;
	__shared__ DOUBLE* FS, *H_moi, *t_u, *t_v, *VTH, *Htdu, *Htdv, *Fw, *Kx, *Ky;
	__shared__ int width, tridiag_coeff_width;
	__shared__ DOUBLE Ks, g, wss, Dxr, Ufr, dY, dX, dT, H_TINH; 

	DOUBLE *AA, *BB, *CC, *DD;

	width = arg-> M + 3;
	tridiag_coeff_width = arg->N + 3;
	FS = arg->FS;
	H_moi = arg->H_moi;
	t_v = arg->t_v;
	t_u = arg->t_u;
	Htdv = arg->Htdv;
	Htdu = arg->Htdu;
	VTH = arg->VTH;
	Kx = arg->Kx;
	Ky = arg->Ky;
	Fw = arg->Fw;
	AA = &(arr->AA[j * tridiag_coeff_width]);
    BB = &(arr->BB[j * tridiag_coeff_width]);
    CC = &(arr->CC[j * tridiag_coeff_width]);
    DD = &(arr->DD[j * tridiag_coeff_width]);
	int pos = i * width + j;

	Ks = coeffs->Ks;
    g = coeffs->g;
    wss = coeffs->wss;
    Dxr = coeffs->Dxr;
    Ufr = coeffs->Ufr;
    dY = coeffs->dY;
    dX = coeffs->dX;
    dT = coeffs->dT;
    H_TINH = coeffs->H_TINH;

	DOUBLE c = 18 * log(12 * H_moi[pos] / Ks);
	DOUBLE Uf = sqrt(g) * abs(VTH[pos]) / c;
	DOUBLE wsm = wss *  pow((1 - FS[pos]), 4);
	DOUBLE Zf = 0;
	DOUBLE S = 0;

	DOUBLE gamav = 0.98 - 0.198 * Zf + 0.032 * Zf * Zf;

	if (FS[pos] > 0.00001)
		Zf = wsm / (0.4 * (Uf + 2 * wsm));
	// if (j == 3){
	// 	printf("tai j = 2 FS[%d %d] %.15lf\n", i, j, FS[i * width + j]);
	// }
	// if (FS[pos] != 0)
	// 	printf("%d %d %.15lf\n", i, j, FS[pos]);
	if (ketdinh){
		S = source_chs(coeffs, i, j, wsm, Zf, FS[pos], Fw[pos], VTH[pos], H_moi[pos] );
		// if (i == 463 && (j == 173 || j == 172))
  //   	printf("S = %d %.15lf\n",j, S);
	}
	else
		S = source_nchs(coeffs, wsm, Zf, Dxr, Ufr, Uf, H_moi[pos], FS[pos]);        
    
    if ( t < s_start) 
    	S = 0;
    
    AA[i] = -gamav * 0.5 * (t_u[pos] + t_u[pos - width]) / (dX * 2) - Htdu[pos - width] * Kx[pos - width] / (H_moi[pos] * (dX * dX));
    CC[i] = gamav * 0.5 * (t_u[pos] + t_u[pos - width]) / (dX * 2) - Htdu[pos] * Kx[pos] / (H_moi[pos] * (dX * dX));
    BB[i] = 2.0 / dT + (Htdu[pos] * Kx[pos] + Htdu[pos - width] * Kx[pos - width]) / (H_moi[pos] * (dX * dX));

    DOUBLE p, q;
    p = q = 0;
    if ((H_moi[pos - 1] > H_TINH) && (H_moi[pos + 2] > H_TINH))
        p = (FS[pos + 1] - FS[pos]);
    if ((H_moi[pos - 2] > H_TINH) && (H_moi[pos + 1] > H_TINH))
        q = (FS[pos] - FS[pos - 1]);

    p = (1 / (H_moi[pos] * (dY * dY))) * (Htdv[pos] * Ky[pos] * p - Htdv[pos - 1] * Ky[pos - 1] * q);

    if ((H_moi[pos - 2] > H_TINH) && (H_moi[pos + 2] > H_TINH)){
        q = (FS[pos + 1] - FS[pos - 1]) / (dY * 2);
        q = (t_v[pos] + t_v[pos - 1]) * 0.5 * q * gamav;
    }
    else q = 0;

    

    DD[i] = FS[pos] / (dT * 0.5) - q + p + (S / H_moi[pos]);

    // if (DD[i] != 0) 
    // 	printf("i = %d j = %d, %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %d\n", i, j, AA[i], BB[i], CC[i], DD[i],p, q, S, ketdinh);

    if (i == first + 1){

	    if ((bienran1) || (t_u[first * width + j] == 0))
	    	BB[i] += AA[i];
	    else
	    	DD[i] -= AA[i] * FS[first * width + j];

	}   

	if (i == last - 1){

		if ((bienran2) && (t_u[last * width + j] == 0)){
	    	BB[i] += CC[i];
	    } else
	    	DD[i] -= CC[i] * FS[last * width + j];
		arr->SN[j * segment_limit + seg_no] = last - first - 2;

	}


}


__device__ void _FSj_extract_solution(int i, int j, int first, int last, bool bienran1, bool bienran2, DOUBLE NDnen, Argument_Pointers* arg, Array_Pointers* arr){

	__shared__ DOUBLE* FS, *x, *t_u, *tFS;
	__shared__ int width;
	width = arg->M + 3;
	FS = arg->FS;
	tFS = arg->tFS;
	t_u = arg->t_u;
	int pos = i * width + j;
	if (i > last - 1 || i < first + 1 || last < first + 1)
		return;
	x = arr->x + j * (arg->N + 3);
	tFS[pos] = x[i];
	if (x[i] < 0) 
		tFS[pos] = NDnen;

	if (i == first + 1){
	    if ((bienran1) || (t_u[pos - width] == 0))
	    	// this is for song luy only
	    	tFS[first * width + j] = tFS[pos];
	    else
	    	tFS[first * width + j] = FS[first * width + j];
	} 

	if (i == last - 1){
		if ((bienran2) || (t_u[pos + width] == 0)){

			tFS[last * width + j] = tFS[pos];
	    } else
	    // this case means FS [i, last] is boundary condition
	    	tFS[last * width + j] = FS[last * width + j];
	}

}


__global__ void Calculate_Qb(bool ketdinh, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
	int i = blockIdx.y* blockDim.y + threadIdx.y + 2;
	int j = blockIdx. x* blockDim.x + threadIdx.x + 2;
	if (i > arg->N || j > arg->M || i < 2 || j < 2)
		return;
	__shared__  DOUBLE Toe, ro, ghtoe, hstoe, Sx, g, dm, Dxr; 
	DOUBLE Tob;
	DOUBLE Toee = Toe;
	DOUBLE Tx = 0;
	__shared__ int* khouot;
	__shared__ DOUBLE *H_moi, *t_u, *t_v, *VTH, *Qbx, *Qby, *Fw;
	__shared__ int width, M, N;
	M = arg->M;
	N = arg->N;
	width = M + 3;
	H_moi = arg->H_moi;
	Fw = arg->Fw;
	Qbx = arg->Qbx;
	Qby = arg->Qby;
	t_v = arg->t_v;
	t_u = arg->t_u;
	VTH = arg->VTH;	
	khouot = arg->khouot;
	Toe = coeffs->Toe;
	ro = coeffs->ro;
	ghtoe = coeffs->ghtoe;
	hstoe = coeffs->hstoe;
	Sx = coeffs->Sx;
	g = coeffs->g;
	dm = coeffs->dm;
	Dxr = coeffs->Dxr;

	int pos = i * width + j;
	if ((!VTH[pos]) && (khouot[pos] == 0)){
      
	    Tob = ro * Fw[pos] * (VTH[pos] * VTH[pos]);
	    if (H_moi[pos] > ghtoe)
	        Toee = Toe * (1 + hstoe * ((H_moi[pos] - ghtoe)));
	    
	    if (Tob > Toee)
	        Tx = Tob / Toee - 1;
    	Qbx[pos] = 0.053 * powf( (Sx - 1) * g , 0.5) * powf(dm , 1.5) * powf(Tx , 2.1) * powf(Dxr , -0.3) * (t_u[pos] + t_u[pos - width]) * 0.5 / VTH[pos];
    	Qby[pos] = 0.053 * powf( (Sx - 1) * g , 0.5) * powf(dm , 1.5) * powf(Tx , 2.1) * powf(Dxr , -0.3) * (t_v[pos] + t_v[pos - 1]) * 0.5 / VTH[pos];
    }
    else {
    	Qbx[pos] = 0;
    	Qby[pos] = 0;
    }
    if (i == 2){
    	Qbx[pos - width] = Qbx[pos];
    	Qby[pos - width] = Qby[pos];
    	
    }
    if (i == N){
    	Qbx[pos + width] = Qbx[pos];
    	Qby[pos + width] = Qby[pos];
    	
    }
    if (j == 2){
    	Qbx[pos - 1] = Qbx[pos];
    	Qby[pos - 1] = Qby[pos];
    	if (i == 2)
    		Qbx[1 * width + 1] = Qbx[pos];
    		Qby[1 * width + 1] = Qby[pos];
    	if (i == N)
    		Qbx[(N + 1 ) * width + 1] = Qbx[pos];
    		Qby[(N + 1 ) * width + 1] = Qby[pos];
    }
    if (j == M){
    	Qbx[pos + 1] = Qbx[pos];
    	Qby[pos + 1] = Qby[pos];
    	if (i == 2){
    		Qbx[1 * width + M + 1] = Qbx[pos];
    		Qby[1 * width + M + 1] = Qby[pos];
    	}
    	if (i == N)
    	{
    		Qbx[(N + 1) * width + M + 1] = Qbx[pos];
    		Qby[(N + 1) * width + M + 1] = Qby[pos];
    	}
    }

}

__device__ void _bed_load(DOUBLE t, bool ketdinh, int i, int j, int first, int last, bool bienran1, bool bienran2, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
	if (last - first < 2 || j < first || j > last)
		return;
	__shared__ DOUBLE* FS, *H_moi, *t_u, *t_v, *VTH, *Htdu, *Htdv, *Fw, *Kx, *Ky, *Qbx, *Qby, *dH ;
	__shared__ int width, *khouot, M, N;
	__shared__ DOUBLE Ks, g, wss, Dxr, Ufr, dX, dY, dT, Dorong;
	M = arg->M;
	N = arg->N;
	width = arg-> M + 3;
	FS = arg->FS;
	H_moi = arg->H_moi;
	t_u = arg->t_u;
	t_v = arg->t_v;
	Htdv = arg->Htdv;
	Htdu = arg->Htdu;
	VTH = arg->VTH;
	Kx = arg->Kx;
	Ky = arg->Ky;
	Fw = arg->Fw;
	Qbx = arg->Qbx;
	Qby = arg->Qby;
	khouot = arg->khouot;
	dH = arg->dH;

	Ks = coeffs->Ks;
    g = coeffs->g;
    wss = coeffs->wss;
    Dxr = coeffs->Dxr;
    Ufr = coeffs->Ufr;
    dY = coeffs->dY;
    dX = coeffs->dX;
    dT = coeffs->dT;
    Dorong = coeffs->Dorong;

	int pos = i * width + j;
	DOUBLE p, q;
	p = q = 0;
	if (khouot[pos - width] == 1){
		if (t_u[pos] < 0)
            p = (-3 * Qbx[pos] + 4 * Qbx[pos + width] - Qbx[pos + 2 * width]) / (dX * 2) ;

	} else {
		if (khouot[pos + width] == 1){
			if (t_u[pos > 0])
                p = (3 * Qbx[pos] - 4 * Qbx[pos - width] + Qbx[pos - 2 * width]) / (dX * 2);
		} else 
            p = (Qbx[pos + width] - Qbx[pos - width]) / (dX * 2);
	}
    
    if (khouot[pos - 1] == 1) {
        if (t_v[pos] < 0) 
            q = (-3 * Qby[pos] + 4 * Qby[pos + 1] - Qby[pos + 2]) / (dY * 2);
    } else {
        if (khouot[pos + 1] == 1) {
            if (t_v[pos] > 0) 
                q = (3 * Qby[pos] - 4 * Qby[pos - 1] + Qby[pos - 2]) / (dY * 2);
        }  else
            q = (Qby[pos + 1] - Qby[pos - 1]) / (dY * 2);
    }
        
    DOUBLE tH = p + q;
    p = FS[pos + width] - FS[pos];

    if ((khouot[pos + 2 * width] == 1 || khouot[pos - width] == 1) && (i < N)) 
        p = 0;
    q = FS[pos] - FS[pos - width];
    if ((i > 2) && (khouot[pos - 2 * width] == 1 || khouot[pos + width] == 1)) 
        q = 0;
        
    p = 1 / (dX * dX) * (Htdu[pos] * Kx[pos] * p - Htdu[pos - width] * Kx[pos - width] * q);
        
    tH = tH + p;
	p = FS[pos + 1] - FS[pos];

	if (((last < M) && (khouot[pos - 1] == 1 || khouot[pos + 2] == 1)) ||
		(last == M && (khouot[pos - 1] == 1)))
			p = 0;
    q = FS[pos] - FS[pos - 1];
    if (((first > 2) && (khouot[pos - 2] == 1 || khouot[pos + 1] == 1)) || (first == 2 && khouot[pos + 1] == 1))
    	q = 0;
                 
    p = 1 / ((dY * dY)) * (Htdv[pos] * Ky[pos] * p - Htdv[pos - 1] * Ky[pos - 1] * q);
    
    tH = tH + p;

    
    DOUBLE c = 18 * log(12 * H_moi[pos] / Ks);
	DOUBLE Uf = sqrt(g) * abs(VTH[pos]) / c;
	DOUBLE wsm = wss * pow((1 - FS[pos]), 4);
	DOUBLE Zf = 0;
	DOUBLE S = 0;
	if (FS[pos] > 0.0)
		Zf = wsm / (0.4 * (Uf + 2 * wsm));
	
	if (ketdinh)
		S = source_chs(coeffs, i, j, wsm, Zf, FS[pos], Fw[pos], VTH[pos], H_moi[pos] );
	else
		S = source_nchs(coeffs, wsm, Zf, Dxr, Ufr, Uf, H_moi[pos], FS[pos]); 

    dH[pos] += dT / (1 - Dorong) * (S + tH);   
    
}

__global__ void hesoK(Constant_Coeffs* coeffs, Argument_Pointers* arg){
	
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x  + threadIdx.x;
	int M, N;
	M = arg->M;
	N = arg->N;
	if (i > N || j > M)
		return;
	DOUBLE Cz = 0;
	__shared__ DOUBLE *Kx, *Ky, *Htdu, *Htdv, *t_u, *t_v, *h;
	__shared__ DOUBLE Ks, g;
	Kx = arg->Kx;
	Ky = arg->Ky;
	Htdu = arg->Htdu;
	Htdv = arg->Htdv;
	t_u = arg->t_u;
	t_v = arg->t_v;
	h = arg->h;
	Ks = coeffs->Ks;
	g = coeffs->g;
	__shared__ int width;
	width = M + 3;
	int pos = i * width + j;
	if (Htdu[pos] > 0){
		Cz = 18 * log(12 * Htdu[pos] / Ks);
		Kx[pos] = 5.93 * sqrt(g) * (h[pos - 1] + h[pos])* 0.25 * abs(t_u[pos]) / Cz;
		Kx[pos] = min(100.0, max(5.0, Kx[pos]));
	}
	if (Htdv[pos] > 0){
		Cz = 18 * log(12 * Htdv[pos] / Ks);
		Ky[pos] = 5.93 * sqrt(g) * (h[pos - width] + h[pos])* 0.25 * abs(t_v[pos]) / Cz;
		Ky[pos] = min(100.0, max(5.0, Ky[pos]));
	}
}

__global__ void Find_VTH(Constant_Coeffs* coeffs, Argument_Pointers* arg){
	// printf("%d %d\n", blockDim.x, blockDim.y );
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x  + threadIdx.x;
	int M, N;
	M = arg->M;
	N = arg->N;
	if ( i < 2 || j < 2 || j > M || i > N)
		return;
	__shared__ DOUBLE * H_moi, *htaiz, *VTH, *t_u, *t_v;
	__shared__ int width;
	width = M + 3;
	H_moi = arg->H_moi;
	htaiz = arg->htaiz;
	t_u = arg->t_u;
	t_v = arg->t_v;
	VTH = arg->VTH;
	int pos = i * width + j;
	// if (i == 463 && j == 5)
	// 	printf("htaiz = %.15lf, h_moi %.15lf\n", htaiz[pos], H_moi[pos] );
	if ((htaiz[pos] > coeffs->NANGDAY) && (H_moi[pos] > coeffs->H_TINH)){
		DOUBLE ut = (t_u[pos] + t_u[pos - width]) * 0.5;
		DOUBLE vt = (t_v[pos] + t_v[pos - 1]) * 0.5;
		VTH[pos] = sqrt(ut * ut + vt * vt);
		// if (i == 463 && j == 5)
		// 	printf("ut = %.15lf vt = %.15lf\n", ut, vt);
	} else VTH[pos] = 0;
}

__global__ void Scan_FSi(DOUBLE t, DOUBLE s_start, bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
// find first, last, bienran1, bienran2, i, j, pass argument
	int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if (i > endidx ) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    int seg_no = locate_segment_v(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->daui, arg->cuoii, arg->moci, arg->h, coeffs->NANGDAY);
    _FSi_calculate__mactrix_coeff(coeffs, t, s_start, ketdinh, i, j, first, last, seg_no, bienran1, bienran2, arg, arr);
}
__global__ void FSi_extract_solution( bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs*coeffs){
	int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if (i > endidx ) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    locate_segment_v(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->daui, arg->cuoii, arg->moci, arg->h, coeffs->NANGDAY);
    _FSi_extract_solution(i, j, first, last, bienran1, bienran2, coeffs->NDnen, arg, arr);

}


__global__ void Scan_FSj(DOUBLE t, DOUBLE s_start, bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
// find first, last, bienran1, bienran2, i, j, pass argument
	int i = blockIdx.y * blockDim.y + threadIdx.y + 2;
    int j = blockIdx.x * blockDim.x + threadIdx.x + startidx;
    // if (i == 2 && j == startidx)
    //     printf("add0 %x\n", arg->Htdu);
    if (j > endidx ) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    int seg_no = locate_segment_u(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->dauj, arg->cuoij, arg->mocj, arg->h, coeffs->NANGDAY);
    _FSj_calculate__mactrix_coeff(coeffs, t, s_start, ketdinh, i, j, first, last, seg_no, bienran1, bienran2, arg, arr);
}

__global__ void FSj_extract_solution(bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
	int i = blockIdx.y * blockDim.y + threadIdx.y + 2;
    int j = blockIdx.x * blockDim.x + threadIdx.x + startidx;
    if (j > endidx ) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    locate_segment_u(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->dauj, arg->cuoij, arg->mocj, arg->h, coeffs->NANGDAY);
    _FSj_extract_solution(i, j, first, last, bienran1, bienran2,coeffs->NDnen, arg, arr);
}


__global__ void Update_FS(Argument_Pointers* arg){
	int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
	int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (i  > arg->N || j > arg->M) return;
	DOUBLE* tFS, *FS;
	FS = arg->FS;
	tFS = arg->tFS;
	int width = arg->M + 3;
	FS[i * width + j] = tFS[i * width + j];
}
__global__ void BedLoad(DOUBLE t, bool ketdinh, int startidx, int endidx, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
	int i = blockIdx.y * blockDim.y + threadIdx.y + startidx;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if (i >= endidx ) return;
    bool bienran1 = false;
    bool bienran2 = false;
    int first = 0; int last = 0;
    locate_segment_v(arg->N, arg->M, &bienran1, &bienran2, &first, &last, i, j, arg->daui, arg->cuoii, arg->moci, arg->h, coeffs->NANGDAY);
    _bed_load(t, ketdinh, i, j, first, last, bienran1, bienran2, arg, arr, coeffs);
}