/**
ULSAN NATIONAL INSTIUTE OF SCIENCE AND TECHNOLOGY
Copyright (c) 2019 HVCL lab
Created by Huong Nguyen
**/


#include "support_funcs.h"

// __device__ void calculate_index(int *i, int *j, int M){
//     int thrx = blockIdx.x * blockDim.x + threadIdx.x;
//     int thry = blockIdx.y * blockDim.y + threadIdx.y;

//     int thrnu =  thrx * (blockDim.y * gridDim.y) + thry;
//     *i = thrnu / M;
//     *j = thrnu % M;
// }

// checked Mar-30

__global__ void Onetime_init( Argument_Pointers *arg, Constant_Coeffs* coeffs){
	__shared__ DOUBLE NANGDAY;
	__shared__ DOUBLE g; 
	__shared__ DOUBLE mu_mn;
	NANGDAY = coeffs->NANGDAY;
	g = coeffs->g;
	mu_mn = coeffs->mu_mn;
	//cout << "test" << endl;
	int M = arg->M;
	int N = arg->N;
	int* khouot = arg->khouot;
	DOUBLE* h = arg->h;
	DOUBLE* H_moi = arg->H_moi;
	DOUBLE* Kx1 = arg->Kx1;
	DOUBLE* Ky1 = arg->Ky1;
	DOUBLE* hsnham = arg->hsnham;
	DOUBLE* htaiz = arg->htaiz;
	DOUBLE* htaiz_bd = arg->htaiz_bd;
	int width = M + 3;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if (( i >= N + 3) || (j >= M + 3)) return;
	// ATTENTION
	khouot [i * width] = khouot [j] = 2;
	i++; j++;
	khouot[i * width +  j] = 2;
	
	if ((i > N + 1) || (j > M + 1)) return;
	// khouot
	if ((h[(i - 1) * width + j - 1] + h[(i - 1) * width + j] + h[i * width + j - 1] + h[i * width + j]) * 0.25 > NANGDAY){
		khouot[i * width + j] = 0;
		H_moi[i * width + j] = 0;
		// htaiz[i * width + j];
	}
	// printf("khouot[%d][%d] = %d\n",i, j, khouot[i * width + j] );
	// giatriHtaiZ
	if (i > N || j > M)  return;
	htaiz[i * width + j] = (h[(i - 1) * width + j - 1] + h[(i - 1) * width + j] + h[i *width + j - 1] + h[i * width + j]) * 0.25;
	htaiz_bd[i * width + j] = htaiz[i * width + j];

	// hesok
	if ( h[i * width + j - 1 ] + h[i * width + j] != 0 )
		Kx1[i * width + j] = g * powf( (h[i * width + j - 1] + h[i * width + j]) * 0.5, -2 * mu_mn) * powf( (hsnham[i * width + j - 1] + hsnham[i * width + j] ) * 0.5, 2);

	if (h[(i - 1) * width + j] + h[i * width + j] != 0)
		Ky1[i * width + j] = g * powf((h[(i - 1) * width + j] + h[i * width + j]) * 0.5, -2 * mu_mn) * powf((hsnham[(i - 1) * width + j] + hsnham[i * width + j]) * 0.5, 2);


}
__global__ void update_h_moi(Argument_Pointers* arg){
	int M = arg-> M ;
	int N = arg-> N;
	int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
	int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (( i > N) || (j > M)) return;
	int* khouot = arg-> khouot;
	DOUBLE* H_moi = arg-> H_moi;
	DOUBLE* t_z = arg-> t_z;
	DOUBLE* htaiz = arg-> htaiz;
	int grid_pos = i * (M + 3) + j;
	if (khouot[grid_pos] == 0){
		H_moi[grid_pos] = htaiz[grid_pos] + t_z[grid_pos];
		
	}
}

__global__ void Reset_states_horizontal(Argument_Pointers* arg, Constant_Coeffs* coeffs){
	int M = arg-> M;
	int N = arg-> N;
	DOUBLE* H_moi = arg->H_moi;
	DOUBLE* htaiz = arg-> htaiz;
	int* khouot = arg->khouot;
	DOUBLE* z = arg->z;
	DOUBLE* t_z = arg-> t_z;
	DOUBLE* t_u = arg-> t_u;
	DOUBLE* t_v = arg-> t_v;
	DOUBLE* FS = arg->FS;

	int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
	// int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
	int offset = M + 3;
	if (i > N) return;
	// if (( i > N) || (j > M)) return;
	for (int j = 2; j <= M; j++)
	{
		if (t_z[i * offset + j] > z[i * offset + j]){
			if (khouot[i * offset + j - 1] == 1){
				t_z[i * offset + j - 1] = t_z[i * offset + j];
				H_moi[i * offset + j - 1] = htaiz[i * offset + j - 1] + t_z[i * offset + j];
	            khouot[i * offset + j - 1] = 0;
	            FS[i * offset + j - 1] = FS[i * offset + j];
			}
			if (khouot[i * offset + j + 1] == 1){
				t_z[i * offset + j + 1] = t_z[i * offset + j];
	            H_moi[i * offset + j + 1] = htaiz[i * offset + j + 1] + t_z[i * offset + j];
	            khouot[i * offset + j + 1] = 0;
	            FS[i * offset + j  + 1] = FS[i * offset + j ];
			}
		}
	}

	for (int j = 2; j <= M; j++){
		
		if ((khouot[i * offset + j] == 0) && (H_moi[i * offset + j] <= coeffs->H_TINH) ){
			
			t_u[(i - 1) * offset + j] = 0;
			t_u[i * offset + j] = 0;
			
			t_v[i * offset + j - 1] = 0;
			t_v[i * offset + j] = 0;
			khouot[i * offset + j] = 1;
			FS[i * offset + j ] = 0;
			
		}
	}

}


// doc
__global__ void Reset_states_vertical(Argument_Pointers* arg, Constant_Coeffs* coeffs){

	int M = arg-> M;
	int N = arg-> N;
	DOUBLE* H_moi = arg->H_moi;
	DOUBLE* htaiz = arg-> htaiz;
	int* khouot = arg->khouot;
	DOUBLE* z = arg->z;
	DOUBLE* t_z = arg-> t_z;
	DOUBLE* t_u = arg-> t_u;
	DOUBLE* t_v = arg-> t_v;
	DOUBLE* FS = arg->FS;

	// int i = blockIdx.y * blockDim.y + threadIdx.y + 2;
	int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (j > M) return;
	int offset = M + 3;

	for (int i = 2; i <= N; i++) {
		if (t_z[i * offset + j] > z[i * offset + j]){
			if (khouot[(i - 1) * offset + j] == 1){
				t_z[(i - 1) * offset + j] = t_z[i * offset + j];
				H_moi[(i - 1) * offset + j] = htaiz[(i - 1) * offset + j] + t_z[i * offset + j];
	            khouot[(i - 1) * offset + j] = 0;
	            FS[(i - 1) * offset + j] = FS[i * offset + j];
			}
			if (khouot[(i + 1) * offset + j] == 1){
				t_z[(i + 1) * offset + j] = t_z[i * offset + j];
	            H_moi[(i + 1) * offset + j] = htaiz[(i + 1) * offset + j] + t_z[i * offset + j];
	            khouot[(i + 1) * offset + j] = 0;
	            FS[(i + 1) * offset + j] = FS[i * offset + j];
			}
		}
	}

	for (int i = 2; i <= N; i++) {
		if ((khouot[i * offset + j] == 0) && (H_moi[i * offset + j] <= coeffs->H_TINH)){
			t_u[(i - 1) * offset + j] = 0;
			t_u[i * offset + j] = 0;

			t_v[i * offset + j - 1] = 0;
			t_v[i * offset + j] = 0;
			khouot[i * offset + j] = 1;
			FS[i * offset + j] = 0;
		}
	}
}


__device__ void Interpolate_FS_ng(int location, int offset, int sign, Argument_Pointers* arg, Constant_Coeffs* coeffs){

	int M = arg->M;
	int width = M + 3;
	DOUBLE* t_u = arg->t_u;
	DOUBLE* u = arg->u;
	DOUBLE* FS = arg->FS;
	DOUBLE deltaX;
	DOUBLE dT, dX;
	dT = coeffs->dT;
	dX = coeffs->dX;
	int i, i_deltaX, i_dX;
	for (int j = 1; j < M; j++){
		if (t_u[location  * width + j] > 0){
			deltaX = (u[location * width + j ] + t_u[location * width + j]) * 0.25 * dT;
			i = offset + sign * (deltaX / dX);
			i_deltaX = __double2int_rn(deltaX);
			i_dX = __double2int_rn(dX);
			FS[location * width  + j] = FS[i * width + j] + sign * (FS[(i + 1) * width + j] - FS[i * width + j] )* (deltaX / dX - i_deltaX / i_dX);
		}
	}

}

__device__ void Interpolate_FS_d(int location, int offset, int sign, Argument_Pointers* arg, Constant_Coeffs* coeffs){

	int M = arg->M;
	int width = M + 3;
	DOUBLE* t_v = arg->t_v;
	DOUBLE* v = arg->v;
	DOUBLE* FS = arg->FS;
	DOUBLE deltaY;
	DOUBLE dT, dY;
	dT = coeffs->dT;
	dY = coeffs->dY;
	int j, j_deltaY, j_dY;
	for (int i = 1; i < M; i++){
		if (t_v[i  * width + location] > 0){
			deltaY = (v[i * width + location ] + t_v[i * width + location]) * 0.25 * dT;
			j = offset + sign * (deltaY / dY);
			j_deltaY = __double2int_rn(deltaY);
			j_dY = __double2int_rn(dY);
			FS[i * width  + location] = FS[i * width + j + 1] + (FS[i * width + j + 1] - FS[i * width + j] )* (deltaY / dY - j_deltaY / j_dY);
		}
	}

}
__global__ void Find_Calculation_limits_Horizontal( Argument_Pointers *arg, Constant_Coeffs* coeffs){

	// int thrx = blockIdx.x * blockDim.x + threadIdx.x;
 //    int thry = blockIdx.y * blockDim.y + threadIdx.y;
 //    int i = thrx * (blockDim.y * gridDim.y) + thry + 2;
	int i = blockDim.y * blockIdx.y + threadIdx.y + 2;
	int M = arg->M;
	int N = arg->N;
    if (i > N) return;
    int* khouot = arg->khouot;
    int* moci = arg->moci;
    int* cuoii = arg->cuoii;
    int* daui = arg->daui;
	int number_of_seg = 0;
	int start = 2;
	int end = 0;
	int offset = M + 3;
	if (i == N)
		Interpolate_FS_ng(N, N - 1, -1, arg, coeffs);
	if (i == 2)
		Interpolate_FS_ng(2, 2, 1, arg, coeffs);

	while (start < M){
		//printf("i: %d, start %d \n",i, start );
		if (khouot[i * offset + start] != 0){
			while ((khouot[i * offset + start]) && (start < M)) start++;
		} 
		if (start + 1 == M) start = M;

		if (khouot[i * offset + start] == 0 && start + 1 < M){
			daui[i * segment_limit + number_of_seg] = start;
			// if (threadIdx.x == 0)
				// printf("start: %d, i: %d\n", start, i );
			end  = start;
			while((khouot[i * offset + end] == 0) && (end < M)) end++;

			if ((khouot[i * offset + end] != 0) && (end <= M)){
				cuoii[i * segment_limit + number_of_seg] = end - 1;
				start = end;
				number_of_seg++;
			} else{
				cuoii[i * segment_limit + number_of_seg] = M ;
				//printf("i: %d, cuoii : %d\n", i, M);
				start = M;
				number_of_seg++;
			}
			
		} 

	}
	moci[i] = number_of_seg;

}


__global__ void Find_Calculation_limits_Vertical(Argument_Pointers *arg, Constant_Coeffs* coeffs){

	// int thrx = blockIdx.x * blockDim.x + threadIdx.x;
 //    int thry = blockIdx.y * blockDim.y + threadIdx.y;
 //    int j = thrx * (blockDim.y * gridDim.y) + thry + 2;
	int j = blockDim.y * blockIdx.y + threadIdx.y + 2;
	int M = arg->M;
	int N = arg->N;
    if (j > M) return;
    int* khouot = arg->khouot;
    int* mocj = arg->mocj;
    int* cuoij = arg->cuoij;
    int* dauj = arg->dauj;
	int number_of_seg  = 0;
	int start = 2;
	int end = 0;
	int offset = M + 3;
	
	if (j == M)
		Interpolate_FS_d(M, M - 1, -1, arg, coeffs);
	while (start < N){
		if (khouot[start * offset + j] != 0  ){
			while ((khouot[start * offset + j]) && (start < N)) start++;
		}
		if (start + 1 == N) start = N;
		if (khouot[start * offset + j] == 0 && start + 1 < N){
			dauj[j * segment_limit + number_of_seg] = start;
			end = start;
			while ( (khouot[end * offset + j] == 0) && (end < N) ) {end++;}
			
			if ((khouot[end * offset + j] != 0) && ( end <= N)){
				
				// if (threadIdx.x == 0 && j == 3)
				// 	printf(" khouot[%d %d], %d\n", end, j, khouot[end * offset + j]);
				cuoij[j * segment_limit + number_of_seg] = end - 1;
				number_of_seg++;
				// if (j == 3 && threadIdx.x == 0) 
				// 	printf("%d %d\n",end, number_of_seg );
				start = end;
			} else{
				cuoij[j * segment_limit + number_of_seg] = N;
				start = N;
				// if (j == 3 && threadIdx.x == 0) 
				// 	printf("%d %d\n",end, number_of_seg );
				number_of_seg++;
			}
		}
	}
	mocj[j] = number_of_seg;

}


__global__ void Htuongdoi(Argument_Pointers* arg){

	int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
	int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int M = arg->M;
	int N = arg->N; 
	DOUBLE* Htdu = arg->Htdu;
	DOUBLE* Htdv = arg->Htdv;
	DOUBLE* h = arg->h;
	DOUBLE* z = arg->z;

	if ((i > N) || (j > M)) return;
	int width = M + 3;
    Htdu[i * width + j] = (h[i * width + j - 1] + h[i * width + j] + z[(i + 1) * width + j] + z[i * width + j]) * 0.5;
    Htdv[i * width + j] = (h[(i - 1) * width + j] + h[i * width + j] + z[i * width + j + 1] + z[i * width + j]) * 0.5;
    
}

// __global__ void boundary_up(DOUBLE t, int M, int N, bool* bienQ, DOUBLE* t_z, int* daui, int* cuoii, 
// 	int* dauj, int* cuoij, DOUBLE* vbt, DOUBLE* vbd, DOUBLE* ubt, DOUBLE* ubp ){
__global__ void boundary_up (DOUBLE t, Argument_Pointers* arg, Constant_Coeffs* coeffs){
	DOUBLE *t_z;
	int M,  *dauj, *cuoij, *bienQ;
	DOUBLE dY = coeffs->dY;
	bienQ = arg->bienQ;
	t_z = arg->t_z;
	M = arg->M;
	dauj = arg->dauj;
	cuoij = arg->cuoij;

	int thrx = blockIdx.x * blockDim.x + threadIdx.x;
    int thry = blockIdx.y * blockDim.y + threadIdx.y;
    int i = thrx * (blockDim.y * gridDim.y) + thry + dauj[M * segment_limit];
	int offset = M + 3;
	//printf("i = %d, cuoij: %d\n", i, cuoij[M * offset] );
	if (i > cuoij[M * segment_limit]) return;
	// if (bienQ[0])
	// 	{vbt[i] = 0;
	// 			printf("here\n");
	// 	}
	else{
		t_z[i * offset + M] = 0.01 * cos(2 * PI / 27.750 * t) * cos(2 * (PI / 100) * (100 - dY / 2));
		t_z[i * offset + M + 1] = 0.01 * cos(2 * PI / 27.75 * t)  * cos(2 * (PI / 100) *  (100 + dY / 2));
		//printf("tz[%d, %d] = %.15f\n",i, M, t_z[i * offset + M]);
	}
}

// __global__ void boundary_down(DOUBLE t,  int M, int N, bool* bienQ, DOUBLE* t_z, int* daui, int* cuoii, 
// 	int* dauj, int* cuoij, DOUBLE* vbt, DOUBLE* vbd, DOUBLE* ubt, DOUBLE* ubp ){

__global__ void boundary_down(DOUBLE t, Argument_Pointers* arg, Constant_Coeffs* coeffs){
	DOUBLE *t_z, *vbd;
	int M,  *dauj, *cuoij, *bienQ;
	DOUBLE dY = coeffs->dY;
	bienQ = arg->bienQ;
	t_z = arg->t_z;
	M = arg->M;
	dauj = arg->dauj;
	cuoij = arg->cuoij;
	vbd = arg->vbd;
	
	int thrx = blockIdx.x * blockDim.x + threadIdx.x;
    int thry = blockIdx.y * blockDim.y + threadIdx.y;
    int i = thrx * (blockDim.y * gridDim.y) + thry + dauj[2 * segment_limit];
	int offset = M + 3;
	if (i > cuoij[2 * segment_limit]) return;
	if (bienQ[1])
		vbd[i] = 0;
	else{
		t_z[i * offset + 2] = 0.01 * cos(2 * PI / 27.75 * t ) * cos(2 * (PI / 100) * dY / 2);
        t_z[i * offset + 1] = 0.01 * cos(2 * PI / 27.75 * t ) * cos(2 * (PI / 100) * (-dY) / 2);
        //if (t >= 6.75) printf(" tz[%d, %d] = %.15f\n",i, 2, t_z[i * offset + 2]);
	}
}

// __global__ void boundary_left(DOUBLE t, int M, int N, bool* bienQ, DOUBLE* t_z, int* daui, int* cuoii, 
// 	int* dauj, int* cuoij, DOUBLE* vbt, DOUBLE* vbd, DOUBLE* ubt, DOUBLE* ubp){

__global__ void boundary_left(DOUBLE t, Argument_Pointers* arg, Constant_Coeffs* coeffs){
	DOUBLE *t_z, *ubt;
	int M, *daui, *cuoii;
	int *bienQ;
	DOUBLE dX = coeffs->dX;
	bienQ = arg->bienQ;
	t_z = arg->t_z;
	M = arg->M;
	daui = arg->daui;
	cuoii = arg->cuoii;
	ubt = arg->ubt;

	int thrx = blockIdx.x * blockDim.x + threadIdx.x;
    int thry = blockIdx.y * blockDim.y + threadIdx.y;
    int i = thrx * (blockDim.y * gridDim.y) + thry + daui[2 * segment_limit];
	int offset = M + 3;
	//printf("i: %d\n", cuoii[2 * offset]);
	if (i > cuoii[2 * segment_limit]) return;
	if (bienQ[2])
		ubt[i] = 0;
	else{
		t_z[2 * offset + i] = 0.01 * cos(2 * PI / 27.75 * t ) * cos(2 * (PI / 100) * dX / 2);
        t_z[1 * offset + i] = 0.01 * cos(2 * PI / 27.75 * t ) * cos(2 * (PI / 100) * (- dX) / 2);

	}

}


__global__ void boundary_right(DOUBLE t, Argument_Pointers* arg, Constant_Coeffs* coeffs){
	DOUBLE *t_z, *ubp;
	int M, N, *daui, *cuoii;
	int *bienQ;
	DOUBLE dX = coeffs->dX;
	bienQ = arg->bienQ;
	t_z = arg->t_z;
	M = arg->M;
	N = arg->N;
	daui = arg->daui;
	cuoii = arg->cuoii;
	ubp = arg->ubp;

	int thrx = blockIdx.x * blockDim.x + threadIdx.x;
    int thry = blockIdx.y * blockDim.y + threadIdx.y;
    int i = thrx * (blockDim.y * gridDim.y) + thry + daui[2 * segment_limit];
	int offset = M + 3;
	if (i > cuoii[N * segment_limit]) return;
	if (bienQ[3])
		ubp[i] = 0;
	else{
		t_z[N * offset + i] = 0.01 * cos(2 * PI / 27.75 * t)  * cos(2 * (PI / 100) *  (100 - dX / 2));
        t_z[(N + 1) * offset + i] = 0.01 * cos(2 * PI / 27.75 * t)  * cos(2 * (PI / 100) *  (100 + dX / 2));
        //if (t >= 6.75) printf("tz[%d, %d] = %.15f\n",N, i, t_z[N * offset + i]);
	}
}



__global__ void preprocess_data(Argument_Pointers* arg, Constant_Coeffs* coeffs){

	// hi is a so-called "scratch pad" array that is used to store some temporary values which
	// are later used to calculate boundary values based on the CCHE formula
	// pre-calculate values and stored in hi is simply to reduce redundant calculation every iterations
	DOUBLE* hi = arg->hi;
	DOUBLE* h = arg->h;
	// DOUBLE* hsnham = arg->hsnham;
	int* daui = arg->daui;
	int* cuoii = arg->cuoii;
	int* dauj = arg->dauj;
	int* cuoij = arg->cuoij;
	int* moci = arg->moci;
	int* mocj = arg->mocj;
	int N = arg->N; int M = arg->M;
	int width = M + 3;

	DOUBLE NANGDAY, dX, ros;
	NANGDAY = coeffs ->NANGDAY;
	dX = coeffs->dX;
	ros = coeffs->ros;

	// int* daui, *cuoii, *dauj, *cuoij, *moci, *mocj;
	// h[i, M]  - hi_up
	for (int k = 0; k < mocj[M]; k++){
		DOUBLE sum  = 0;
		for (int i = dauj[M * segment_limit + k]; i <= cuoij[M * segment_limit + k]; i++){
			hi[i] = (h[i * width + M] + h[(i - 1) * width + M]) / 2.0 - NANGDAY;
			sum += powf(hi[i], 5.0/3.0);// / hsnham[i * width + M];
		}

		for (int i = dauj[M * segment_limit + k]; i <= cuoij[M * segment_limit + k]; i++){
			hi[i] = powf(hi[i], 2.0 / 3.0) / (sum * dX ) ;//* hsnham[i * width + M]);
			
		}

	}

	// h[i, 2] - hi_down

	hi += N + 3;

	for (int k = 0; k < mocj[2]; k++){
		DOUBLE sum  = 0;
		// int sum = 0;
		for (int i = dauj[2 * segment_limit + k]; i <= cuoij[2 * segment_limit + k]; i++){
			hi[i] = (h[i * width + 1] + h[(i - 1) * width + 1]) / 2.0 - NANGDAY;

			sum += powf(hi[i], 5.0/3.0);// / hsnham[i * width + 2];
		}

		for (int i = dauj[2 * segment_limit + k]; i <= cuoij[2 * segment_limit + k]; i++){
			hi[i] = powf(hi[i], 2.0 / 3.0) / (sum * dX) ; // * hsnham[i * width + 2]);
			
		}

	}

	// // h[2, i] - hi_left
	hi += N + 3;
	for (int k = 0; k < moci[2]; k++){
		DOUBLE sum = 0;
		// DOUBLE sum  = 0;
		for (int i = daui[2 * segment_limit + k]; i <= cuoii[2 * segment_limit + k]; i++){
			hi[i] = (h[1 * width + i] + h[1 * width + i - 1]) / 2.0 - NANGDAY;
			sum += powf(hi[i], 5.0/3.0);// / hsnham[2 * width + i];
		}

		for (int i = daui[2 * segment_limit + k]; i <= cuoii[2 * segment_limit + k]; i++){
			hi[i] = powf(hi[i], 2.0 / 3.0) / (sum * dX); // * hsnham[2 * width + i]);
			
		}

	}
	// h[N, i] - hi_right

	hi += M + 3;
	for (int k = 0; k < moci[N]; k++){
		DOUBLE sum  = 0;
		for (int i = daui[N * segment_limit + k]; i <= cuoii[N * segment_limit + k]; i++){
			hi[i] = (h[N * width + i] + h[N * width + i - 1]) / 2.0 - NANGDAY;
			sum += powf(hi[i], 5.0/3.0);// / hsnham[N * width + i];
		}
		for (int i = daui[N * segment_limit + k]; i <= cuoii[N * segment_limit + k]; i++){
			hi[i] = powf(hi[i], 2.0 / 3.0) / (sum * dX);// * hsnham[N * width + i]);
			
		}

	}

	DOUBLE * htaiz_bd = arg->htaiz_bd;
	arg->hmax_u = htaiz_bd[2 * width + M];
	arg->hmax_d = htaiz_bd[2 * width + 2];
	for (int i = 3; i < N; i ++){
		if (htaiz_bd[i * width + M] > arg->hmax_u)
			arg->hmax_u = htaiz_bd[i * width + M];
		if (htaiz_bd[i * width + 2] > arg->hmax_d)
			arg->hmax_d = htaiz_bd[i * width + 2];
	}

	arg->hmax_l = htaiz_bd[2 * width + 2];
	arg->hmax_r = htaiz_bd[N * width + 2];
	for (int i = 3; i < M; i ++){
		if (htaiz_bd[2 * width + i] > arg->hmax_l)
			arg->hmax_l = htaiz_bd[2 * width + i];
		if (htaiz_bd[N * width + i] > arg->hmax_r)
			arg->hmax_r = htaiz_bd[N * width + i];
	}
	if (arg->hmax_u != 0)
	arg->hmax_u = 1.0 / (ros * arg->hmax_u);
	if (arg->hmax_d != 0)
	arg->hmax_d = 1.0 / (ros * arg->hmax_d);
	if (arg->hmax_l != 0)
	arg->hmax_l = 1.0 / (ros * arg->hmax_l);
	if (arg->hmax_r != 0)
	arg->hmax_r = 1.0 / (ros * arg->hmax_r);
	// printf("hmax: %d %lf %lf %lf %lf\n", threadIdx.x, arg->hmax_u, arg->hmax_d, arg->hmax_l, arg->hmax_r );


}

__device__ void Boundary_value(bool isU, DOUBLE t, int location, int location_extension, int width, int total_time,
	int* boundary_type, DOUBLE* hi, DOUBLE* boundary_array, DOUBLE* t_z, DOUBLE* boundary_condition, int* moc, int* dau, int* cuoi){

	int t1 = t / 3600;
	// DOUBLE t2 = (t - (3600.0 * (DOUBLE) t1) ) / 3600.0 ;
	DOUBLE t2 = (t - (3600.0f * t1) ) / 3600.0f ;
	
	// locate which segment the thread is in charge of
	int i = blockIdx.y * blockDim.y + threadIdx.y + 2;

	int seg_no = - 1;
	for (int k = 0; k < moc[location]; k++){
		if ((dau[location * segment_limit +  k] <= i) && (i <= cuoi[location * segment_limit + k])) 
			seg_no = k;
		break;
	}

	// if there is no boundary in a certain edge
	if (seg_no == -1) return;
	

	DOUBLE boundary_value = boundary_condition[seg_no * total_time + t1] * (1.0f - t2) + boundary_condition[seg_no * total_time + t1 + 1] * t2;
	
	if (boundary_type[seg_no] == 2){
		boundary_array[i] = boundary_value * hi[i];
		
	} else{
	// if boudary condition is given in Z
		if (isU){
		t_z[location * width + i] = boundary_value;
		t_z[location_extension * width  + i] = boundary_value;
		} else{
		 t_z[i * width + location] = boundary_value;
		 t_z[i * width + location_extension] = boundary_value;
		}
	}
	
}

__device__ void FS_boundary(bool isU, DOUBLE t, int width, int total_time, int location, DOUBLE hmax,  
	int* boundary_type, DOUBLE* htaiz_bd, DOUBLE* FS, DOUBLE* CC, int* moc, int* dau, int*cuoi ){



	int t1 = t / 3600;
	// DOUBLE t2 = (t - (3600.0 * (DOUBLE) t1) ) / 3600.0 ;
	DOUBLE t2 = (t - (3600.0 * t1) ) / 3600.0 ;
	

	// locate segment
	int i, j;
	i = blockIdx.y * blockDim.y + threadIdx.y + 2;

	int seg_no = - 1;
	for (int k = 0; k < moc[location]; k++){
		if ((dau[location * segment_limit +  k] <= i) && (i <= cuoi[location * segment_limit + k])) 
			seg_no = k;
		break;
	}
	if (seg_no == - 1 || boundary_type[seg_no] != 2) return;

	if (isU){
		i = location;
		j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	} else{     
		j = location;
		i = blockIdx.y * blockDim.y + threadIdx.y + 2;

	}

	DOUBLE boundary_value = CC[seg_no * total_time + t1] * (1.0 - t2)  + CC[seg_no * total_time + t1 + 1] * t2; 
	// hmax = 1 / (ros * hmax)
	// FS[i * width + j] = boundary_value * htaiz_bd[i * width + j] * 1 / (ros  * hmax);
	FS[i * width + j] = boundary_value * htaiz_bd[i * width + j] * hmax;
}




__global__ void Update_Boundary_Value(DOUBLE t, int total_time, Argument_Pointers* arg ){
	
	DOUBLE* hi_up, *hi_down, *hi_left, *hi_right;
	__shared__ int M, N;
	int* boundary_type;

	M = arg->M; N = arg->N;


	if ( (blockIdx.y * blockDim.y + threadIdx.y + 2 > M) && (blockDim.y * blockIdx.y + threadIdx.y + 2  > N)) return;
	// up

	// printf("%d %d %d\n", arg->mocj[M], arg->dauj[M * segment_limit], arg->cuoij[M * segment_limit]);	
	hi_up = arg->hi;
	boundary_type = arg->boundary_type;

	// printf("%d %d\n", arg->dauj[120 * segment_limit], arg->cuoij[120 * segment_limit] );
	Boundary_value(false, t, M, M + 1, M + 3, total_time, boundary_type, hi_up, arg->vbt, arg->t_z, arg->bc_up, arg->mocj, arg->dauj,  arg->cuoij);
	FS_boundary(false, t, M + 3, total_time, M, arg->hmax_u, 
		boundary_type, arg->htaiz_bd, arg->FS, arg->CC_u, arg->mocj, arg->dauj,  arg->cuoij);


	// down
	boundary_type += segment_limit;
	hi_down = hi_up + (N + 3);
	Boundary_value(false, t, 2, 1, M + 3, total_time, boundary_type, hi_down, arg->vbd,arg->t_z, arg->bc_down, arg->mocj, arg->dauj,  arg->cuoij);
	FS_boundary(false, t, M + 3, total_time, 2, arg->hmax_d, 
		boundary_type, arg->htaiz_bd, arg->FS, arg->CC_d, arg->mocj, arg->dauj,  arg->cuoij);

	// left
	boundary_type += segment_limit;
	hi_left = hi_down + (N + 3);

	Boundary_value(true, t, 2, 1, M + 3, total_time, boundary_type, hi_left, arg->ubt, arg->t_z, arg->bc_left, arg->moci, arg->daui, arg->cuoii);
	FS_boundary(true, t, M + 3, total_time, 2, arg->hmax_l, 
		boundary_type, arg->htaiz_bd, arg->FS, arg->CC_l, arg->moci, arg->daui,  arg->cuoii);

	// right

	boundary_type += segment_limit;
	hi_right = hi_left + (M + 3);
	Boundary_value(true, t, N, N + 1, M + 3, total_time, boundary_type, hi_right, arg->ubp, arg->t_z, arg->bc_right, arg->moci, arg->daui, arg->cuoii);
	FS_boundary(true, t, M + 3, total_time, N, arg->hmax_r, 
		boundary_type, arg->htaiz_bd, arg->FS, arg->CC_r, arg->moci, arg->daui,  arg->cuoii);
	

}


// __global__ void update_uvz(int M, int N, DOUBLE* u, DOUBLE* v, DOUBLE* z,  DOUBLE* t_u, DOUBLE* t_v, DOUBLE* t_z, DOUBLE* tmp_u, DOUBLE* tmp_v, int kenhhepd=0, int kenhhepng=0){
__global__ void update_uvz(Argument_Pointers* arg, Constant_Coeffs* coeffs){

	int M, N, kenhhepng, kenhhepd;
	DOUBLE *u, *v, *z, *t_u, *t_v, *t_z;

	M = arg->M; 
	N = arg->N;
	u = arg->u;
	t_u = arg->t_u;
	v = arg->v;
	t_v = arg->t_v;
	z = arg->z;
	t_z = arg->t_z;
	kenhhepd = coeffs->kenhhepd;
	kenhhepng = coeffs->kenhhepng;

	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	if ((i > N) || ( j > M)) return;
	int offset = M + 3;

	// if(kenhhepd + kenhhepng != 0)
	// 	printf("kenhhepng = %d kenhhepd = %d\n", kenhhepd, kenhhepng);

	z[i * offset + j] = t_z[i * offset + j];
	u[i * offset + j] = t_u[i * offset + j] * (1 - kenhhepd);
	v[i * offset + j] = t_v[i * offset + j] * (1 - kenhhepng);
}

__device__ void _normalize (DOUBLE coeff, int N, int M, int closenb_dist, int farnb_dist, DOUBLE* tmp_buffer, DOUBLE* val_buff, int* khouot){
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int width = M + 3;
	int grid_pos = i * width + j;
	tmp_buffer[grid_pos] = val_buff[grid_pos];

	if (i > N || j > M || i < 2 || j < 2) return;
	DOUBLE neigbor = 0;
	if (val_buff[grid_pos] != 0){
		int count = 2;
		neigbor = val_buff[grid_pos - closenb_dist] + val_buff[grid_pos + closenb_dist];

		if (khouot[grid_pos - farnb_dist] == 0){ neigbor += val_buff[grid_pos - farnb_dist]; count++;}

		if (khouot[grid_pos + farnb_dist] == 0){ neigbor += val_buff[grid_pos + farnb_dist]; count++;}

		tmp_buffer[grid_pos] = tmp_buffer[grid_pos] * coeff + (1 - coeff) * neigbor / count;
	}

	// this one is for debugging only. need to change afterward.
}

__global__ void update_buffer(bool updateU, Argument_Pointers* arg, Array_Pointers* arr){
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int width = arg->M + 3;
	int grid_pos = i * width + j;
	if (i > arg->N || j > arg->M || i < 2 || j < 2) return;
	DOUBLE* tmp_buffer;
	DOUBLE* val_buffer;
	if (updateU){
		val_buffer = arg->t_u;
		tmp_buffer = arr->AA;
	} else{
		val_buffer = arg->t_v;
		tmp_buffer = arr->BB;
	}
	val_buffer[grid_pos] = tmp_buffer[grid_pos];
}

__global__ void Normalize(bool isU, Argument_Pointers* arg, Array_Pointers* arr, Constant_Coeffs* coeffs){
	if (isU)
		_normalize(coeffs->heso, arg->N, arg->M, arg->M + 3, 1, arr->AA, arg->t_u, arg->khouot);
	else
		_normalize(coeffs->heso, arg->N, arg->M, 1, arg->M + 3, arr->BB, arg->t_v, arg->khouot);

}

__device__ int locate_segment_v(int N, int M, bool* bienran1, bool* bienran2, int* first, int* last, 
	int row, int col,  int* daui, int* cuoii, int* moci, DOUBLE* h, DOUBLE NANGDAY){
    
    for (int k = 0; k < moci[row]; k++){
        int width = segment_limit;
        if ((daui[row * width +  k] <= col) && (col <= cuoii[row * width + k])) 
        {
            *first = daui[row * width + k];
            *last = cuoii[row * width + k];
            //printf("thread: %d A: dau: %d, cuoi: %d\n", threadIdx.x, *first, *last);
            //printf("first %d\n", *first);
            
            width = M + 3;
            if ((*first > 2) || ((*first == 2) && ((h[row * width + *first - 1] + h[(row - 1) * width + *first - 1]) * 0.5 == NANGDAY))) 
               *bienran1 = true;
            if ((*last < M) || ( (*last == M) && ((h[row * width +  *last] + h[(row - 1) * width + *last]) * 0.5 == NANGDAY) ) )
               *bienran2 = true;
            return k;
        }
    }
}

__device__ int locate_segment_u(int N, int M, bool* bienran1, bool* bienran2, int* first, int* last, 
			int row, int col,  int* dauj, int* cuoij, int* mocj, DOUBLE* h, DOUBLE NANGDAY){
    
    // printf("mocj[%d] = %d\n",col, mocj[col]);
    for (int k = 0; k < mocj[col]; k++){
        int width = segment_limit;
        if ((dauj[col * width +  k] <= row) && (row <= cuoij[col * width + k])) 
        {
            *first = dauj[col * width +  k];
            *last = cuoij[col * width + k];
    
            width = M + 3;
    
            if ((*first > 2) || ( (*first == 2) && ((h[1 * width + col] + h[1 * width + col - 1]) * 0.5 == NANGDAY )) ){
                *bienran1 = true;
                
            }
    
            if ((*last < N) || ((*last == N) && ((h[N * width + col] + h[N * width + col - 1]) * 0.5 == NANGDAY)))
                *bienran2 = true;


        	return k;
        }
    }

}
