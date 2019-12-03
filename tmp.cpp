
#include "Engine.h"

template <typename T>
pair<int, int> load_file(string filename, vector<T> arr){

	ifstream ifs;
	ifs.open(filename);
	if (!ifs) {
		cout << "cannot open" <<  filename << endl;
		exit(1);
	}
	string line;
	int M, N;
	N = 0;
	while (getline(ifs, line)){
		N++;
		stringstream linestream(line);
		T x;
		M = 0;
		while (linestream >> x){
			M++;
			arr.push_back(x);
		}
		arr.push_back(arr[arr.size() - 1]); // pad edge
	}
	return pair(N, M)
}

pair<int, int> load_input(string dir,
				vector<DOUBLE> &h,
				vector<DOUBLE> &hsnham,
				vector<DOUBLE> &VISCOINDX, 
				vector<DOUBLE> &bc_up,
				vector<DOUBLE> &bc_right,
				vector<DOUBLE> &bc_left, 
				vector<DOUBLE> &bc_down, 
				vector<DOUBLE> &CC_u, 
				vector<DOUBLE> &CC_d,
				vector<DOUBLE> &CC_l,
				vector<DOUBLE> &CC_r,
				vector<int> &bienQ)
{
	pair<int, int > size;
	size = load_file<DOUBLE> (dir + "bandodosau.txt", h);
	load_file<DOUBLE> (dir + "hsnham.txt", hsnham);
	load_file<DOUBLE> (dir + "hsnhotroiA.txt", VISCOINDX);
	load_file<DOUBLE> (dir + "bientren.txt", bc_up);
	load_file<DOUBLE> (dir + "bienduoi.txt", bc_down);
	load_file<DOUBLE> (dir + "bientrai.txt", bc_left);
	load_file<DOUBLE> (dir + "bienphai.txt", bc_right);
	load_file<DOUBLE> (dir + "FSbientren.txt", CC_u);
	load_file<DOUBLE> (dir + "FSbienduoi.txt", CC_d);
	load_file<DOUBLE> (dir + "FSbientrai.txt", CC_l);
	load_file<DOUBLE> (dir + "FSbienphai.txt", CC_r);
	load_file<int> (dir + "boundary_type.txt", bienQ);
	return size;
}

void load_initial_condition(string dir, 
							vector<DOUBLE> &u, 
							vector<DOUBLE> &v,
							vector<DOUBLE> &z,
							vector<DOUBLE> &FS,
							vector<int> &khouot)
{

	load_file<DOUBLE> (dir + "u.txt", u);
	load_file<DOUBLE> (dir + "v.txt", v);
	load_file<DOUBLE> (dir + "z.txt", z);
	load_file<DOUBLE> (dir + "FS.txt", FS);
	load_file<int> (dir + "khouot.txt", khouot);
}

template <typename T>
T* device_alloc_and_copy(vector<T> &h_array){
	int nBytes = h_array.size() * sizeof(T);
	T* d_array;
	cudaError_t malloc_status = cudaMalloc((void**) &d_array, nBytes);
	// assert success here
	cudaError_t copy_status = cudaMemcpy(T, h_array.data(), nBytes, cudaMemcpyHostToDevice);
	// asser success here
	// returna address of array on device here
	return d_array;

}

template <typename T>
T* device_alloc(int nBytes){
	T* d_array;
	cudaError_t status = cudaMalloc((void**) &d_array, nBytes);
	// assert success
	return d_array;
}

template <typename T>
void device_copy(vector<T> &source, T* des)
{
	cudaError_t status cudaMemcpy(des, source.data(), sizeof(T) * source.size(), cudaMemcpyHostToDevice);
	// assert success here
}

// this function allocate memory on GPU, store those pointers values on a local variables
// then copy those value to device, save arguemnt pointer struct's pointer on device to device_arg_ptr
// and return the argument pointer struct that store address of pointers on device
Argument_Pointers attribute_arrays_memory_alloc(int device, Host_arrays &ap, Argument_Pointers** device_arg_ptr)
{
	// cuda mem alloc corresponding arrays
	int M = ap.M;
	int N = ap.N;
	Argument_Pointers d_ap ;
	d_ap.M = M;
	d_ap.N = N;

	// maps
	d_ap.h = device_alloc_and_copy<DOUBLE> (&ap.h);
	d_ap.hsnham = device_alloc_and_copy<DOUBLE> (&ap.hsnham);
	d_ap.VISCOINDX = device_alloc_and_copy<DOUBLE> (&ap.VISCOINDX);
	d_ap.Fw = device_alloc_and_copy<DOUBLE> (&ap.Fw);

	// boundary condition for FS
	d_ap.CC_u = device_alloc_and_copy<DOUBLE> (&ap.CC_u);
	d_ap.CC_d = device_alloc_and_copy<DOUBLE> (&ap.CC_d);
	d_ap.CC_l = device_alloc_and_copy<DOUBLE> (&ap.CC_l);
	d_ap.CC_r = device_alloc_and_copy<DOUBLE> (&ap.CC_r);

	// boundary conditions
	d_ap.bienQ = device_alloc_and_copy<int> (&ap.bienQ);
	d_ap.boundary_type = device_alloc_and_copy<int> (&ap.boundary_type);

	d_ap.bc_up = device_alloc_and_copy<DOUBLE> (&ap.bc_up);
	d_ap.bc_down = device_alloc_and_copy<DOUBLE> (&ap.bc_down);
	d_ap.bc_left = device_alloc_and_copy<DOUBLE> (&ap.bc_left);
	d_ap.bc_right = device_alloc_and_copy<DOUBLE> (&ap.bc_right);

	// intitial conditions
	d_ap.u = device_alloc_and_copy<DOUBLE> (&ap.u);
	d_ap.v = device_alloc_and_copy<DOUBLE> (&ap.v);
	d_ap.z = device_alloc_and_copy<DOUBLE> (&ap.z);
	d_ap.FS = device_alloc_and_copy<DOUBLE> (&ap.FS);
	// d_ap.khouot = device_alloc_and_copy<int> (&ap.khouot);

	// cuda mem alloc arrays on device only

	int nBytes = ap.u.size() * sizeof (DOUBLE);
	d_ap.t_u = device_alloc<DOUBLE> (nBytes);
	d_ap.t_v = device_alloc<DOUBLE> (nBytes);
	d_ap.t_z = device_alloc<DOUBLE> (nBytes);
	d_ap.Kx1 = device_alloc<DOUBLE> (nBytes);
	d_ap.Ky1 = device_alloc<DOUBLE> (nBytes);
	d_ap.Hdtu = device_alloc<DOUBLE> (nBytes);
	d_ap.Hdtv = device_alloc<DOUBLE> (nBytes);
	d_ap.htaiz = device_alloc<DOUBLE> (nBytes);
	d_ap.htaiz_bd = device_alloc<DOUBLE> (nBytes);
	d_ap.H_moi = device_alloc<DOUBLE> (nBytes);

	d_ap.Tsxw = device_alloc<DOUBLE> (nBytes);
	d_ap.Tsyw = device_alloc<DOUBLE> (nBytes);

	d_ap.VTH = device_alloc<DOUBLE> (nBytes);
	d_ap.Qbx = device_alloc<DOUBLE> (nBytes);
	d_ap.Qby = device_alloc<DOUBLE> (nBytes);
	d_ap.FS  = device_alloc<DOUBLE> (nBytes);
	d_ap.tFS = device_alloc<DOUBLE> (nBytes);
	d_ap.Kx  = device_alloc<DOUBLE> (nBytes);
	d_ap.Ky  = device_alloc<DOUBLE> (nBytes);
	d_ap.dH  = device_alloc<DOUBLE> (nBytes);

	d_ap.ubt = device_alloc<DOUBLE> (sizeof(DOUBLE) * (M + 2));
	d_ap.ubt = device_alloc<DOUBLE> (sizeof(DOUBLE) * (M + 2));
	d_ap.vbt = device_alloc<DOUBLE> (sizeof(DOUBLE) * (N + 2));
	d_ap.vbt = device_alloc<DOUBLE> (sizeof(DOUBLE) * (N + 2));

	d_ap.khouot = device_alloc<int> (sizeof (int) * ap.u.size());
	d_ap.moci = device_alloc<int> (sizeof(int) * (N + 2));
	d_ap.mocj = device_alloc<int> (sizeof(int) * (M + 2));
	d_ap.daui = device_alloc<int> (sizeof(int) * segment_limit * (N + 2));
	d_ap.cuoi = device_alloc<int> (sizeof(int) * segment_limit * (N + 2));
	d_ap.dauj = device_alloc<int> (sizeof(int) * segment_limit * (M + 2));
	d_ap.cuoij = device_alloc<int> (sizeof(int) * segment_limit * (M + 2));
	// attention
	d_ap.hi = device_alloc<DOUBLE> (sizeof(DOUBLE) * (2 * (M + N + 6)));
	cudaError_t status = cudaMalloc((void**) device_arg_ptr, sizeof(Argument_Pointers));
	// assert sucess here
	cudaError_t copy_status = cudaMemcpy(*device_arg_ptr, &d_ap, sizeof(Argument_Pointers), cudaMemcpyHostToDevice);
	// assesrt success here
	// cuda copy arrays

	return d_ap;
}

Array_Pointers supporting_arrays_alloc(int M, int N, Array_Pointers** device_arr_ptr)
{
	Argument_Pointers d_ap;
	int nBytes = sizeof (DOUBLE) * (M * N + 2 * max(M, N));
	d_ap.a1 = device_alloc<DOUBLE> (nBytes);
	d_ap.b1 = device_alloc<DOUBLE> (nBytes);
	d_ap.c1 = device_alloc<DOUBLE> (nBytes);
	d_ap.d1 = device_alloc<DOUBLE> (nBytes);
	d_ap.a2 = device_alloc<DOUBLE> (nBytes);
	d_ap.c2 = device_alloc<DOUBLE> (nBytes);
	d_ap.ad = device_alloc<DOUBLE> (nBytes);
	d_ap.f1 = device_alloc<DOUBLE> (nBytes);
	d_ap.f2 = device_alloc<DOUBLE> (nBytes);
	d_ap.f3 = device_alloc<DOUBLE> (nBytes);
	d_ap.f5 = device_alloc<DOUBLE> (nBytes);

	int nBytes = sizeof (DOUBLE) * (2 * M * N + 4 * max(M, N));
	d_ap.AA = device_alloc<DOUBLE> (nBytes);
	d_ap.BB = device_alloc<DOUBLE> (nBytes);
	d_ap.CC = device_alloc<DOUBLE> (nBytes);
	d_ap.DD = device_alloc<DOUBLE> (nBytes);
	d_ap.Ap = device_alloc<DOUBLE> (nBytes);
	d_ap.Bp = device_alloc<DOUBLE> (nBytes);
	d_ap.ep = device_alloc<DOUBLE> (nBytes);	
	d_ap.x = device_alloc<DOUBLE> (nBytes);	

	d_ap.SN = device_alloc<int> (sizeof(int) * segment_limit * max(M, N));
	cudaError_t status = cudaMalloc((void**) device_arr_ptr, sizeof(Array_Pointers));
	// assert success here
	cudaError_t copy_status = cudaMemcpy(*device_arr_ptr, &d_ap, sizeof(Array_Pointers));
	// assert success here

	return d_ap;

}


void initialize(Host_arrays &host_arg_p, 
				Argument_Pointers* d_arg;
				Argument_Pointers** device_arg_ptr,
				Array_Pointers* d_arr;
				Array_Pointers** device_arr_ptr,
				bool initial_condition)
{
	if (initial_condition){
		load_initial_condition(&host_arr_p.u, &host_arg_p.v, &host_arg_p.z, &host_arg_p.FS, &host_arg_p.khouot);
	}
	*d_arg = attributes_arrays_memory_alloc(host_arg_p, device_arg_ptr);
	
	*d_arr = supporting_arrays_alloc(host_arg_p.M, host_arg_p.N, device_arr_ptr);

	// launches some kernels
	dim3 block_2d (min(M, 1024), 1, 1);
	dim3 grid_2d( (M + 2) / block_2d.x + 1, N + 3, 1);
	// launch onetime init kernel
	Onetime_init <<<block_2d, grid_2d>>>>(*device_arg_ptr);
	if (initial_condition)
		device_copy<DOUBLE>(&host_arg_p.khouot, d_arg->khouot);
	// launch some other kernels
	
	dim3 block_1d(32, 1, 1);
	dim3 grid_1dx(1, N, 1);
	dim3 grid_1dy(1, M, 1);
	dim3 grid_1d(1, 1, 1);

	Find_Calculation_limits_Horizontal <<<grid_1dx, block_1d>>> (device_arg_ptr);
	Find_Calculation_limits_Vertical <<<grid_1dy, block_1d>>> (device_arg_ptr);
	Htuongdoi <<<grid_2d , block_2d>>> (device_arg_ptr);
	preprocess <<<grid_1d, block_1d>>>(device_arg_ptr);
}


void update_boundary(DOUBLE t, int total_time, Argument_Pointers* device_arg_ptr)
{
	Update_Boundary_Value<<< >>> (t, total_time, device_arg_ptr),
}


// void simulation(int Tmax, bool ketdinh, bool channel=false,
// 	bool debug=false, bool plot=false, int interval=1, int sediment_start=10, int bed_change_start=20 ){
void simulation(Option ops, Argument_Pointers* device_arg_ptr, Array_Pointers* device_arr_ptr){
	// declare block size, grid size

	DOUBLE t = 0;
	int total_time = 0; // NEED TO CHANGE!!!
	while (t < Tmax){
		t += dT / 2;
		// boudary update
		update_boundary(t, total_time, device_arg_ptr);

		// set running limit for start, end



		// set blocksize, grid size


		// launch a bunch of kernels


		t+= dT / 2;

		// boudary update

		// set running limit for start, end


		// set blocksize, grid size


		// launch a bunch of kernels		

		// check if it's time to export resul

	}



}


// m is size of the matrix, which is SN
// __device__ void findBestGrid( int m, int tile_marshal, int *p_m_pad, int *p_b_dim, int *p_s, int *p_stride)
// {
//     int b_dim, m_pad, s, stride;
//     int B_DIM_MAX, S_MAX;
    
//     if ( sizeof(DOUBLE) == 4) {
//         B_DIM_MAX = 256;
//         S_MAX     = 512;    
//     }
//     else if (sizeof(DOUBLE) == 8){ /* double and complex */
//         B_DIM_MAX = 128;
//         S_MAX     = 256;     
//     }
    
//     /* b_dim must be multiple of 32 */
//     if ( m < B_DIM_MAX * tile_marshal ) {
//         b_dim = max( 32, (m/(32*tile_marshal))*32);
//         s = 1;
//         m_pad = ((m + b_dim * tile_marshal -1)/(b_dim * tile_marshal)) * (b_dim * tile_marshal);
//         stride = m_pad/(s*b_dim);    
//     }
//     else {
//         b_dim = B_DIM_MAX;
        
//         s = 1;
//         do {       
//             int s_tmp = s * 2;
//             int m_pad_tmp = ((m + s_tmp*b_dim*tile_marshal -1)/(s_tmp*b_dim*tile_marshal)) * (s_tmp*b_dim*tile_marshal);           
//             float diff = (float)(m_pad_tmp - m)/float(m);
//             /* We do not want to have more than 20% oversize */
//             if ( diff < .2 ) {
//                 s = s_tmp;      
//             }
//             else {
//                 break;
//             }
//         } while (s < S_MAX);
                       
//         m_pad = ((m + s*b_dim*tile_marshal -1)/(s*b_dim*tile_marshal)) * (s*b_dim*tile_marshal);        
//         stride = m_pad/(s*b_dim);
//     }
      
//     *p_stride = stride;
//     *p_m_pad  = m_pad;
//     *p_s      = s;
//     *p_b_dim  = b_dim;        
// }

// __global__ void copy_b_to_x(DOUBLE* b, DOUBLE*x, int m){
//     int i = blockIdx.x * blockDim.x + threadIdx.x ;
//     if (i >= m ) return;
//     x[i] = b[i];

// }

// __global__ void  tridiagSolver_v2(bool print, bool isU, int startidx, int endidx, int jumpstep, int tridiag_coeff_width, Argument_Pointers* arg, Array_Pointers * arr){

//     int i = blockIdx.x * blockDim.x + threadIdx.x + startidx;
//     if (i > endidx) return;

//     int number_of_segments;
//     int* dau, *cuoi;
//     if (isU){
//         number_of_segments = arg->mocj[i];
//         dau = arg->dauj;
//     }
//     else{ 
//         number_of_segments = arg->moci[i];
//         dau = arg->daui;

//         // if (blockIdx.y + 2 == 154 && threadIdx.x == 0)
//         //     printf("seg_no = %d\n", number_of_segments);
//     }
    
//     for (int j = 0; j < 1; j++){
//         int first = dau[i * segment_limit + j];
//         int pos = i * tridiag_coeff_width + first * jumpstep + jumpstep % 2; 
//         printf("j = %d %d\n", j, i);
//         // if ( j == 1 && i == 156)
//         //         printf("here %d\n", first);
//         DOUBLE* dl = &(arr->AA[pos]);
//         DOUBLE* d = &(arr->BB[pos]);
//         DOUBLE* du = &(arr->CC[pos]);
//         DOUBLE* b = &(arr->DD[pos]);
//         DOUBLE* x = &(arr->x[pos]);
//         int m = arr->SN[i * segment_limit + j] + 1;

        
//         // if ()
//         int s; //griddim.x
//         int stride;
//         int b_dim;
//         int m_pad;

//         int tile_marshal = 16;
//         int T_size = sizeof(DOUBLE);
        
//         findBestGrid( m, tile_marshal, &m_pad, &b_dim, &s, &stride);
       

        
//         int local_reduction_share_size = 2*b_dim*3*T_size;
//         int global_share_size = 2*s*3*T_size;
//         int local_solving_share_size = (2*b_dim*2+2*b_dim+2)*T_size;
//         int marshaling_share_size = tile_marshal*(tile_marshal+1)*T_size;
//         // printf("threadIdx=%d m=%d m_pad=%d s=%d b_dim=%d stride=%d\n", i, m, m_pad, s, b_dim, stride);
        
        
//         dim3 g_data(b_dim/tile_marshal,s);
//         dim3 b_data(tile_marshal,tile_marshal);
        

//         DOUBLE* dl_buffer;   //dl buffer
//         DOUBLE* d_buffer;    //b
//         DOUBLE* du_buffer; 
//         DOUBLE* b_buffer;
//         DOUBLE* w_buffer;
//         DOUBLE* v_buffer;
//         DOUBLE* c2_buffer;
        
//         DOUBLE* x_level_2;
//         DOUBLE* w_level_2;
//         DOUBLE* v_level_2;
        
        
//         //buffer allocation
//         cudaMalloc((void **)&c2_buffer, T_size*m_pad); 
//         cudaMalloc((void **)&dl_buffer, T_size*m_pad); 
//         cudaMalloc((void **)&d_buffer, T_size*m_pad); 
//         cudaMalloc((void **)&du_buffer, T_size*m_pad); 
//         cudaMalloc((void **)&b_buffer, T_size*m_pad); 
//         cudaMalloc((void **)&w_buffer, T_size*m_pad); 
//         cudaMalloc((void **)&v_buffer, T_size*m_pad); 
        
//         cudaMalloc((void **)&x_level_2, T_size*s*2); 
//         cudaMalloc((void **)&w_level_2, T_size*s*2); 
//         cudaMalloc((void **)&v_level_2, T_size*s*2); 

//         //kernels 
//         if (v_buffer == NULL){
//             printf("v_buffer is NULL, threadID = %d\n", i);
//             // return;
//         }
        
//         //data layout transformation
//         // foward_marshaling_bxb<<< g_data ,b_data, marshaling_share_size >>>(dl_buffer, dl, stride, b_dim, m, cuGet(0));
//         // foward_marshaling_bxb<<< g_data ,b_data, marshaling_share_size >>>(d_buffer,  d,  stride, b_dim, m, cuGet(1));
//         // foward_marshaling_bxb<<< g_data ,b_data, marshaling_share_size >>>(du_buffer, du, stride, b_dim, m, cuGet(0));
//         // foward_marshaling_bxb<<< g_data ,b_data, marshaling_share_size >>>(b_buffer,  b,  stride, b_dim, m, cuGet(0));
        
//         // cudaDeviceSynchronize();
//         // partitioned solver
//         // thomas_v1<<<s,b_dim>>>(b_buffer, w_buffer, v_buffer, c2_buffer, dl_buffer, d_buffer, du_buffer, stride);


//         // //SPIKE solver
//         // spike_local_reduction_x1<<<s,b_dim,local_reduction_share_size>>>(b_buffer,w_buffer,v_buffer,x_level_2, w_level_2, v_level_2,stride);
//         // spike_GPU_global_solving_x1<<<1,32,global_share_size>>>(x_level_2,w_level_2,v_level_2,s);
//         // spike_GPU_local_solving_x1<<<s,b_dim,local_solving_share_size>>>(b_buffer,w_buffer,v_buffer,x_level_2,stride);
//         // spike_GPU_back_sub_x1<<<s,b_dim>>>(b_buffer,w_buffer,v_buffer, x_level_2,stride);

//         // back_marshaling_bxb<<<g_data ,b_data, marshaling_share_size >>>(b,b_buffer,stride,b_dim, m);
//         // b_data = dim3(min(1024, (int) ceilf((float) m / 32) * 32), 1, 1);
//         // g_data = dim3( (int) ceilf((float) m / b_data.x), 1, 1);
//         // copy_b_to_x<<<g_data, b_data>>>(b, x, m);
        
//         //free
//         // cudaDeviceSynchronize();

//         __syncthreads(); // this is to make sure that all kernels have done their work on data before freeing them
//         cudaFree(dl_buffer);
//         cudaFree(d_buffer);
//         cudaFree(du_buffer);
//         cudaFree(b_buffer);
//         cudaFree(w_buffer);
//         cudaFree(v_buffer);
//         cudaFree(c2_buffer);
//         cudaFree(x_level_2);
//         cudaFree(w_level_2);
//         cudaFree(v_level_2);
//         // __syncthreads();
        
//     }
// }


// left over: total time and hmax u, d, l, r is not yet on GPU