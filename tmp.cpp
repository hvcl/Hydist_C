
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


// left over: total time and hmax u, d, l, r is not yet on GPU