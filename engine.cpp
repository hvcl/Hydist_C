#include "engine.h"
using namespace std;


void Load_coeffs(Constant_Coeffs& ce){

ce.dX = 10.0;
ce.dY = 10.0;
ce.dT = 2;
ce.dXbp = ce.dX  * ce.dX;
ce.dX2 = 2 * ce.dX;
ce.dYbp = ce.dY * ce.dY;
ce.dY2 = 2 * ce.dY;
ce.dTchia2dX = ce.dT / ce.dX2;
ce.dTchia2dY = ce.dT / ce.dY2;
ce.QuyDoiTime = 1.0 / 3600;
ce.QuyDoiPi = 1.0 / PI;
ce.HaiChiadT = 2.0 / ce.dT;
ce.kenhhepng = 0;
ce.kenhhepd = 0;

// gioi han tinh
ce.NANGDAY = 0.0; // thong so nang day
ce.H_TINH = 0.02; // do sau gioi han (m)

// thong so ve gio
ce.Wind = 0.0;  // van toc gio (m/s)
ce.huonggio = 0.0;  // huong gio (degree)

// Zban dau
ce.Zbandau = 0.0;

// He so lam tron va he so mu manning
ce.heso = 0.94; //1.0,
ce.mu_mn = 0.2;
// ND number (kg/m3)
ce.NDnen = 0.00;
ce.NDbphai = 0.5;
ce.NDbtrai = 0.5;
ce.NDbtren = 0.5;
ce.NDbduoi = 0.5;

// tod(toi han boi), toe(toi han xoi)
// hstoe (he so tinh ung suat tiep toi han xoi theo do sau)
// ghtoe (gioi han do sau tinh toe(m))
// Mbochat (kha nang boc hat M(kg/m2/s))
ce.Tod = 0.06;
ce.Toe = 0.15;
ce.hstoe = 0;
ce.ghtoe = 3;
ce.Mbochat = 0.00001;

// khoi luong rieng cua nuoc (ro) va khoi luong rieng cua hat (ros) (kg/m3)
ce.ro = 1000;
ce.ros = 2000;

// duong kinh trung binh cua hat 50% (m) (dm)
ce.dm = 0.00001;
// duong kinh hat trung binh 90% (m) 
ce.d90 = 0.002;

// he so nhot dong hoc cua nuoc sach
ce.muy = 1.01e-06;

// Do rong cua hat (Dorong) va Ty trong (KLR cua hat va nuoc) (Sx)
ce.Dorong = 0.5;
ce.Sx = 2;

//tong so do sau de tinh he so nham
ce.sohn = 8;
//tong so do sau de tinh he so nhot
ce.soha = 3;
// tong so do sau de tinh Fw
ce.sohn1 = 3;
// //luc coriolis
ce.CORIOLIS_FORCE = 0.0;
ce.g = 9.81; //gia toc trong truong = 9.81 m2/s

// //lan truyen
ce.Ks = 2.5 * ce.dm; // = 2.5 * dm
// hoac bang Ks = 3 * d 90% 
Windx = 0.0013 * (0.00075 + 0.000067 * abs(Wind)) * abs(Wind) * Wind * cos(huonggio * (PI / 180));
Windy = 0.0013 * (0.00075 + 0.000067 * abs(Wind)) * abs(Wind) * Wind * sin(huonggio * (PI / 180));


}


template <typename T>
pair<int, int> load_file(const char* filename, vector<T> &arr){

	ifstream ifs;
	ifs.open(filename);
	if (!ifs) {
		cout << "cannot open" <<  filename << endl;
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
	return pair<int, int>(N, M);
}

template pair<int, int> load_file<double>(const char* filename, std::vector<double> &arr);
template pair<int, int> load_file<float>(const char* filename, std::vector<float> &arr);


pair<int, int> load_inputs(string dir,
				vector<DOUBLE> &h,
				vector<DOUBLE> &hsnham,
				vector<DOUBLE> &VISCOINDX, 
				vector <DOUBLE> &Fw,
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
	size = load_file<DOUBLE> ((dir + "bandodosau.txt").c_str(), h);
	cout << size.first << " " << size.second << endl;
	load_file<DOUBLE> ((dir + "hsnham.txt").c_str(), hsnham);

	cout << size.first << " " << size.second << endl;
	load_file<DOUBLE> ((dir + "hsnhotroiA.txt").c_str(), VISCOINDX);

	load_file<DOUBLE> ((dir + "Fw_map.txt").c_str(), Fw);

	// load_file<DOUBLE> (dir + "bientren.txt", bc_up);
	load_file<DOUBLE> ((dir + "bienduoi.txt").c_str(), bc_down);
	// load_file<DOUBLE> (dir + "bientrai.txt", bc_left);
	load_file<DOUBLE> ((dir + "bienphai.txt").c_str(), bc_right);
	// load_file<DOUBLE> (dir + "FSbientren.txt", CC_u);
	// load_file<DOUBLE> (dir + "FSbienduoi.txt", CC_d);
	// load_file<DOUBLE> (dir + "FSbientrai.txt", CC_l);
	// load_file<DOUBLE> (dir + "FSbienphai.txt", CC_r);
	load_file<int> ((dir + "boundary_type.txt").c_str(), bienQ);
	return size;
}

void load_initial_condition(string dir, 
							vector<DOUBLE> &u, 
							vector<DOUBLE> &v,
							vector<DOUBLE> &z,
							vector<DOUBLE> &FS,
							vector<int> &khouot)
{

	load_file<DOUBLE> ((dir + "u.txt").c_str(), u);
	load_file<DOUBLE> ((dir + "v.txt").c_str(), v);
	load_file<DOUBLE> ((dir + "z.txt").c_str(), z);
	load_file<DOUBLE> ((dir + "FS.txt").c_str(), FS);
	load_file<int> ((dir + "khouot.txt").c_str(), khouot);
}


template <typename T>
T* device_alloc_and_copy(vector<T> &h_array){
	int nBytes = h_array.size() * sizeof(T);
	T* d_array;
	cudaError_t malloc_status = cudaMalloc((void**) &d_array, nBytes);
	// assert success here
	cudaError_t copy_status = cudaMemcpy(d_array, h_array.data(), nBytes, cudaMemcpyHostToDevice);
	// asser success here
	assert(copy_status == cudaSuccess);
	// returna address of array on device here
	return d_array;

}

template <typename T>
T* device_alloc(int nBytes){
	T* d_array;
	cudaError_t status = cudaMalloc((void**) &d_array, nBytes);
	// assert success
	assert(status == cudaSuccess);
	return d_array;
}

template <typename T>
void device_copy(vector<T> &source, T* des)
{
	cudaError_t status =  cudaMemcpy(des, source.data(), sizeof(T) * source.size(), cudaMemcpyHostToDevice);
	// assert success here
	assert(status == cudaSuccess);

}


/*
	This function 
*/
Argument_Pointers attribute_arrays_memory_alloc(int device, Host_arrays &ap, Argument_Pointers** device_arg_ptr)
{
	// cuda mem alloc corresponding arrays
	int M = ap.M;
	int N = ap.N;
	Argument_Pointers d_ap ;
	d_ap.M = M;
	d_ap.N = N;
	d_ap.hmax_u = 0;
	d_ap.hmax_d = 0;
	d_ap.hmax_l = 0;
	d_ap.hmax_r = 0;

	// maps
	d_ap.h = device_alloc_and_copy<DOUBLE> (ap.h);
	d_ap.hsnham = device_alloc_and_copy<DOUBLE> (ap.hsnham);
	d_ap.VISCOIDX = device_alloc_and_copy<DOUBLE> (ap.VISCOIDX);
	d_ap.Fw = device_alloc_and_copy<DOUBLE> (ap.Fw);

	// boundary condition for FS
	d_ap.CC_u = device_alloc_and_copy<DOUBLE> (ap.CC_u);
	d_ap.CC_d = device_alloc_and_copy<DOUBLE> (ap.CC_d);
	d_ap.CC_l = device_alloc_and_copy<DOUBLE> (ap.CC_l);
	d_ap.CC_r = device_alloc_and_copy<DOUBLE> (ap.CC_r);

	// boundary conditions
	d_ap.bienQ = device_alloc_and_copy<int> (ap.bienQ);
	d_ap.boundary_type = device_alloc_and_copy<int> (ap.boundary_type);

	d_ap.bc_up = device_alloc_and_copy<DOUBLE> (ap.bc_up);
	d_ap.bc_down = device_alloc_and_copy<DOUBLE> (ap.bc_down);
	d_ap.bc_left = device_alloc_and_copy<DOUBLE> (ap.bc_left);
	d_ap.bc_right = device_alloc_and_copy<DOUBLE> (ap.bc_right);

	// intitial conditions
	// need to double check this, to do this only if initial condition is given
	d_ap.u = device_alloc_and_copy<DOUBLE> (ap.u);
	d_ap.v = device_alloc_and_copy<DOUBLE> (ap.v);
	d_ap.z = device_alloc_and_copy<DOUBLE> (ap.z);
	d_ap.FS = device_alloc_and_copy<DOUBLE> (ap.FS);
	d_ap.khouot = device_alloc_and_copy<int> (ap.khouot);

	d_ap.khouot = device_alloc_and_copy<int> (ap.khouot);

	// cuda mem alloc arrays on device only

	int nBytes = ap.u.size() * sizeof (DOUBLE);
	d_ap.t_u = device_alloc<DOUBLE> (nBytes);
	d_ap.t_v = device_alloc<DOUBLE> (nBytes);
	d_ap.t_z = device_alloc<DOUBLE> (nBytes);
	d_ap.Kx1 = device_alloc<DOUBLE> (nBytes);
	d_ap.Ky1 = device_alloc<DOUBLE> (nBytes);
	d_ap.Htdu = device_alloc<DOUBLE> (nBytes);
	d_ap.Htdv = device_alloc<DOUBLE> (nBytes);
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

	
	d_ap.moci = device_alloc<int> (sizeof(int) * (N + 2));
	d_ap.mocj = device_alloc<int> (sizeof(int) * (M + 2));
	d_ap.daui = device_alloc<int> (sizeof(int) * segment_limit * (N + 2));
	d_ap.cuoii = device_alloc<int> (sizeof(int) * segment_limit * (N + 2));
	d_ap.dauj = device_alloc<int> (sizeof(int) * segment_limit * (M + 2));
	d_ap.cuoij = device_alloc<int> (sizeof(int) * segment_limit * (M + 2));

	// attention
	d_ap.hi = device_alloc<DOUBLE> (sizeof(DOUBLE) * (2 * (M + N + 6)));
	cudaError_t status = cudaMalloc((void**) device_arg_ptr, sizeof(Argument_Pointers));
	// assert sucess here
	assert(status == cudaSuccess);

	cudaError_t copy_status = cudaMemcpy(*device_arg_ptr, &d_ap, sizeof(Argument_Pointers), cudaMemcpyHostToDevice);
	// assesrt success here
	assert(status == cudaSuccess);

	return d_ap;
}


Array_Pointers supporting_arrays_alloc(int M, int N, Array_Pointers** device_arr_ptr)
{
	Array_Pointers d_ap;
	int nBytes = sizeof (DOUBLE) * (M * N + 2 * max(M, N));
	d_ap.a1 = device_alloc<DOUBLE> (nBytes);
	d_ap.b1 = device_alloc<DOUBLE> (nBytes);
	d_ap.c1 = device_alloc<DOUBLE> (nBytes);
	d_ap.d1 = device_alloc<DOUBLE> (nBytes);
	d_ap.a2 = device_alloc<DOUBLE> (nBytes);
	d_ap.c2 = device_alloc<DOUBLE> (nBytes);
	d_ap.d2 = device_alloc<DOUBLE> (nBytes);
	d_ap.f1 = device_alloc<DOUBLE> (nBytes);
	d_ap.f2 = device_alloc<DOUBLE> (nBytes);
	d_ap.f3 = device_alloc<DOUBLE> (nBytes);
	d_ap.f5 = device_alloc<DOUBLE> (nBytes);

	nBytes = sizeof (DOUBLE) * (2 * M * N + 4 * max(M, N));
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
	// cout << "device copy: " << cudaSuccess << endl;
	assert(status == cudaSuccess);


	cudaError_t copy_status = cudaMemcpy(*device_arr_ptr, &d_ap, sizeof(Array_Pointers), cudaMemcpyHostToDevice);
	// assert success here
	assert(status == cudaSuccess);


	return d_ap;

}





