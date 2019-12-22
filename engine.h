#ifndef ENGINE_H__
#define ENGINE_H__
#include <vector>
#include <algorithm>
#include <cmath>
#include "cuda.h"
#include "cuda_runtime.h"
using namespace std;

#define DOUBLE double
#define	segment_limit 20
#define PI 3.14159265358979323846

struct Options{
	DOUBLE Tmax, t_start, sediment_start, bed_change_start;
	int interval;
	bool kenhhepd, kenhhepng;
	int M, N;
	bool ketdinh, channel, debug, plot;
	int total_time;
	Options(DOUBLE Tmx, DOUBLE ts, int itv, DOUBLE sds, DOUBLE bcs, bool cohesive, bool kenhng,
			bool kenhd, bool db, bool plt){
		Tmax = Tmx;
		t_start = ts;
		interval = itv;
		sediment_start = sds;
		bed_change_start = bcs;
		kenhhepng = kenhng;
		kenhhepd = kenhd;
		ketdinh = cohesive;
		channel = kenhd ^ kenhng;
		debug = db;
		plot = plt;
	}
};

struct Constant_Coeffs{
	DOUBLE dX, dY, dT, dXbp, dYbp, dX2, dY2, dTchia2dX, dTchia2dY, HaiChiadT;
	DOUBLE QuyDoiTime, QuyDoiPi;
	DOUBLE NANGDAY, H_TINH, heso, mu_mn;
	DOUBLE kenhhepng, kenhhepd;
	DOUBLE Wind, huonggio, Zbandau;
	DOUBLE NDnen, NDbphai, NDbtrai, NDbtren, NDbduoi;
	DOUBLE Tod, Toe, hstoe, ghtoe , Mbochat;
	DOUBLE ro, ros, dm, muy, d90, Dorong, Sx, sohn, soha, sohn1;
	DOUBLE CORIOLIS_FORCE, g, Ks;
	DOUBLE Windx, Windy; 


};


struct Argument_Pointers {
	int M, N;
	DOUBLE hmax_u, hmax_d, hmax_l, hmax_r;
	int *bienQ;
	int* daui, *dauj, *cuoii, *cuoij, *moci, *mocj, *khouot, *boundary_type;
	DOUBLE* h, *v, *u, *z, *t_u, *t_v, *t_z, *Htdu, *Htdv, *H_moi, *htaiz, *htaiz_bd;
	DOUBLE* ubt, *ubp, *vbt, *vbd; //*vt, *ut;
	DOUBLE* hsnham, *VISCOIDX, *Kx1, *Ky1, *Tsyw, *Tsxw;
	DOUBLE* bc_up, *bc_down, *bc_left, *bc_right;
	DOUBLE* hi;
	DOUBLE* FS, *tFS, *CC_u, *CC_d, *CC_l, *CC_r;
	DOUBLE *VTH, *Kx, *Ky, *Fw;
	DOUBLE* Qbx, *Qby;
	DOUBLE* dH;
};

struct Array_Pointers {
	DOUBLE* a1, *b1, *c1, *d1, *a2, *c2, *d2;
	DOUBLE* f1, *f2, *f3, *f5;
	DOUBLE* AA, *BB, *CC, *DD;
	DOUBLE* x;
	DOUBLE* Ap, *Bp, *ep;
	int *SN;
};

struct  Host_arrays
{	
	int M, N, total_time;
	vector<DOUBLE> h, hsnham, VISCOIDX, Fw, FS, dH;
	vector<DOUBLE> bc_up, bc_down, bc_left, bc_right;
	vector<DOUBLE> CC_u, CC_d, CC_l, CC_r;
	vector<DOUBLE> u, v, z;
	vector<int> khouot, bienQ, boundary_type;
};


void update_boundary_at_t(int M, int N, float t, bool channel, int total_time, Argument_Pointers* d_arg_ptr, Constant_Coeffs* coeffs);

void synch_and_check();

template <typename T>
void save_file(T* array, int width, int height, const char* filename);

void Hydraulic_Calculation(DOUBLE dT, DOUBLE NANGDAY, Argument_Pointers* d_arg_ptr, Array_Pointers* d_arr_ptr, Constant_Coeffs* coeffs,
						 Array_Pointers h_arg_ptr,  Options ops);


#endif