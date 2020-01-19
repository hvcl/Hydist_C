#include "spike_kernel.hxx"

__device__
void findBestGrid( int m, int tile_marshal, int *p_m_pad, int *p_b_dim, int *p_s, int *p_stride);

__device__
void gtsv_spike_partial_diag_pivot_v1(const DOUBLE* dl, const DOUBLE* d, const DOUBLE* du, DOUBLE* b,const int m);