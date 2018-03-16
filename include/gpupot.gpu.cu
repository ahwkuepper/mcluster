//#include "gpupot.h"
// #include <iostream>
#include <cstdio>
// #include <cutil.h>
#ifdef WITH_CUDA5
#  include <helper_cuda.h>
#  define CUDA_SAFE_CALL checkCudaErrors
#else
#  include <cutil.h>
#endif
#include "cuda_pointer.h"
#define NTHREAD 128

#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

static float2 float2_split(double x){
	const int shift = 20;
	float2 ret;
	x *= (1<<shift);
	double xi = (int)x;
	double xf = x - xi;
	ret.x = xi * (1./(1<<shift));
	ret.y = xf * (1./(1<<shift));
	return ret;
}
__device__ float2 float2_accum(float2 acc, float x){
	float tmp = acc.x + x;
	acc.y -= (tmp - acc.x) - x;
	acc.x = tmp;
	return acc;
}

__device__ float2 float2_regularize(float2 acc){
	float tmp = acc.x + acc.y;
	acc.y = acc.y -(tmp - acc.x);
	acc.x = tmp;
	return acc;
}

__device__ float2 float2_add(float2 a, float2 b){
  float tmp = a.x + b.x;
  a.y -= (tmp - a.x) - b.x - b.y;
  a.x = tmp;
  // a.x = a.x + b.x;
  // a.y = a.y + b.y;
  return a;
}

struct Particle{
	float2 pos[3];
	float mass;
	float pad;

	Particle(double x[3], double m){
		pos[0] = float2_split(x[0]);
		pos[1] = float2_split(x[1]);
		pos[2] = float2_split(x[2]);
		mass = (float)m;
	}
	Particle(int){
		pos[0].x = pos[0].y = pos[1].x = pos[1].y = pos[2].x = pos[2].y = mass = pad = 0.f;
	}
	__device__ Particle() {}
};

__global__ void pot_kernel(int n, int nblock, Particle *ptcl, float2 *phi){
	__shared__ Particle jpbuf[NTHREAD];
	int i = NTHREAD * blockIdx.x + threadIdx.x;
	Particle ip = ptcl[i];
	float2 phii = make_float2(0.f, 0.f);
	for(int j=0; j<n; j+= NTHREAD){
		__syncthreads();
		jpbuf[threadIdx.x] = ptcl[j + threadIdx.x];
		__syncthreads();
#pragma unroll 4
		for(int jj=0; jj<NTHREAD; jj++){
			// if(j+jj == i) continue;
			Particle &jp = jpbuf[jj];
			float dx = (jp.pos[0].x - ip.pos[0].x) + (jp.pos[0].y - ip.pos[0].y);
			float dy = (jp.pos[1].x - ip.pos[1].x) + (jp.pos[1].y - ip.pos[1].y);
			float dz = (jp.pos[2].x - ip.pos[2].x) + (jp.pos[2].y - ip.pos[2].y);
			float r2 = dx*dx + dy*dy + dz*dz;
			// if(r2==0.f) continue;
			float pij = jp.mass * rsqrtf(r2);
			// phii = float2_accum(phii, pij);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
		phii = float2_regularize(phii);
	}
	phi[i] = phii;

    // for(int j = nblock/2 + nblock%2; j>1; j = j/2 + j%2) {
    //   int offset = j%2;
    //   if(blockIdx.x<j-offset) phi[i] = float2_add(phi[i],phi[i + j*NTHREAD]);
    // }
    
    // for(int j = NTHREAD/2; j>1; j/= 2) {
    //   if(threadIdx.x<j) phi[i] = float2_add(phi[i],phi[i + j]);
    // }
}

extern "C" void gpupot(int n, double **star, double *pot) {

    int numGPU=0;
    cudaGetDeviceCount(&numGPU);
    assert(numGPU>0);
    cudaSetDevice(0);

	double t0 = get_wtime();

	cudaPointer <float2> phi;
	cudaPointer <Particle> ptcl;

	int ng = NTHREAD * (n/NTHREAD + (n%NTHREAD ? 1 : 0));

	phi.allocate(ng);
	ptcl.allocate(ng);

	// std::cout << n << " " << ng << std::endl;
	for(int i=0; i<n; i++){
		// ptcl_h[i] = Particle(x[i], m[i]);
		ptcl[i] = Particle(&star[i][1], star[i][0]);
	}
	for(int i=n; i<ng; i++){
		// ptcl_h[i] = Particle(0);
		ptcl[i] = Particle(0);
	}

	// cudaMemcpy(ptcl_d, ptcl_h, ng * sizeof(Particle), cudaMemcpyHostToDevice);
	ptcl.htod(ng);
	
	dim3 grid(ng/NTHREAD, 1, 1);
	dim3 threads(NTHREAD, 1, 1);
	int sharedMemSize = NTHREAD * sizeof(Particle);
	// pot_kernel <<<grid, threads, sharedMemSize >>> (n, ptcl_d, phi_d);
    pot_kernel <<<grid, threads, sharedMemSize >>> (n, ng/NTHREAD, ptcl, phi);

    // phi.dtoh(1);
    // double pot = (double)phi[0].x + (double)phi[0].y;
	// cudaMemcpy(phi_h, phi_d, n * sizeof(float2), cudaMemcpyDeviceToHost);
    phi.dtoh(n);
    *pot = 0;
    for(int i=0; i<n; i++){
		// pot[i] = (double)phi_h[i].x + (double)phi_h[i].y;
      *pot -= star[i][0]*((double)phi[i].x + (double)phi[i].y);
    }
    *pot /=2.;
    
	phi.free();
	ptcl.free();

	double t1 = get_wtime();
#ifdef PROFILE
	fprintf(stderr, "gpupot: %f sec\n", t1 - t0);
#endif
}

