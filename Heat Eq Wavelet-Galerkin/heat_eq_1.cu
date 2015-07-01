/*
	Solves the heat eq in 2D using CUDA.
*/

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <assert.h>
#include <iostream>

#define N 5
#define PI 3.1415926f


__device__ __host__
float phi(float x, float y, float xc, float yc, int lvl){
	/*
		This is a wavelet function (or basis function in FEM-speak).
		x & y are the variables.
		xc & yc defines the point around which it is centered.
		lvl defines the resolution. (Zero is the most course resolution)
		
		For now this is the mexican hat function.
	*/
	float sigma = 0.05f; //FIX!
	float tmp = -(x*x+y*y)/(2.0f*sigma*sigma);
	return -(1.0f+tmp)*exp(tmp)/(PI*sigma*sigma*sigma*sigma);
}

__device__ __host__
float phi_xx(float x, float y, float xc, float yc, int lvl){
	/*
		This is the second derivative wrt x of the phi-function.
	*/
	return 0.0f; //FIX!s
}

__device__ __host__
float phi_yy(float x, float y, float xc, float yc, int lvl){
	/*
		This is the second derivative wrt y of the phi-function.
	*/
	return 0.0f; //FIX!
}

void lowerTriGpu(const float *A, float *L, const int n){
	for(int x=0; x<n; x++){
		for(int y=0; y<n; y++){
			if(y<=x){
				L[x*n+y] = A[x*n+y];
			}
			else{
				L[x*n+y] = 0.0f;
			}
		}
	}
}

void upperTriGpu(const float *A, float *U, const int n){
	for(int x=0; x<n; x++){
		for(int y=0; y<n; y++){
			if(y>x){
				U[x*n+y] = A[x*n+y];
			}
			else{
				U[x*n+y] = 0;
			}
		}
	}
}

__host__
void matSolvGpu(float *A, float *x, float *b, int n){
	/*
		Solves the eq Ax=b using Gauss-Seidel on the
		GPU using CUBLAS.
	*/
	
	//A = L + U
	float *L;
	float *U;
	
	(float *) malloc(N * N * sizeof (*L));
	(float *) malloc(N * N * sizeof (*U));
	
	lowerTriGpu(A, L, n);
	upperTriGpu(A, U, n);
	
}

int main(){
	
	float *A=0;
	float *devPtrA;
	
	A = (float *) malloc(N * N * sizeof (*A));
	
	
	cudaError_t cudaStat;    
    cublasStatus_t stat;
    cublasHandle_t handle;

	cudaStat = cudaMalloc ((void**)&devPtrA, N*N*sizeof(*A));
	if (cudaStat != cudaSuccess) { 
		std::cout<<"device memory allocation failed"<<std::endl;
		return EXIT_FAILURE;
	}
	stat = cublasCreate(&handle); 
	if (stat != CUBLAS_STATUS_SUCCESS) { 
		std::cout<<"CUBLAS initialization failed"<<std::endl;
		return EXIT_FAILURE;
	} 
	stat = cublasSetMatrix(N, N, sizeof(*A), A, N, devPtrA, N); 
	if (stat != CUBLAS_STATUS_SUCCESS) { 
		std::cout<<"data download failed"<<std::endl;
		cudaFree (devPtrA);
		cublasDestroy(handle); 
		return EXIT_FAILURE;
	}


	cudaFree (devPtrA);
    cublasDestroy(handle);
    
    free(A);
	
	return 0;
}

