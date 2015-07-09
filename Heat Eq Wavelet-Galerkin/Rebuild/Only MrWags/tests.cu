/*

	Tests the status = cublasSgemv-function. It behaves weirdly in MrWags.
	
*/


#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <assert.h>
#include <iostream>
#include <stdio.h>








static const char *_cudaGetErrorEnum(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED\n";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}






int main(){

   
    cublasStatus_t status;
    cublasHandle_t handle;
    cudaError_t cudaStat;

	cublasCreate(&handle);
	
	int n = 4;
	
	float one = 1.0f;
	float negative_one = -1.0f;
	float alpha = -0.000001f;
	
	float *devA, *devB, *devX, *A, *b, *x;
	
	cudaMalloc ((void**)&devB, n*sizeof(*devB));
	cudaMalloc ((void**)&devA, n*n*sizeof(*devA));
	cudaMalloc ((void**)&devX, n*sizeof(*devX));
	
	
	A = (float*) calloc (n*n,sizeof(*A));
	x = (float*) calloc (n,sizeof(*x));
	b = (float*) calloc (n,sizeof(*b));
	
	A[0] = 1;
	A[5] = 1;
	A[10] = 1;
	A[15]= 1;
	
	b[0] = 3;
	
	x[0] = 0.5;
	x[2] = 1;
	
	status = cublasSetMatrix(n,n, sizeof(float), A, n, devA,n);
	if(status != CUBLAS_STATUS_SUCCESS){
		printf(_cudaGetErrorEnum(status));
		std::cout<<"1"<<std::endl;
		return EXIT_FAILURE;
	}
	status = cublasSetVector(n, sizeof(float), b, 1, devB,1);
	if(status != CUBLAS_STATUS_SUCCESS){
		printf(_cudaGetErrorEnum(status));
		std::cout<<"2"<<std::endl;
		return EXIT_FAILURE;
	}
	status = cublasSetVector(n, sizeof(float), x, 1, devX,1);
	if(status != CUBLAS_STATUS_SUCCESS){
		printf(_cudaGetErrorEnum(status));
		std::cout<<"3"<<std::endl;
		return EXIT_FAILURE;
	}
	
	status = cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_one, devA,n,devX,1,&one,devB, 1);

	if(status != CUBLAS_STATUS_SUCCESS){
		printf(_cudaGetErrorEnum(status));
		return EXIT_FAILURE;
	}

	assert(status == CUBLAS_STATUS_SUCCESS);
	
	

	cublasGetMatrix(n,n, sizeof(float), devA, n, A,n);
	cublasGetVector(n, sizeof(float), devB, 1, b,1);
	cublasGetVector(n, sizeof(float), devX, 1, x,1);
	
	for(int i= 0; i<n; i++){
		std::cout<<"b: "<<b[i]<<std::endl;
	}

	return 0;
};
