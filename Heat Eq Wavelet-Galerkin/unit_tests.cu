#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <assert.h>
#include <iostream>

void prvec(float *A, int n){
	for(int x=0; x<n; x++){
		std::cout<<A[x]<<"   ";
	}
}

void prmat(float *A, int n){
	for(int x=0; x<n; x++){
		for(int y=0; y<n; y++){
			std::cout<<A[x+y*n]<<"   ";
		}
		std::cout<<std::endl;
	}
}

//Solves Ax=b for x. (A is n-by-n and also symmetric).
//Based on CG.
void matSolve(const float *A, float *x, const float *b, const int n, const float absTol, cublasHandle_t handle){
	//Stolen from here: https://en.wikipedia.org/wiki/Conjugate_gradient_method
	std::cout<<"matSolve()"<<std::endl;
	float *r, *p, *Ap;
	float r_norm_old, r_norm;
	
	//For debugging. Remove. FIX!
	float *Q;
	Q = (float*) calloc (n*n,sizeof(*Q));
	
	
	cudaMalloc ((void**)&r, n*sizeof(*r));
	cudaMalloc ((void**)&p, n*sizeof(*p));
	cudaMalloc ((void**)&Ap, n*sizeof(*Ap));
	
	
	//Kolla om man får nån speed-up av att flytta dessa till GPUn. FIX!
	float alpha, beta, negative_alpha;
	
	//Kolla om man får nån speed-up av att flytta dessa till GPUn. FIX!
	float one = 1.0f;
	float negative_one = -1.0f;
	float zero = 0.0f;
	
	int counter = 0;
	
	//r = b-A*x
	cublasScopy(handle,n, b, 1, r,1);
	cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_one, A,n,x,1,&one,r, 1);
	
	//define r_norm_old
	cublasSdot(handle,n,r,1,r,1,&r_norm_old);
	//cublasSnrm2(handle,n,r,1, &r_norm_old);
	
	//p = r
	cublasScopy(handle, n, r, 1, p,1);
	
	//For debugging. Remove. FIX!
		std::cout<<std::endl;
		cublasGetVector(n, sizeof(float), p, 1, Q, 1);
		std::cout<<"p:"<<std::endl;
		prvec(Q,n);
		std::cout<<std::endl;
	
	while(counter<n){
		//std::cout<<"++counter"<<std::endl;
		++counter;
		
		//Ap = A*p
		//cublasSgbmv(handle,CUBLAS_OP_N,n,n,(n+1)/2,(n+1)/2,&one, A,n,p,1,&zero,Ap, 1);
		cublasSgemv(handle,CUBLAS_OP_N,n,n,&one, A,n,p,1,&zero,Ap, 1);
		
		//update alpha
		cublasSdot(handle,n,p,1,Ap,1,&alpha);
		alpha = r_norm_old/alpha;
		std::cout<<"alpha: "<<alpha<<std::endl;
		assert(alpha>=0.0f);
		
		//x = x + alpha*p
		cublasSaxpy(handle,n,&alpha, p,1, x, 1);
		
		//r = r - alpha*A*p (=r-alpha*Ap. FIX! cublasSaxpy() borde funka)
		negative_alpha = - alpha;
		//cublasSgbmv(handle,CUBLAS_OP_N,n,n,(n+1)/2,(n+1)/2,&negative_alpha, A,n,p,1,&one,r, 1);
		cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_alpha, A,n,p,1,&one,r, 1);
	
		//if r small: break
		cublasSdot(handle,n,r,1,r,1,&r_norm);
		//cublasSnrm2(handle,n,r,1, &r_norm);
		std::cout<<"r_norm: "<<r_norm<<std::endl;
		if(r_norm < absTol){
			break;
		}
		else{
			assert(r_norm != 0.f);
		}
		
		//update beta
		beta = r_norm/r_norm_old;
		
		//p = r + beta*p
		 cublasSscal(handle, n,&beta, p,1);
		 cublasSaxpy(handle, n, &one, r, 1, p, 1);
		 
		 r_norm_old = r_norm;
	}
	//assert(0);
		
	cudaFree(r);
	cudaFree(p);
	cudaFree(Ap);
}


int main(){

	//cudaError_t cudaStat;    
    cublasStatus_t stat;
    cublasHandle_t handle;
    
    stat = cublasCreate(&handle); 
	if (stat != CUBLAS_STATUS_SUCCESS) { 
		std::cout<<"CUBLAS initialization failed (blame main in unit_tests.cu)"<<std::endl;
		return EXIT_FAILURE;
	} 
	
	float *A,*devA, *x, *devX, *b, *devB;
	int n = 3;
	A = (float*) calloc(n*n,sizeof(*A));
	x = (float*) calloc(n,sizeof(*x));
	x[0] = 2;
	x[1] = 1;
	b = (float*) calloc(n,sizeof(*b));
	
	cudaMalloc((void**)&devA, n*n*sizeof(*A));
	cudaMalloc((void**)&devX, n*sizeof(*x));
	cudaMalloc((void**)&devB, n*sizeof(*b));
	
	/*
	A[0] = 2;
	A[1] = -1;
	A[2] = 0;
	A[3] = -1;
	A[4] = 2;
	A[5] = -1;
	A[7] = -1;
	A[8] = 2;
	
	b[1]=1;
	b[2]=0;
	b[3]=2;
	*/
	
	n=2;
	A[0] = 4;
	A[1] = 1;
	A[2] = 1;
	A[3] = 3;
	
	b[0] = 1;
	b[1] = 2;
	
	
	bool failure = false;
	//Sets failure=true (and keeps it that way) if any function call fails.
	failure = cublasSetMatrix(n,n, sizeof(*A), A, n, devA, n)!=CUBLAS_STATUS_SUCCESS || failure; 
	failure = cublasSetVector (n, sizeof(*x),x,1, devX, 1)!=CUBLAS_STATUS_SUCCESS || failure;
	failure = cublasSetVector (n, sizeof(*b),b,1, devB, 1)!=CUBLAS_STATUS_SUCCESS || failure;
	if (failure) { 
		std::cout<<"data download failed (1)"<<std::endl; 
		cudaFree (devA);
		cudaFree (devB);
		cudaFree (devX);
		cublasDestroy(handle); 
		return EXIT_FAILURE;
	}
	else{
		std::cout<<"data download successful (1)"<<std::endl;
	}
	
	matSolve(devA, devX, devB, n,0.01,handle);
	cublasGetVector(n, sizeof(*x), devX, 1, x,1);
	
	//This should show (7, 4/5, 1/15, 3).
	std::cout<<std::endl;
	for(int i=0; i<n; i++){
		std::cout<<x[i]<<std::endl;
	}
	
	return 0;
}

