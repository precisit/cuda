#define datatype float

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <cmath>

__device__ int myMax(const int i, const int j){
	if(i>j){
		return i;
	}
	else{
		return j;
	}
}

__device__ datatype getLaplacian(const int i, const int j, const datatype h, const int n){
	int dist = abs(i-j);
	if( dist == 0){
		return -4.0f/(h*h);
	}
	else if( dist == 1 ){
		if (myMax(i,j) % n == 0){
			return 0.0f;
		}
		else{
			return 1.0f/(h*h);
		}
	}
	else if( dist == n ){
		return 1.0f/(h*h);
	}
	else{
		return 0.0f;
	}
}

__global__ void jacobiSmootherLaplacianDev(const datatype* u, const datatype* b, const int n, datatype* u_new){
			
			//	Smoothes the error to D u = b.
			

	const datatype h = 1.0f/( (datatype) n-1.0f);
	datatype tmpSum;
	int i = threadIdx.x;
	if(i< n*n){
		tmpSum = 0.0f;
		for(int j=0; j < n*n; j++){
			if(j != i){
				tmpSum += getLaplacian(i,j,h,n)*u[j];
			}
		}
		u_new[i] = (b[i]-tmpSum) / getLaplacian(i,i,h,n);
	}
}


__global__ void calculateErrorLaplacian(const datatype* u, const datatype* b, datatype* d, const int n){

	int i = threadIdx.x;
	if(i<n*n){
		const datatype h = 1.0f/((datatype) n-1.0f);
		datatype sum = 0.0f;
		for(int j=0; j < n*n; j++){
			sum += getLaplacian(j,i,h,n)*u[j];
		}
		d[i] = b[i]-sum;
	}
};

__global__ void restrictDtoB(){



};

__global__ void interpolate(){

};

int main(){

	const int n = 5;

	datatype *u, *b, *u_new;

	//Allocate memory
	u = (datatype*) malloc(n*n*sizeof(datatype));
	u_new = (datatype*) malloc(n*n*sizeof(datatype));
	b = (datatype*) malloc(n*n*sizeof(datatype));

	datatype *u_dev, *b_dev, *u_dev_new;
	cudaMalloc( (void**) &u_dev, n*n*sizeof(datatype));
	cudaMalloc( (void**) &u_dev_new, n*n*sizeof(datatype));
	cudaMalloc( (void**) &b_dev, n*n*sizeof(datatype));

	//initiate data
	for(int i=0; i<n*n; i++){
		u[i] = 0.0f;
		b[i] = 1.0f;
	}

	//Copy to GPU
	cudaMemcpy(u_dev, u, n*n*sizeof(datatype), cudaMemcpyHostToDevice);
	cudaMemcpy(u_dev_new, u_new, n*n*sizeof(datatype), cudaMemcpyHostToDevice);
	cudaMemcpy(b_dev, b, n*n*sizeof(datatype), cudaMemcpyHostToDevice);

	dim3 block_size(n*n);

	for(int i=0; i<150; i++){
		jacobiSmootherLaplacianDev<<<1,block_size>>>(u_dev, b_dev, n, u_dev_new);
		jacobiSmootherLaplacianDev<<<1,block_size>>>(u_dev_new, b_dev, n, u_dev);
	}

	calculateErrorLaplacian<<<1,block_size>>>(u_dev, b_dev, u_dev_new,n);

	//Copy data from GPU
	cudaMemcpy(u_new, u_dev_new, n*n*sizeof(datatype), cudaMemcpyDeviceToHost);
	cudaMemcpy(u, u_dev, n*n*sizeof(datatype), cudaMemcpyDeviceToHost);
	cudaMemcpy(b, b_dev, n*n*sizeof(datatype), cudaMemcpyDeviceToHost);

	//Print the result
	for(int i=0; i<n*n; i++){
		std::cout<<u_new[i]<<std::endl;
	}
	for(int i=0; i<n*n; i++){
		std::cout<<u[i]<<std::endl;
	}

	//Free memory
	free(b);
	free(u);
	free(u_new);

	cudaFree(u_dev);
	cudaFree(u_dev_new);
	cudaFree(b_dev);

	return 0;
}
