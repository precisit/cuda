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
	const int dist = abs(i-j);
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
	const int i = threadIdx.x;
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

	const int i = threadIdx.x;
	if(i<n*n){
		const datatype h = 1.0f/((datatype) n-1.0f);
		datatype sum = 0.0f;
		for(int j=0; j < n*n; j++){
			sum += getLaplacian(j,i,h,n)*u[j];
		}
		d[i] = b[i]-sum;
	}
};

__global__ void restrictDtoB(const datatype* d, const int d_len, datatype* b, const int b_len){

	const int x = threadIdx.x + blockIdx.x*blockDim.x;
	const int y = threadIdx.y + blockIdx.y*blockDim.y;

	if(x < d_len && y < d_len){
		if( x%2==0 && y%2==0){
			b[ (x/2) * b_len + (y/2)] = d[x*d_len + y];
		}
	}


};

__global__ void interpolate(const datatype* from, const int from_len, datatype* to, const int to_len){

	//to_len > from_len, obvs!

	const int x = threadIdx.x;
	const int y = threadIdx.y;

	if(x< to_len && y<to_len){
		if(x<from_len && y<from_len){
			//Copy the nodal values
			to[(2*x)*to_len+2*y] = from[x*from_len+y];
		}
		syncthreads();
		
		//Interpolate in the x-dir (might be y-dir. It's possible I fucked up here. But it doesn't really matter.)
		if(y >= 1 && x<from_len && y<from_len ){
			to[2*x*to_len+2*y-1] = (from[x*from_len+y]+from[x*from_len+y-1])/2.0f;
		}

		syncthreads();

		//Interpolate in the y-dir
		if(x >= 1 && x%2==1){
			to[x*to_len+y] = (to[(x-1)*to_len+y]+to[(x+1)*to_len+y])/2.0f;
		}
	}

};

datatype* firstPonterFromGridLevel(const int* nodesPerLevel, const int level, const int numberOfLevels, datatype* u){
	/*
		Level = 0 is the coarsest grid.
	*/
	assert(level<numberOfLevels);
	assert(level>=0);

	//backward level is coarsest when = numberOfLevels -1 and finest when = 0
	int backwardLevel = numberOfLevels - level -1;
	datatype* ptr = u;
	for(int i = 0; i<backwardLevel; i++){
		if( i = backwardLevel ){
			return ptr;
		}
		ptr += nodesPerLevel[i];
	}
	assert(0);
	return ptr;
	
}
__global__ void addVectors(datatype* y, const datatype* x, const int n){
	// y = y + x
	// Make sure that they are all of length n
	const int i= threadIdx.x + blockIdx.x*blockDim.x;

	if(i<n){
		y[i] += x[i];
	}
}


void multigrid(const int k, datatype* u, datatype* b, datatype* v, const int n ,
	 const int preSmoothingIterations, const int solvingSmoothingIterations, 
	 const int postSmoothingIterations, const int* nodesPerLevel, const int numberOfLevels,
	 datatype* u_current, datatype* b_current, datatype* v_current){

	/*
		This is a multigrid V-cycle that solves the Poisson eq (It might be the
		 Laplacian, I always get them mixed up). It uses cuda and works on the GPU.

		 The smoother is Jacobi-type (since that's parallel and stuff). Be carefull
		 with that one, since it is very much not robust and can (will) diverge at
		 the first sight of trouble. Yellow-bellied bastard!
	*/

	dim3 block_size_2d(n,n);
	dim3 block_size_1d(n*n);
	//Pre-smoothing.
	for (int iter = 0; iter < preSmoothingIterations/2; iter++){
		//Calculates the error of the guess u and puts it in v. 
		jacobiSmootherLaplacianDev<<<1,block_size_1d>>>(u_current, b_current, n, v_current);

		//If we do it twice we don't have to copy data and stuff.
		jacobiSmootherLaplacianDev<<<1,block_size_1d>>>(v_current, b_current, n, u_current);
	}

	//Calculate error
	//v = f-L*u;
	calculateErrorLaplacian<<<1,block_size_1d>>>(u_current, b_current, v_current, n);

	//Initiate a coarser grid.
	datatype* u_coarse, *b_coarse, *v_coarse;
	u_coarse = firstPonterFromGridLevel(nodesPerLevel, k-1, numberOfLevels, u);
	b_coarse = firstPonterFromGridLevel(nodesPerLevel, k-1, numberOfLevels, b);
	v_coarse = firstPonterFromGridLevel(nodesPerLevel, k-1, numberOfLevels, v);

	//Restrict the grid
	restrictDtoB<<<1, block_size_2d>>>(v_current, n, b_coarse, (n-1)/2+1);

	//Stoping condition to the recursive call.
	if(k==1){
		//"Solves" the coarse system. This should probably be something
		//better than a Jacobi smoother. OPT!
	    for (int iter = 0; iter < solvingSmoothingIterations/2; iter++){
			//Calculates the error of the guess u and puts it in v. 
			jacobiSmootherLaplacianDev<<<1,block_size_1d>>>(u_coarse, b_coarse, (n-1)/2+1, v_coarse);
	
			//If we do it twice we don't have to copy data and stuff.
			jacobiSmootherLaplacianDev<<<1,block_size_1d>>>(v_coarse, b_coarse, (n-1)/2+1, u_coarse);
		}
	}
	else{
		//Recursive call to itself
	    multigridCycleLaplacian(k-1, u,b,v, (n-1)/2+1, preSmoothingIterations, solvingSmoothingIterations, postSmoothingIterations, nodesPerLevel, numberOfLevels, u_coarse, b_coarse, v_coarse);
	}

	interpolate(&tmpGrid, coarseGrid.lengthOfFinerGrid()); //Den här raden e inte klar. FIX!

	//Remove (some of) the error from u_current
	addVectors<<<1,block_size_1d>>>(u_current,v_current, n);

	//Post-smoothing.
	for (int iter = 0; iter < postSmoothingIterations/2; iter++){
		//Calculates the error of the guess u and puts it in v. 
		jacobiSmootherLaplacianDev<<<1,block_size_1d>>>(u_current, b_current, n, v_current);

		//If we do it twice we don't have to copy data and stuff.
		jacobiSmootherLaplacianDev<<<1,block_size_1d>>>(v_current, b_current, n, u_current);
	}

}

int main(){

	const int numberOfLevels = 3;
	const int n = 1>>numberOfLevels;

	//Trust me on this. Im a mathematician. We don't do mistakees.
	const int numberOfDatapoints = numberOfLevels +1>>(numberOfLevels+1)+2; 

	datatype *u, *b, *u_new, *v;

	//Allocate memory
	u = (datatype*) malloc(numberOfDatapoints*sizeof(datatype));
	u_new = (datatype*) malloc(numberOfDatapoints*sizeof(datatype));
	b = (datatype*) malloc(numberOfDatapoints*sizeof(datatype));

	int *nodesPerLevel;
	nodesPerLevel = (int*) malloc(numberOfLevels*sizeof(int));

	datatype *u_dev, *b_dev, *u_dev_new;
	cudaMalloc( (void**) &u_dev, numberOfDatapoints*sizeof(datatype));
	cudaMalloc( (void**) &u_dev_new, numberOfDatapoints*sizeof(datatype));
	cudaMalloc( (void**) &b_dev, numberOfDatapoints*sizeof(datatype));

	int *nodesPerLevel_dev;
	cudaMalloc( (void**) &nodesPerLevel_dev,numberOfLevels*sizeof(int));

	//initiate data
	for(int i=0; i<n*n; i++){
		u[i] = 0.011111111111111;
		b[i] = 1.0f;
	}
	for(int i=0; i<9; i++){
		u_new[i] = i;
	}

	nodesPerLevel[0] = 9*9;
	nodesPerLevel[1] = 5*5;
	nodesPerLevel[2] = 3*3;

	//Copy to GPU
	cudaMemcpy(u_dev, u, numberOfDatapoints*sizeof(datatype), cudaMemcpyHostToDevice);
	cudaMemcpy(u_dev_new, u_new, numberOfDatapoints*sizeof(datatype), cudaMemcpyHostToDevice);
	cudaMemcpy(b_dev, b, numberOfDatapoints*sizeof(datatype), cudaMemcpyHostToDevice);
	cudaMemcpy(nodesPerLevel_dev, nodesPerLevel, numberOfLevels*sizeof(int), cudaMemcpyHostToDevice);

	dim3 block_size(n,n);

	//restrictDtoB<<<1,block_size>>>(u_dev, n, u_dev_new, 3);
	interpolate<<<1,block_size>>>(u_dev_new, 3, u_dev, 5);

	//Copy data from GPU
	cudaMemcpy(u_new, u_dev_new, numberOfDatapoints*sizeof(datatype), cudaMemcpyDeviceToHost);
	cudaMemcpy(u, u_dev, numberOfDatapoints*sizeof(datatype), cudaMemcpyDeviceToHost);
	cudaMemcpy(b, b_dev, numberOfDatapoints*sizeof(datatype), cudaMemcpyDeviceToHost);

	//Print the result
	for(int i=0; i<n*n; i++){
		std::cout<<u[i]<<std::endl;
	}
	std::cout<<std::endl;
	for(int i=0; i<9; i++){
		std::cout<<u_new[i]<<std::endl;
	}

	//Free memory
	free(b);
	free(u);
	free(u_new);
	free(nodesPerLevel);

	cudaFree(u_dev);
	cudaFree(u_dev_new);
	cudaFree(b_dev);

	return 0;
}
