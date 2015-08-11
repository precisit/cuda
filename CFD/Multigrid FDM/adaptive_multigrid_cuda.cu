#include "adaptive_grid.cpp"
#include <cuda.h>

__global__ void test2(Node * arr){
	int id = threadIdx.x;
	printf("On thread %d stream=%f \n",id, arr[id].stream);
	arr[id].stream = 0.5;
	printf("Second time: on thread %d stream=%f \n",id, arr[id].stream);
}

__device__ int abs_dev(const int x){
	if(x < 0){
		return -x;
	}
	else{
		return x;
	}
}
__device__ datatype getLaplacianStream(const Node * u, const int index1, const int index2, const datatype h){
	const int dist = abs_dev(u[index1].x_index - u[index2].x_index) +abs_dev(u[index1].y_index - u[index2].y_index);

	if( dist > 1){
		 	return 0.0f;
		 }
		 else{
		 	if( dist == 0 ){
		 		return -4.0f/(h*h);
		 	}
			else{
			 	//if(isInBoundary(u[index2].x_index, u[index2].y_index)){
			 	//	return getFromBoundary(u[index2].x_index, u[index2].y_index)/h;
			 	//}
			 	//else{
			 		return 1.0f/(h*h);
			 	//}
			}
		}

}

__global__ void calculateErrorLaplacian(const Node* b, const Node* u, Node* d, const int len, const datatype h){
	int id = threadIdx.x + blockIdx.x*blockDim.x;

	if(id < len){
		datatype sum = 0.0f;
		for(int j=0; j<len; j++){
			sum += getLaplacianStream(u, j, id, h)*u[j].stream;
		}
		d[id].stream = b[id].stream - sum;
	}

}

__global__ void restrictMat(Node* u_coarse, const int len_coarse, Node* u_fine, const int len_fine ){
	int j = threadIdx.x + blockIdx.x*blockDim.x;

	if( j < len_fine ){

		for (int i = 0; i < len_coarse; ++i)
		{
			if( u_fine[j].x_index_global == u_coarse[i].x_index_global 
					&& u_fine[j].y_index_global == u_coarse[i].y_index_global )
				{
					u_coarse[i].stream = u_fine[j].stream;
				}
		}
	}
}

__global__ void calculateRHS(const Node* u, const Node* d, Node* b, const int len, const datatype h){
	const int id = threadIdx.x + blockIdx.x*blockDim.x;
	if( id < len){

		datatype sum = 0.0f;
		for (int j = 0; j < len; ++j)
		{
			sum += getLaplacianStream(u, j, id, h )*u[j].stream;
		}
		b[id].stream = d[id].stream + sum;
	}
}

__global__ void copy(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		to[i].stream = from[i].stream;
	}
}

__global__ void subtract_gpu(Node* a, const Node* b, const Node* c, const int len){
	//a = b - c
	const int i = threadIdx.x + blockIdx.x*blockDim.x;
	if( i < len ){
		a[i].stream = b[i].stream - c[i].stream;
	}
}

__global__ void add_gpu(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		to[i].stream += from[i].stream;
	}
}


__device__ datatype findNodeIdx(const Node* arr, const int x, const int y, const int n){
	for (int i = 0; i < n; ++i)
	{
		if (arr[i].x_index == x && arr[i].y_index == y)
		{
			return i;
		}
	}
	return -1;
}

__device__ datatype findNodeVal(const Node* arr, const int x, const int y, const int n){
	for (int i = 0; i < n; ++i)
	{
		if (arr[i].x_index == x && arr[i].y_index == y)
		{
			return arr[i].stream;
		}
	}
	printf("			SOMETHING HAS FUCKED UP!");
	return 0.0f;
}

__global__ void interpolate(Node* u_fine, const Node* u_coarse, const int len_fine, const int len_coarse){

	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	//const int thr = threadIdx.x;

	//extern __shared__ datatype valdata[];
	//extern __shared__ int idxdata[];

	if (i < len_fine)
	{

		bool isOnCoarserGrid = false;

		//Take the values from the coarser grid.
		for (int j = 0; j < len_coarse; ++j)
			{
			if( u_coarse[j].x_index_global == u_fine[i].x_index_global && u_coarse[j].y_index_global == u_fine[i].y_index_global )
			{
				u_fine[i].stream = u_coarse[j].stream;
				isOnCoarserGrid = true;
			}
		}

		if (u_fine[i].nodeRight != NULL && isOnCoarserGrid == true)
		{


			datatype tmp_val = findNodeVal(u_fine, u_fine[i].x+2, u_fine[i].y, len_fine);
			int tmp_idx = findNodeIdx(u_fine, u_fine[i].x+1, u_fine[i].y, len_fine);
			u_fine[tmp_idx].stream = (tmp_val + u_fine[i].stream)/2.0f;

		}

		__syncthreads();

		int tmp_idx = -1;
		tmp_idx = findNodeIdx(u_fine, u_fine[i].x, u_fine[i].y-2, len_fine);
		if(tmp_idx != -1 && u_fine[i].y_index % 2 == 0){
			datatype tmp_val = u_fine[tmp_idx].stream;

			tmp_idx = findNodeIdx(u_fine, u_fine[i].x, u_fine[i].y-1, len_fine);

			u_fine[tmp_idx].stream = (u_fine[i].stream+tmp_val)/2.0f;
		}
	}

}

#ifdef OLDCODE
__global__ void interpolate_old( Node* u_fine, Node* u_coarse, const int len_fine, const int len_coarse){
	//This isn't done yet. FIX!
	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < len_fine){


		bool isOnCoarserGrid = false;

		//Take the values from the coarser grid.
		for (int j = 0; j < len_coarse; ++j)
			{
			if( u_coarse[j].x_index_global == u_fine[i].x_index_global && u_coarse[j].y_index_global == u_fine[i].y_index_global )
			{
				u_fine[i].stream = u_coarse[j].stream;
				isOnCoarserGrid = true;
			}
		}

		syncthreads();

		printf("threadIdx: %d, isCoarse: %d, x: %d, y: %d, stream: %f \n", i, isOnCoarserGrid, u_fine[i].x_index, u_fine[i].y_index, u_fine[i].stream);

		if( isOnCoarserGrid == true ){
			Node *middleNode;
			middleNode = &u_fine[i];

			datatype tmp = -7.0f;




			//Interpolate right
			if(middleNode->nodeRight != NULL){
				printf("node right: %p \n", middleNode->nodeRight);
				printf("node right right: %p \n", middleNode->nodeRight->nodeRight);
				//return;
				if(middleNode->nodeRight->nodeRight != NULL){
					//return;
					//middleNode->nodeRight->stream =  6;//middleNode->stream;// + middleNode->nodeRight->nodeRight->stream ) / 2.0f;
					tmp = 0.2;//middleNode->nodeRight->stream;
				}
			}

			printf("tmp = %f \n", tmp);

			//return;

			//Interpolate left
			if(middleNode->nodeLeft != NULL){
				if(middleNode->nodeLeft->nodeLeft != NULL){
					middleNode->nodeLeft->stream =  ( middleNode->stream + middleNode->nodeLeft->nodeLeft->stream ) / 2.0f;
				}
			}

			//Interpolate up
			if(middleNode->nodeAbove != NULL){
				if(middleNode->nodeAbove->nodeAbove != NULL){
					middleNode->nodeAbove->stream =  ( middleNode->stream + middleNode->nodeAbove->nodeAbove->stream ) / 2.0f;
				}
			}

			//Interpolate down
			if(middleNode->nodeBelow != NULL){
				if(middleNode->nodeBelow->nodeBelow != NULL){
					middleNode->nodeBelow->stream =  ( middleNode->stream + middleNode->nodeBelow->nodeBelow->stream ) / 2.0f;
				}
			}

			//And at last, the middle point.
			if(middleNode->nodeBelow != NULL){
				if(middleNode->nodeBelow->nodeRight != NULL){
					middleNode->nodeBelow->nodeRight->stream =  ( middleNode->nodeBelow->stream + middleNode->nodeBelow->nodeRight->nodeRight->stream ) / 2.0f;
				}
			}

		}

		
	}
}
#endif

__global__ void findNeighbours(Node* array, const int len ){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < len){
		for (int j = 0; j < len; ++j)
		{
			//check if below
			if(array[i].x_index_global == array[j].x_index_global && array[i].y_index_global-1 == array[j].y_index_global){
				array[i].nodeBelow = &array[j];
			}
			//check if it's above
			if(array[i].x_index_global == array[j].x_index_global && array[i].y_index_global +1 == array[j].y_index_global){
				array[i].nodeAbove = &array[j];
			}
			//to the right
			if(array[i].x_index_global+1 == array[j].x_index_global && array[i].y_index_global == array[j].y_index_global){
				array[i].nodeRight = &array[j];
			}
			//and to the left
			if(array[i].x_index_global-1 == array[j].x_index_global && array[i].y_index_global == array[j].y_index_global){
				array[i].nodeLeft = &array[j];
			}
		}
	}
	//if(array[i].nodeRight == NULL){
	//	printf("i: %d \n", i);
	//	printf("threadIdx: %d, x: %d, y: %d, stream: %f \n", i, array[i].x_index, array[i].y_index, array[i].stream);
	//}
}

__global__ void jacobiSmootherLaplacianStream(Node *from, Node * to, Node * b, const datatype h, const int len){
	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < len){

		datatype sum = 0.0f;
		for (int j = 0; j < len; ++j)
		{
			if(j != i){
				sum += getLaplacianStream(from, i, j, h) * from[j].stream;
			}
		}
		to[i].stream = (b[i].stream - sum) / getLaplacianStream(from, i,i, h);
	}
}

void ERRORCHECK(){
	cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("Error (!): %s\n", cudaGetErrorString(err));
    exit(-1);
  }
}

__global__ void setVectorsToZero(Node * arr, const int len){
	int i = threadIdx.x + blockIdx.x*blockDim.x;

	if(i < len){
		arr[i].stream = 0.0f;
		//grid->d[i].stream = 0.0f;
		//grid->w[i].stream = 0.0f;
	}
}

void multigrid_gpu(int k, AdaptiveGrid* grid, int pre, int sol, int post, func_t externalFunc){

	//Define the size of the current grid in a CUDA format
	dim3 block_size_1d(grid->len);

	//And the same for the coarser grid.
	dim3 block_size_1d_coarse(grid->coarserGrid->len);

	


	//Pre-smoothing
	for (int i = 0; i < k*pre/2 + 1; ++i)
	{
		jacobiSmootherLaplacianStream<<<1, block_size_1d>>>( grid->u, grid->d, grid->b, grid->h, grid->len);
		jacobiSmootherLaplacianStream<<<1, block_size_1d>>>( grid->d, grid->u, grid->b, grid->h, grid->len);
	}

	//std::cout<<"  1  "<<std::endl;
	ERRORCHECK();

	//Calculate error
	//d = f-L*u;
	calculateErrorLaplacian<<<1, block_size_1d>>>(grid->b, grid->u, grid->d, grid->len, grid->h);

	//std::cout<<"  2  "<<std::endl;
	ERRORCHECK();

	//Restrict d
	restrictMat<<<1, block_size_1d>>>(grid->coarserGrid->d, grid->coarserGrid->len, grid->d , grid->len ); 
	
	//Restrict u
	restrictMat<<<1, block_size_1d>>>(grid->coarserGrid->u, grid->coarserGrid->len, grid->u , grid->len ); 

	//Calculate the right hand side b for the next layer.
	calculateRHS<<<1, block_size_1d_coarse>>>(grid->coarserGrid->u, grid->coarserGrid->d, grid->coarserGrid->b, grid->coarserGrid->len, grid->coarserGrid->h);

	//Copy u to w.
	copy<<<1, block_size_1d>>>(grid->coarserGrid->w, grid->coarserGrid->u, grid->coarserGrid->len);

	//std::cout<<"  3  "<<std::endl;
	ERRORCHECK();


	//Solve the system either using jacobi, or multigrid.
	if(k == 1){
		//"Solves" the coarse system. This should probably be something
		//better than a Jacobi smoother. OPT!
	    for (int iter = 0; iter < sol/2 + 1 ; iter++){
			jacobiSmootherLaplacianStream<<<1, block_size_1d_coarse>>>( grid->coarserGrid->u, grid->coarserGrid->d, grid->coarserGrid->b, grid->h, grid->coarserGrid->len );
			jacobiSmootherLaplacianStream<<<1, block_size_1d_coarse>>>( grid->coarserGrid->d, grid->coarserGrid->u, grid->coarserGrid->b, grid->h, grid->coarserGrid->len );
		}
	}
	else{
		multigrid_gpu(k-1, grid->coarserGrid, pre, sol, post, externalFunc);
	}

	//std::cout<<"  4  "<<std::endl;
	ERRORCHECK();

	subtract_gpu<<<1, block_size_1d_coarse>>>(grid->coarserGrid->d, grid->coarserGrid->u, grid->coarserGrid->w, grid->coarserGrid->len);

	//std::cout<<"  5  "<<std::endl;
	ERRORCHECK();

	interpolate<<<1, block_size_1d>>>( grid->d, grid->coarserGrid->d, grid->len, grid->coarserGrid->len);

	//std::cout<<"  6  "<<std::endl;
	ERRORCHECK();
	//std::cout<<"  7  "<<std::endl;

	add_gpu<<<1, block_size_1d>>>(grid->w, grid->d, grid->len);

	//Post-smoothing
	for (int i = 0; i < k*post/2 + 1; ++i)
	{
		jacobiSmootherLaplacianStream<<<1, block_size_1d>>>( grid->w, grid->u, grid->b, grid->h, grid->len);
		jacobiSmootherLaplacianStream<<<1, block_size_1d>>>( grid->u, grid->w, grid->b, grid->h, grid->len);
	}
	jacobiSmootherLaplacianStream<<<1, block_size_1d>>>( grid->w, grid->u, grid->b, grid->h, grid->len);
}

__global__ void printAll(Node * arr, const int len){
	int i = threadIdx.x + blockDim.x*blockIdx.x;

	if(i<len){
		printf("threadIdx: %d, x: %d, y: %d, stream: %f , ptr: %p \n", i, arr[i].x_index, arr[i].y_index, arr[i].stream, &arr[i]);
	}
}

__global__ void gpuPrint(Node* arr, const int len){
	int i = threadIdx.x + blockDim.x*blockIdx.x;

	if(i<len){
		printf("threadIdx: %d, x: %d, y: %d, stream: %f \n", i, arr[i].x_index, arr[i].y_index, arr[i].stream);
	}
}

void move2gpu(AdaptiveGrid * grid){

	Node *u_dev, *b_dev, *w_dev, *d_dev;

	dim3 block_size_1d(grid->len);

	cudaMalloc( (void**) &u_dev, grid->len*sizeof(Node));
	cudaMalloc( (void**) &b_dev, grid->len*sizeof(Node));
	cudaMalloc( (void**) &w_dev, grid->len*sizeof(Node));
	cudaMalloc( (void**) &d_dev, grid->len*sizeof(Node));

	cudaMemcpy( u_dev, grid->u, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	cudaMemcpy( b_dev, grid->b, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	cudaMemcpy( w_dev, grid->w, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	cudaMemcpy( d_dev, grid->d, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);

	

	free(grid->u);
	free(grid->b);
	free(grid->w);
	free(grid->d);

	grid->u = u_dev;
	grid->w = w_dev;
	grid->d = d_dev;
	grid->b = b_dev;

	findNeighbours<<<1, block_size_1d>>>( grid->u, grid->len );
	findNeighbours<<<1, block_size_1d>>>( grid->b, grid->len );
	findNeighbours<<<1, block_size_1d>>>( grid->w, grid->len );
	findNeighbours<<<1, block_size_1d>>>( grid->d, grid->len );

	ERRORCHECK();

	std::cout<<"Moved data to the gpu."<<std::endl;

}

void move2host(AdaptiveGrid * grid){

	ERRORCHECK();

	Node *u_host, *b_host, *w_host, *d_host;

	dim3 block_size_1d(grid->len);
	//gpuPrint<<<1,block_size_1d>>>(grid->d,grid->len);

	u_host = (Node*) malloc( grid->len*sizeof(Node));
	b_host = (Node*) malloc( grid->len*sizeof(Node));
	w_host = (Node*) malloc( grid->len*sizeof(Node));
	d_host = (Node*) malloc( grid->len*sizeof(Node));

	cudaMemcpy( u_host, grid->u, sizeof(Node)*grid->len, cudaMemcpyDeviceToHost);
	cudaMemcpy( b_host, grid->b, sizeof(Node)*grid->len, cudaMemcpyDeviceToHost);
	cudaMemcpy( w_host, grid->w, sizeof(Node)*grid->len, cudaMemcpyDeviceToHost);
	cudaMemcpy( d_host, grid->d, sizeof(Node)*grid->len, cudaMemcpyDeviceToHost);

	//cudaFree(grid->u);
	//cudaFree(grid->b);
	//cudaFree(grid->w);
	//cudaFree(grid->d);

	grid->u = u_host;
	grid->w = w_host;
	grid->d = d_host;
	grid->b = b_host;
	ERRORCHECK();


	std::cout<<"Moved data from the gpu."<<std::endl;

}

#ifndef UNITTESTING
int main(int argc, char const *argv[])
{

	Node* fromWavelet;
	fromWavelet = (Node*) malloc((5+3+2)*sizeof(Node));

	//datatype h = 0.125f;

	fromWavelet[0] = Node(0.0f, 0.0f, 0,0, 1.0f, 1.0f, 1.0f);
	fromWavelet[1] = Node(0.0f, 1.0f, 0,2, 1.0f, 1.0f, 1.0f);
	fromWavelet[2] = Node(1.0f, 0.0f, 2,0, 1.0f, 1.0f, 1.0f);
	fromWavelet[3] = Node(1.0f, 1.0f, 2,2, 1.0f, 1.0f, 1.0f);
	fromWavelet[4] = Node(0.5f, 0.5f, 1,1, 1.0f, 1.0f, 1.0f);

	AdaptiveGrid grid = AdaptiveGrid(1,3,0,0, &fromWavelet[0], 5, 0.5f); 
	//AdaptiveGrid *dev_grid;
	//cudaMalloc( (void**) &dev_grid, sizeof(AdaptiveGrid));
	//cudaMemcpy(dev_grid, &grid, sizeof(AdaptiveGrid), cudaMemcpyHostToDevice);

	Node *dev_u;
	cudaMalloc( (void**) &(dev_u), 9*sizeof(Node));
	//dev_grid->u = dev_u;
	std::cout<<dev_u<<std::endl;
	grid.u = dev_u;
	std::cout<<grid.u<<std::endl;
	test2<<<1,3>>>(grid.u);

	//test<<<1,4>>>(dev_grid);	

	Node *host_u;
	host_u = (Node*) calloc(9, sizeof(Node));

	cudaMemcpy(host_u, grid.u, sizeof(Node)*9, cudaMemcpyDeviceToHost);
	//assert(0);


	std::cout<<host_u[0].stream<<std::endl;

	cudaFree(dev_u);
	free(fromWavelet);
	return 0;
}
#endif
