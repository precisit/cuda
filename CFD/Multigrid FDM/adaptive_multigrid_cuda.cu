#include "adaptive_grid.cpp"

__global__ void test2(Node * arr){
	int id = threadIdx.x;
	printf("On thread %d stream=%f \n",id, arr[id].stream);
	arr[id].stream = 0.5;
	printf("Second time: on thread %d stream=%f \n",id, arr[id].stream);
}



__global__ void jacobiSmootherLaplacianStream(){}

__device__ int abs(const int x){
	if(x < 0){
		return -x;
	}
	else{
		return x;
	}
}
__device__ datatype getLaplacianStream(const int index1, const int index2, const datatype h){
	const int dist = abs(u[index1].x_index - u[index2].x_index) +abs(u[index1].y_index - u[index2].y_index);

	if( dist > 1){
		 	return 0.0f;
		 }
		 else{
		 	if( dist == 0 ){
		 		return -4.0f/(h*h);
		 	}
			else{
			 	if(isInBoundary(u[index2].x_index, u[index2].y_index)){
			 		return getFromBoundary(u[index2].x_index, u[index2].y_index)/h;
			 	}
			 	else{
			 		return 1.0f/(h*h);
			 	}
			}
		}

}

__global__ calculateErrorLaplacian(const Node* b, const Node* u, Node* d, const int len, const datatype h){
	int id = threadIdx.x + blockIdx.x;

	if(id < len){
		datatype sum = 0.0f;
		for(int j=0; j<len; j++){
			sum += getLaplacianStream(j, id, h)*u[id].stream;
		}
		d[i].stream = b[i].stream - sum;
	}

}




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
