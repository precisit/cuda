#include "wavelet_d_cuda.h"

__device__ float interpolDotEdge(float p1, float p2){

 	float p;

 	p = (p1 + p2)/2.0;

 	return p; 
 }

__device__ int pow(int layer, int steplen){

 	int s = 1;

 	for(int i = 1; i<=layer; i++){
 		
 		s = s*steplen;
 	}	

 	return s;
 }

__global__ void true_kernel_y(Node* matrix, int idx, int idy, int step_in, int step, int row, int colum, int layer, int layers){

	float vort;

	if(layer > layers){
		return;
	}	

		int layer_new;
		int step_new;

		layer_new = layer +1;
		step_new = pow(layer_new, step_in);

		dim3 block = (1);
		dim3 grid = (1);
		
		true_kernel_y<<<grid, block>>>(matrix, idx, idy, step_in, step_new, row, colum, layer_new, layers);		

		if(layer == layers-1){

			cudaDeviceSynchronize();
			//true_kernel_y<<<grid, block>>>(matrix, idx, idy, step_in, step_new, row, colum, layer_new, layers);			
		}

		
		if((matrix[idx*colum + idy].isPicked == 0)){
		//if((matrix[idx*colum + idy].isPicked == true)){
			if((matrix[idx*colum + idy+(step/2)].isPicked == 1) && (matrix[idx*colum + idy-(step/2)].isPicked == 1)){	

				if(matrix[idx*colum + idy].y_index_global < (colum -(step/2))){

					if((matrix[idx*colum + idy+(step/2)].y_index_global % step) == 0){ 

					//if((matrix[idx*colum + idy+(step/2)].isPicked == false)){
					//if((matrix[idx*colum + idy+(step/2)].isPicked == 0)){

						
						vort = interpolDotEdge(matrix[idx*colum + idy+(step/2)].vort, matrix[idx*colum + idy-(step/2)].vort);
						atomicExch(&matrix[idx*colum + idy].vort, vort); 
						//matrix[idx*colum + idy+(step/2)].isPicked = true;
						atomicExch(&matrix[idx*colum + idy].isPicked, 1);
					
				
					}
				}
			}
		}
		__threadfence();
		__syncthreads();

	//cudaDeviceSynchronize();	
}



__global__ void true_kernel_x(Node* matrix, int idx, int idy, int step_in, int step, int row, int colum, int layer, int layers){

	float vort;

	if(layer > layers){
		return;
	}	

		int layer_new;
		int step_new;

		layer_new = layer +1;
		step_new = pow(layer_new, step_in);

		dim3 block = (1);
		dim3 grid = (1);
		
		true_kernel_x<<<grid, block>>>(matrix, idx, idy, step_in, step_new, row, colum, layer_new, layers);		

		if(layer == layers-1){

			cudaDeviceSynchronize();
			//true_kernel_x<<<grid, block>>>(matrix, idx, idy, step_in, step_new, row, colum, layer_new, layers);			
		}

		
		if((matrix[idx*colum + idy].isPicked == 0)){		
			if((matrix[(idx+(step/2))*colum + idy].isPicked == 1) && (matrix[(idx-(step/2))*colum + idy].isPicked == 1)){
				if(matrix[idx*colum + idy].x_index_global < (row -(step/2))){
		
					if ((matrix[(idx+(step/2))*colum + idy].x_index_global % step) == 0){			
				
					//if((matrix[(idx+(step/2))*colum + idy].isPicked == false)){
					//if((matrix[(idx+(step/2))*colum + idy].isPicked == 0)){	

						
						vort = interpolDotEdge(matrix[(idx+(step/2))*colum + idy].vort, matrix[(idx-(step/2))*colum + idy].vort);
						atomicExch(&matrix[idx*colum + idy].vort, vort); 
						//matrix[(idx+(step/2))*colum + idy].isPicked = true;	
						atomicExch(&matrix[idx*colum + idy].isPicked, 1);
										
				
					}
				}
			}			
		}			
		__threadfence();
		__syncthreads();
	//cudaDeviceSynchronize();		
}








__global__ void wavelet_d_kernal1(Node* matrix, Node* array, int row, int colum, int countTrue, int layers, int step_in){
		
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;

	for(int i = 0; i<countTrue; i++){

		if(matrix[idx*colum + idy].x_index_global == array[i].x_index_global && matrix[idx*colum + idy].y_index_global == array[i].y_index_global){

			matrix[idx*colum + idy] = array[i];
		}
	}

	free (array);
		
	__syncthreads();
	
	dim3 block = (1);
	dim3 grid = (1);
		
	true_kernel_y<<<grid, block>>>(matrix, idx, idy, step_in, 2, row, colum, 1, layers);
}

__global__ void wavelet_d_kernal2(Node* matrix, Node* array, int row, int colum, int countTrue, int layers, int step_in){
	
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	//__threadfence();
	__syncthreads();
	dim3 block = (1);
	dim3 grid = (1);

	true_kernel_x<<<grid, block>>>(matrix, idx, idy, step_in, 2, row, colum, 1, layers);
	__threadfence();
	
		
	
		/*if(layer == layers-1){

			cudaDeviceSynchronize();
			true_kernel_y<<<grid, block>>>(matrix, idx, idy, step_in, step_new, row, colum, layer_new, layers);			
		}*/	

		/*dim3 block = (1);
		dim3 grid = (1);
		true_kernel_y<<<grid,block>>>(matrix, idx, idy, step_in, 2, row, colum, 1, layers);*/	
	__syncthreads();
}		


//int main(){
void wavelet_decompression(Node* array, Node* matrix, int *countTrue){

	//Hämta in array med node-värden (skapa array)
	//int countTrue = 0;
	//Node *array;
	//int* countTrue;
	//countTrue = (int*) malloc(1*sizeof(int));
	//array = wavelet_start(countTrue, inMatrix);
	Node *d_array;
	const int size = row*colum* sizeof(Node);
	int len = *countTrue * sizeof (Node);
	//int layers = 4;

	std::cout<<sizeof(Node)<<std::endl;
	std::cout<<len<<std::endl;

	Node *d_matrix;	



	/*//const int row = array[0].x_index_global + 1;
	//const int colum = array[0].y_index_global + 1;
	//int size = row*colum* sizeof(Node);
				
	//Node* matrix;
	Node *d_matrix;	
	//matrix = (Node*) calloc(row*colum,sizeof(Node));

	/*int x = 0;	
	int y = 0;
	
	for(int i=0; i< (row*colum); i++){
		
			matrix[i].x_index_global = x;
			matrix[i].y_index_global = y;				
		
			
		if (y<colum - 1){

			y++;
		}
		else{

			y = 0;
			x++;
		}

			//std::cout<< /*"x = "<<matrix[i].x<<std::endl<< "y = "<<matrix[i].y<<std::endl<< "vort = "<<matrix[i].vort<<std::endl;/*<< "x_index = "<<matrix[i].x_index<<std::endl<< "y_index = "<<matrix[i].y_index<<std::endl<< "isPicked = "<<matrix[i].isPicked<<std::endl;*/
	//}*/


	if (cudaMalloc(&d_matrix, size) != cudaSuccess){
		
		std::cout<< "Can't allocate memory 1!"<<std::endl;
	}

	

	if (cudaMalloc(&d_array, len) != cudaSuccess){

		std::cout<< "Can't allocate memory 2!"<<std::endl;
	}

	std::cout<<"Hej"<<std::endl;

	if(cudaMemcpy(d_matrix, matrix, size, cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 1!"<<std::endl;
		cudaFree(d_matrix);
	}

	if(cudaMemcpy(d_array, array, len, cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 2!"<<std::endl;
		cudaFree(d_array);
	}

	
	//dim3 blockDim(*countTrue);
	dim3 blockDim(17/2, 17/2);
	dim3 gridDim(5, 5);

std::cout<<*countTrue<<std::endl;

for(int i = 0; i<*countTrue; i++){

	std::cout<<array[i].vort<<" ";
}


	wavelet_d_kernal1<<<gridDim, blockDim>>>(d_matrix, d_array, row, colum, *countTrue, layers, step); //storlek på de som ska komma tillbaka, vaktor med sparade värden

	cudaError_t err = cudaThreadSynchronize();
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;

	wavelet_d_kernal2<<<gridDim, blockDim>>>(d_matrix, d_array, row, colum, *countTrue, layers, step);

	err = cudaThreadSynchronize();
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;


	if(cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost) != cudaSuccess){
		err = cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost);
		std::cout<<"Copy to CPU: \n" << cudaGetErrorString(err)<<std::endl;
		delete[] matrix;
		cudaFree(d_matrix);
		std::cout<< "Can't copy back to CPU 1!"<<std::endl;

	}

	err = cudaThreadSynchronize();
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;

	cudaDeviceReset();
	cudaThreadExit();

	float vort;

	for(int y=colum-1; y>=0; y--){
    	
    	for(int x=0; x<row; x++){  
            
            vort = matrix[x*colum + y].vort;
            
            if(vort == 1){
                printf("1        ");
            }                       
                    
            else{
                printf ("%f ", vort);
            }
        }

     	std::cout<<std::endl;
    }

	float printIsPicked;
    
    for(int y=colum-1; y>=0; y--){
    	
    	for(int x=0; x<row; x++){            
                    
            printIsPicked = matrix[x*colum + y].isPicked;
            
            if(printIsPicked == 1){
                printf("1   ");
            }                       
                    
            else{
                printf ("0   ");
            }

        }

     	std::cout<<std::endl;
    }





   *countTrue = 0;

   for (int i=0; i<row; i++){

	   	for (int j=0; j<colum; j++){

	   		if (matrix[i*colum + j].isPicked == 1){

	   			*countTrue += 1;
	   		}	    	
	    }
	}


	std::cout<<"countTrue: "<<*countTrue<<std::endl;

	//cudaFree(d_matrix);
	cudaFree(d_array);


	delete[] array;

	//delete[] countTrue;

	//return 0;
	//std::cout<< d_matrix<<std::endl;
	//return d_matrix;
}