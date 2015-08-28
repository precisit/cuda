#include "wavelet_decompression.h"

//Interpolate between two nodes
__device__ datatype interpol_dot_d(datatype p1, datatype p2){

 	datatype p;

 	p = (p1 + p2)/2.0;

 	return p; 
 }

//Device function for power of
__device__ int pow_d(int layer, int steplen){

 	int s = 1;

 	for(int i = 1; i<=layer; i++){
 		
 		s = s*steplen;
 	}	

 	return s;
}

//Wavelet decompression. y-index. Set values to nodes where isPicked-value = false through intepolation
__global__ void wavelet_d_kernal_y(Node* matrix, int LEN_OF_MATRIX, int layer, int step_in){
		
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;

	if (idx < LEN_OF_MATRIX && idy < LEN_OF_MATRIX){		

		int step = pow_d(layer, step_in); 

		if((matrix[idx*LEN_OF_MATRIX + idy].isPicked == false)){

			if((matrix[idx*LEN_OF_MATRIX + idy].y_index_global < (LEN_OF_MATRIX -step/2)) && (matrix[idx*LEN_OF_MATRIX + idy].y_index_global >= step/2)){

				datatype p1 = matrix[idx*LEN_OF_MATRIX + idy + step/2].vort;
				datatype p2 = matrix[idx*LEN_OF_MATRIX + idy - step/2].vort;

				if(((matrix[idx*LEN_OF_MATRIX + idy].y_index_global + step/2) % step) == 0 && (((matrix[idx*LEN_OF_MATRIX + idy].x_index_global + step/2) % step) == 0) || ((matrix[idx*LEN_OF_MATRIX + idy].x_index_global) % step) == 0){	

					if((matrix[idx*LEN_OF_MATRIX + idy + step/2].isPicked == true) && (matrix[idx*LEN_OF_MATRIX + idy - step/2].isPicked == true)){						
						
						matrix[idx*LEN_OF_MATRIX + idy].vort = interpol_dot_d(p1, p2);
						matrix[idx*LEN_OF_MATRIX + idy].isPicked = true;				
					
					}
				}
			}
		}
	}

	__threadfence();
	__syncthreads();

}

//Wavelet decompression. x-index. Set values to nodes where isPicked-value = false through intepolation
__global__ void wavelet_d_kernal_x(Node* matrix, int LEN_OF_MATRIX, int layer, int step_in){
	
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;

	if (idx < LEN_OF_MATRIX && idy < LEN_OF_MATRIX){		

		int step = pow_d(layer, step_in); 

		if((matrix[idx*LEN_OF_MATRIX + idy].isPicked == false)){

			if((matrix[idx*LEN_OF_MATRIX + idy].x_index_global < LEN_OF_MATRIX -step/2) && (matrix[idx*LEN_OF_MATRIX + idy].x_index_global >= step/2)){

				datatype p1 = matrix[(idx + step/2)*LEN_OF_MATRIX + idy].vort;
				datatype p2 = matrix[(idx - step/2)*LEN_OF_MATRIX + idy].vort;

				if(((matrix[idx*LEN_OF_MATRIX + idy].x_index_global + step/2) % step) == 0 && (((matrix[idx*LEN_OF_MATRIX + idy].y_index_global - step/2) % step) == 0 || ((matrix[idx*LEN_OF_MATRIX + idy].y_index_global) % step) == 0)){	

					if((matrix[(idx+ step/2)*LEN_OF_MATRIX + idy].isPicked == true) && (matrix[(idx- step/2)*LEN_OF_MATRIX + idy].isPicked == true)){				

						matrix[idx*LEN_OF_MATRIX + idy].vort = interpol_dot_d(p1, p2);
						matrix[idx*LEN_OF_MATRIX + idy].isPicked = true;

					}
				}
			}			
		}
	}

	__threadfence();
	__syncthreads();

}


void wavelet_decompression(Node* array, Node* matrix, Node* new_matrix, int *countTrue){
	
	const int size = LEN_OF_MATRIX*LEN_OF_MATRIX* sizeof(Node);	
	
	Node *d_matrix;	

	//Print values of vort in array
	for(int i = 0; i<*countTrue; i++){

	std::cout<<array[i].vort<<" ";

	}

	//Copy node-values from array to matrix
	for(int i = 0; i<*countTrue; i++){

		for(int x = 0; x<LEN_OF_MATRIX; x++){

			for(int y = 0; y<LEN_OF_MATRIX; y++){

				if(matrix[x*LEN_OF_MATRIX + y].x_index_global == array[i].x_index_global && matrix[x*LEN_OF_MATRIX + y].y_index_global == array[i].y_index_global){

					matrix[x*LEN_OF_MATRIX + y] = array[i];
				}
			}
		}
	}

	//Allocate memory for matrix on GPU
	if (cudaMalloc(&d_matrix, size) != cudaSuccess){
		
		std::cout<< "Can't allocate memory 1!"<<std::endl;
	}

	//Copy matrix to GPU
	if(cudaMemcpy(d_matrix, matrix, size, cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 1!"<<std::endl;
		cudaFree(d_matrix);
	}

	//Set dim of blocks and threads on GPU
	dim3 blockDim(LEN_OF_MATRIX, LEN_OF_MATRIX);
	dim3 gridDim((LEN_OF_MATRIX+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK, (LEN_OF_MATRIX+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	//Itterate through LAYERS and run kernels
	for(int i = LAYERS; i>0; i--){

	wavelet_d_kernal_y<<<gridDim, blockDim>>>(d_matrix, LEN_OF_MATRIX, i, STEP); 

	cudaError_t err = cudaThreadSynchronize();
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;

	wavelet_d_kernal_x<<<gridDim, blockDim>>>(d_matrix, LEN_OF_MATRIX, i, STEP);

	err = cudaThreadSynchronize();
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;

	}

	//Copy matrix back to CPU
	if(cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost) != cudaSuccess){
		cudaError_t err = cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost);
		std::cout<<"Copy to CPU: \n" << cudaGetErrorString(err)<<std::endl;
		delete[] matrix;
		cudaFree(d_matrix);
		std::cout<< "Can't copy back to CPU 1!"<<std::endl;

	}


	//Print vort-values in a matrix
	datatype vort;

	for(int y=LEN_OF_MATRIX-1; y>=0; y--){
    	
    	for(int x=0; x<LEN_OF_MATRIX; x++){  
            
            vort = matrix[x*LEN_OF_MATRIX + y].vort;
            
            if(vort == 1){
                printf("1        ");
            }                       
                    
            else{
                printf ("%f ", vort);
            }
        }

     	std::cout<<std::endl;
    }

    //Print isPicked values in a matrix
	int printIsPicked;
    
    for(int y=LEN_OF_MATRIX-1; y>=0; y--){
    	
    	for(int x=0; x<LEN_OF_MATRIX; x++){            
                    
            printIsPicked = matrix[x*LEN_OF_MATRIX + y].isPicked;
            
            if(printIsPicked == true){
                printf("1   ");
            }                       
                    
            else{
                printf ("0   ");
            }

        }

     	std::cout<<std::endl;
    }

	//Free memory on GPU
	cudaFree(matrix);

	//Free memory on CPU	
	delete[] array;

}