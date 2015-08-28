#include "wavelet_compression.h"

// Find corner- and middle-nodes
__device__ bool corner_node(Node* matrix, int idx, int idy, int LEN_OF_MATRIX){

	if (matrix[idx*LEN_OF_MATRIX + idy].x_index_global == 0 && matrix[idx*LEN_OF_MATRIX + idy].y_index_global == 0){
		return true;
	}

	else if (matrix[idx*LEN_OF_MATRIX + idy].x_index_global == (LEN_OF_MATRIX -1) && matrix[idx*LEN_OF_MATRIX + idy].y_index_global == 0){
		return true;
	}

	else if (matrix[idx*LEN_OF_MATRIX + idy].x_index_global == 0 && matrix[idx*LEN_OF_MATRIX + idy].y_index_global == (LEN_OF_MATRIX -1)){
		return true;
	}

	else if (matrix[idx*LEN_OF_MATRIX + idy].x_index_global == (LEN_OF_MATRIX -1) && matrix[idx*LEN_OF_MATRIX + idy].y_index_global == (LEN_OF_MATRIX - 1)){
		return true;
	}
	else if (matrix[idx*LEN_OF_MATRIX + idy].x_index_global == (LEN_OF_MATRIX/2) && matrix[idx*LEN_OF_MATRIX + idy].y_index_global == (LEN_OF_MATRIX/2)){
		return true;
	}
	else{return false;}
}

// Interpolate between two nodes and compare to TOLERANCE
__device__ bool interpol_dot_c(datatype p1, datatype p2, datatype p3){

 	datatype p;

 	p = (p1 + p2)/2.0;
 	p = p3 - p;
 	
 	if (p<0){

 		p = p*(-1);
 	}	

 	if (p > TOLERANCE){

 		return true;
 	} 
 	else{return false;}
}

// Device-function for power of
__device__ int pow_c(int layer, int steplen){

 	int s = 1;

 	for(int i = 1; i<=layer; i++){
 		
 		s = s*steplen;
 	}	

 	return s;
}

// Wavelet decompression. Itterate through nodes and set true if not inside TOLERANCE
__global__ void wavelet_kernal(Node* matrix, int LEN_OF_MATRIX, int step){


	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;

	if(idx < LEN_OF_MATRIX && idy < LEN_OF_MATRIX){

		// Get corner- and middle-node and save as anchor-nodes
		if (corner_node(matrix, idx, idy, LEN_OF_MATRIX) == true){
			matrix[idx*LEN_OF_MATRIX + idy].isPicked = true;
			matrix[idx*LEN_OF_MATRIX + idy].layer = LAYERS;
			
		}
		
		// Itterate through LAYERS and save importand nodes
		for (int i=1; i<=LAYERS; i++){			

			//x-index
			if((matrix[idx*LEN_OF_MATRIX + idy].x_index_global < LEN_OF_MATRIX -step/2) && (matrix[idx*LEN_OF_MATRIX + idy].x_index_global >= step/2)){

				datatype p1 = matrix[idx*LEN_OF_MATRIX + idy].vort;
				datatype p2 = matrix[(idx + step/2)*LEN_OF_MATRIX + idy].vort;
				datatype p3 = matrix[(idx - step/2)*LEN_OF_MATRIX + idy].vort;
			
				if (((matrix[idx*LEN_OF_MATRIX + idy].x_index_global - step/2) % step) == 0 && (((matrix[idx*LEN_OF_MATRIX + idy].y_index_global - step/2) % step) == 0 || ((matrix[idx*LEN_OF_MATRIX + idy].y_index_global) % step) == 0)){ 		
				
					if (interpol_dot_c(p2, p3, p1) == true){

						matrix[(idx)*LEN_OF_MATRIX + idy].isPicked = true;
													
					}
				}
			}	
			
			//y-index			
			if((matrix[idx*LEN_OF_MATRIX + idy].y_index_global < LEN_OF_MATRIX -step/2) && (matrix[idx*LEN_OF_MATRIX + idy].y_index_global >= step/2)){

				datatype p1 = matrix[idx*LEN_OF_MATRIX + idy].vort;
				datatype p2 = matrix[idx*LEN_OF_MATRIX + idy + step/2].vort;
				datatype p3 = matrix[idx*LEN_OF_MATRIX + idy - step/2].vort;	

				if(((matrix[idx*LEN_OF_MATRIX + idy].y_index_global - step/2) % step) == 0 && (((matrix[idx*LEN_OF_MATRIX + idy].x_index_global - step/2) % step) == 0 || ((matrix[idx*LEN_OF_MATRIX + idy].x_index_global) % step) == 0)){		

					if (interpol_dot_c(p2, p3, p1) == true){

						matrix[(idx)*LEN_OF_MATRIX + idy].isPicked = true;
									
					}
				}
			}
		
			step = step*2;
		
		}
	}

	__syncthreads();
}

//Set layer-value to all nodes where isPicked = true
__global__ void layer_kernel(Node* matrix, int LEN_OF_MATRIX, int step){

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;

	if(idx < LEN_OF_MATRIX && idy < LEN_OF_MATRIX){

		if(matrix[(idx)*LEN_OF_MATRIX + idy].isPicked == true){ 	

			int new_step = pow_c(LAYERS, step);

			for(int i = LAYERS; i>0; i--){			

				if (((matrix[idx*LEN_OF_MATRIX + idy].x_index_global % (new_step/2)) == 0) && ((matrix[idx*LEN_OF_MATRIX + idy].y_index_global % (new_step/2)) == 0)){

					if(matrix[(idx)*LEN_OF_MATRIX + idy].layer == 0){
						matrix[(idx)*LEN_OF_MATRIX + idy].layer = i;
					}										
				}

			new_step = new_step/2;

			}
		}
	}
}

void wavelet_compression(Node* matrix, int* countTrue){	

	Node *d_matrix;	

	const int size = LEN_OF_MATRIX*LEN_OF_MATRIX* sizeof(Node);	

	//Allocate memory on GPU
	if (cudaMalloc(&d_matrix, size) != cudaSuccess){

		std::cout<< "Can't allocate memory 1!"<<std::endl;
		assert(0);
	}

	//Copy matrix to GPU
	if(cudaMemcpy(d_matrix, matrix, size, cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 1!"<<std::endl;
		cudaFree(d_matrix);
		assert(0);
	}

	//Set dim of blocks and threads on GPU
	dim3 blockDim(LEN_OF_MATRIX, LEN_OF_MATRIX);
	dim3 gridDim((LEN_OF_MATRIX+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK, (LEN_OF_MATRIX+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	//Run wavelet_kernel
	wavelet_kernal<<<gridDim, blockDim>>>(d_matrix, LEN_OF_MATRIX, STEP); 

	cudaError_t err = cudaThreadSynchronize();
	
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;

	//Run layer_kernel
	layer_kernel<<<gridDim, blockDim>>>(d_matrix, LEN_OF_MATRIX, STEP);

	err = cudaThreadSynchronize();
	
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;

	//Copy matrix back to CPU
	if(cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost) != cudaSuccess){
		err = cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost);
		std::cout<<"Copy to CPU: \n" << cudaGetErrorString(err)<<std::endl;
		delete[] matrix;
		cudaFree(d_matrix);
		std::cout<< "Can't copy back to CPU 1!"<<std::endl;
		assert(0);

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

	//Count nodes where isPicked-value = true
  	*countTrue = 0;

   	for (int i=0; i<LEN_OF_MATRIX; i++){

	   	for (int j=0; j<LEN_OF_MATRIX; j++){

	   		if (matrix[i*LEN_OF_MATRIX + j].isPicked == true){

	   			*countTrue += 1;
	   		}	    	
	    }
	}


	std::cout<<"countTrue: "<<*countTrue<<std::endl;	
	
	//Free allocated memory on GPU
	cudaFree(d_matrix);	

}