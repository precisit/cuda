//#include <iostream>
//#include "node.cpp"
//#include "parameters.h"
//#include "wavelet.cu"
//#include <cmath>
#include "wavelet_d_cuda.h"



__device__ float checkInterPol(float p5, float dot){

 	float p;

 	p = p5 - dot;
 	
 	/*if (p<0){

 		p = p*(-1);
 	}*/

 	return p;	

}


__device__ float interpolDotEdge(float p1, float p2){
//float interpolDotEdge(float p1, float p2){

 	float p;

 	p = (p1 + p2)/2.0;

 	return p;
 
 }

__device__ float interpolDotMiddle(float p1, float p2, float p3, float p4){
//float interpolDotMiddle(float p1, float p2, float p3, float p4){

 	float p;

 	p = (p1 + p2 + p3 + p4)/4.0;

 	return p;
 	
 }

__device__ int pow(int layer, int steplen){

 	int s = 1;

 	for(int i = 0; i<layer; i++){
 		
 		s = s*steplen;
 	}	

 	return s;
 }



__global__ void wavelet_d_kernal(Node* matrix, Node* array, int row, int colum, int countTrue, int layers){

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;

	for(int i = 0; i<countTrue; i++){

		if(matrix[idx*colum + idy].x_index_global == array[i].x_index_global && matrix[idx*colum + idy].y_index_global == array[i].y_index_global){

			matrix[idx*colum + idy] = array[i];
		}
	}

	for(int i = layers; i>0; i--){

		int step = pow(i, 2);

		if(matrix[idx*colum + idy].layer = i){
			
			
						
			if(matrix[idx*colum + idy].isPicked == true && (matrix[idx*colum + idy].x_index_global <= (row-step-1)) && (matrix[idx*colum + idy].y_index_global <= (colum-step-1))){

				float p1 = matrix[idx*colum + idy].vort;
				float p2 = matrix[(idx+step)*colum + idy].vort;
				float p3 = matrix[(idx+step)*colum + idy+step].vort;
				float p4 = matrix[idx*colum + idy+step].vort;
				//float p5 = matrix[(idx + step/2)*colum + idy + step/2].vort;
									
				int x1 = matrix[idx*colum + idy].x_index_global;
				int x2 = matrix[(idx+step)*colum + idy].x_index_global;				
				int x3 = matrix[(idx+step)*colum + idy+step].x_index_global;
				int x4 = matrix[idx*colum + idy+step].x_index_global;	
				
				int y1 = matrix[idx*colum + idy].y_index_global;			
				int y2 = matrix[(idx+step)*colum + idy].y_index_global;
				int y3 = matrix[(idx+step)*colum + idy+step].y_index_global;
				int y4 = matrix[idx*colum + idy+step].y_index_global;			
				
				if (matrix[((x1 + x2)/2)*colum + (y1 + y2)/2].isPicked == false){

					matrix[((x1 + x2)/2)*colum + (y1 + y2)/2].vort = interpolDotEdge(p1, p2);
					matrix[((x1 + x2)/2)*colum + (y1 + y2)/2].isPicked = true;
					matrix[((x1 + x2)/2)*colum + (y1 + y2)/2].layer = i-1;
					
				}
				else{

					float p5 = matrix[((x1 + x2)/2)*colum + (y1 + y2)/2].vort;
					float dot = interpolDotEdge(p1, p2);
					float dist = checkInterPol(p5, dot);
					matrix[((x1 + x2)/2)*colum + (y1 + y2)/2].vort = dot + dist;

				}

				if (matrix[((x1 + x4)/2)*colum + (y1 + y4)/2].isPicked == false){

					matrix[((x1 + x4)/2)*colum + (y1 + y4)/2].vort = interpolDotEdge(p1, p4);
					matrix[((x1 + x4)/2)*colum + (y1 + y4)/2].isPicked = true;
					matrix[((x1 + x4)/2)*colum + (y1 + y4)/2].layer = i-1;
					
				}
				else{

					float p5 = matrix[((x1 + x4)/2)*colum + (y1 + y4)/2].vort;
					float dot = interpolDotEdge(p1, p4);
					float dist = checkInterPol(p5, dot);
					matrix[((x1 + x4)/2)*colum + (y1 + y4)/2].vort = dot + dist;					

				}


				if (matrix[((x3 + x4)/2)*colum + (y3 + y4)/2].isPicked == false){

					matrix[((x3 + x4)/2)*colum + (y3 + y4)/2].vort = interpolDotEdge(p3, p4);
					matrix[((x3 + x4)/2)*colum + (y3 + y4)/2].isPicked = true;
					matrix[((x3 + x4)/2)*colum + (y3 + y4)/2].layer = i-1;
					
				}
				else{

					float p5 = matrix[((x3 + x4)/2)*colum + (y3 + y4)/2].vort;
					float dot = interpolDotEdge(p3, p4);
					float dist = checkInterPol(p5, dot);
					matrix[((x3 + x4)/2)*colum + (y3 + y4)/2].vort = dot + dist;
										
				}

				if (matrix[((x3 + x2)/2)*colum + (y3 + y2)/2].isPicked == false){

					matrix[((x3 + x2)/2)*colum + (y3 + y2)/2].vort = interpolDotEdge(p3, p2);
					matrix[((x3 + x2)/2)*colum + (y3 + y2)/2].isPicked = true;
					matrix[((x3 + x2)/2)*colum + (y3 + y2)/2].layer = i-1;
					
				}
				else{

					float p5 = matrix[((x3 + x2)/2)*colum + (y3 + y2)/2].vort;
					float dot = interpolDotEdge(p3, p2);
					float dist = checkInterPol(p5, dot);
					matrix[((x3 + x2)/2)*colum + (y3 + y2)/2].vort = dot + dist;
										
				}


				if (matrix[((x1 + x2 + x3 + x4)/4)*colum + (y1 + y2 + y3 + y4)/4].isPicked == false){

					matrix[((x1 + x2 + x3 + x4)/4)*colum + (y1 + y2 + y3 + y4)/4].vort = interpolDotMiddle(p1, p2, p3, p4);
					matrix[((x1 + x2 + x3 + x4)/4)*colum + (y1 + y2 + y3 + y4)/4].isPicked = true;
					matrix[((x1 + x2 + x3 + x4)/4)*colum + (y1 + y2 + y3 + y4)/4].layer = i-1;
					
					
				}
				else{

					float p5 = matrix[((x1 + x2 + x3 + x4)/4)*colum + (y1 + y2 + y3 + y4)/4].vort;
					float dot = interpolDotMiddle(p1, p2, p3, p4);
					float dist = checkInterPol(p5, dot);
					matrix[((x1 + x2 + x3 + x4)/4)*colum + (y1 + y2 + y3 + y4)/4].vort = dot + dist;
										
				}				
			}
		}	

		__syncthreads();
		
	}
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

	//const int row = array[0].x_index_global + 1;
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
	//}
	

	if (cudaMalloc(&d_matrix, size) != cudaSuccess){

		std::cout<< "Can't allocate memory 1!"<<std::endl;
	}

	if (cudaMalloc(&d_array, len) != cudaSuccess){

		std::cout<< "Can't allocate memory 2!"<<std::endl;
	}

	if(cudaMemcpy(d_matrix, matrix, size, cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 1!"<<std::endl;
		cudaFree(d_matrix);
	}

	if(cudaMemcpy(d_array, array, len, cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 2!"<<std::endl;
		cudaFree(d_array);
	}


	//dim3 blockDim(*countTrue);
	dim3 blockDim(row, colum);
	dim3 gridDim(2, 2);


	wavelet_d_kernal<<<gridDim, blockDim>>>(d_matrix, d_array, row, colum, *countTrue, layers); //storlek på de som ska komma tillbaka, vaktor med sparade värden

	cudaError_t err = cudaThreadSynchronize();
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;



	if(cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost) != cudaSuccess){
		err = cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost);
		std::cout<<"Copy to CPU: \n" << cudaGetErrorString(err)<<std::endl;
		delete[] matrix;
		cudaFree(d_matrix);
		std::cout<< "Can't copy back to CPU 1!"<<std::endl;

	}

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



	//cudaFree(d_matrix);
	cudaFree(d_array);


	delete[] array;

	//delete[] countTrue;

	//return 0;
	//std::cout<< d_matrix<<std::endl;
	//return d_matrix;
}