#include <iostream>
#include "node.cpp"

 __device__ bool interpolDot(float p1, float p2, float p3, float p4, float p5, float tol){

 	float p;

 	p = (p1 + p2 + p3 + p4)/4.0;
 	p = p5 - p;
 	
 	if (p<0){

 		p = p*(-1);
 	}	

 	if (p > tol){

 		return true;
 	} 
 	else{return false;}
 }

__global__ void waveletkernal(Node* matrix, int row, int colum, float tol, int step, int layers){

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;


	for (int i=0; i<layers; i++){

		if(idx < row -step && idy < colum -step){
		
			if (idx % 2 == 0 && idy % step == 0){

				float p1 = matrix[idx*colum + idy].vort;
				float p2 = matrix[idx*colum + idy + step].vort;
				float p3 = matrix[(idx + step)*colum + idy].vort;
				float p4 = matrix[(idx + step)*colum + idy + step].vort;
				float p5 = matrix[(idx + step/2)*colum + idy + step/2].vort;
			

				if (interpolDot(p1, p2, p3, p4, p5, tol) == true){

					matrix[(idx + step/2)*colum + idy + step/2].isPicked = true;
					matrix[(idx + step/2)*colum + idy + step/2].layer = i+1;
					
				}		
			}
					
		}
		step = step*2;
	}
}

int main(){

	const int row = 6;
	const int colum = 6;
	int step = 2;
	int layers = 2;	
	int size = row*colum* sizeof(Node);
	float tol = 0.1;
		
	Node* matrix;
	Node *d_matrix;	
	matrix = (Node*) calloc(row*colum,sizeof(Node));	

	int x = 0;	
	int y = 0;

	for(int i=0; i< (row*colum); i++){
		
			matrix[i].x_index = x;
			matrix[i].y_index = y;
			matrix[i].x = rand()%5;
			matrix[i].y = rand()%5;
			matrix[i].vort = rand()%10;
			//matrix[i].isPicked = true;

			if (y<colum - 1){

				y++;
			}
			else{

				y = 0;
				x++;
			}

			//std::cout<< /*"x = "<<matrix[i].x<<std::endl<< "y = "<<matrix[i].y<<std::endl<< */"vort = "<<matrix[i].vort<<std::endl;/*<< "x_index = "<<matrix[i].x_index<<std::endl<< "y_index = "<<matrix[i].y_index<<std::endl<< "isPicked = "<<matrix[i].isPicked<<std::endl;*/
		}
	

	if (cudaMalloc(&d_matrix, size) != cudaSuccess){

		std::cout<< "Can't allocate memory 1!"<<std::endl;
	}

	if(cudaMemcpy(d_matrix, matrix, size, cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 1!"<<std::endl;
		cudaFree(d_matrix);
	}

	dim3 blockDim(row, colum);
	dim3 gridDim(1, 1);

	waveletkernal<<<gridDim, blockDim>>>(d_matrix, row, colum, tol, step, layers); //storlek på de som ska komma tillbaka, vaktor med sparade värden

	cudaError_t err = cudaThreadSynchronize();
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;

	if(cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost) != cudaSuccess){
		//err = cudaMemcpy(d_saveNode, saveNode, colum/step, cudaMemcpyDeviceToHost);
		//std::cout<<"Copy to CPU: \n" << cudaGetErrorString(err)<<std::endl;
		delete[] matrix;
		cudaFree(d_matrix);
		std::cout<< "Can't copy back to CPU 3!"<<std::endl;

	}

	/*for(int i=0; i<row*colum; i++){

		std::cout<<"tja"<< /*matrix[i].x_index<<std::endl<< matrix[i].y_index<<std::endl<< matrix[i].x<<std::endl<< matrix[i].y<<std::endl<< matrix[i].vort<<std::endl<<matrix[i].isPicked<<std::endl;
	
	}*/

	float printVort;
    
    for(int x=0; x<row; x++){
                    
        for(int y=0; y<colum; y++){
            
            printVort = matrix[x + y*colum].vort;
            
            if(printVort == 0.0f){
                printf("0         ");
            }
                        
            else if(printVort>0.0f){
                printf ("%3f  ", printVort);
            }
                        
            else{
                printf ("%3f ", printVort);
            }
        }

     	std::cout<<std::endl;
    }

    float printIsPicked;
    
    for(int x=0; x<row; x++){
                    
        for(int y=0; y<colum; y++){
            
            printIsPicked = matrix[x + y*colum].isPicked;
            
            if(printIsPicked == 1){
                printf("true    ");
            }                       
                    
            else{
                printf ("false   ");
            }
        }

     	std::cout<<std::endl;
    }

   int countTrue = 0;

   for (int i=0; i<row; i++){

	   	for (int j=0; j<colum; j++){

	   		if (matrix[i*colum + j].isPicked == true){

	   			countTrue ++;
	   		}	    	
	    }
	}

	std::cout<<"countTrue: "<<countTrue<<std::endl;

	Node* ordedNodelist;
	ordedNodelist = (Node*) calloc(countTrue,sizeof(Node));
    
    int orderPlace = countTrue -1;

    //FIXA!!!!

	for (int m = 1; m <= layers; m++) {

		for (int i=0; i<row; i++){

	   		for (int j=0; j<colum; j++){
	   			
		    	if (matrix[i*colum + j].isPicked == true && matrix[i*colum + j].layer == m){

		    		ordedNodelist[orderPlace] = matrix[i*colum + j];

		    		orderPlace --;
		    	}
		    }
		}
	}

	for(int i=0; i<countTrue; i++){

		std::cout<<ordedNodelist[i].x_index<<std::endl<< ordedNodelist[i].y_index<<std::endl/*<< ordedNodelist[i].x<<std::endl<< ordedNodelist[i].y<<std::endl<< ordedNodelist[i].vort<<std::endl<<ordedNodelist[i].isPicked*/<<std::endl;
	}


	cudaFree(d_matrix);
	delete[] matrix;	
	return 0;
}