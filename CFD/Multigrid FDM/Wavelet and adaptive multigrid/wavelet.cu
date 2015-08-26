//#include "node.h"
//#include "parameters.h"
#include "wavelet.h"


__device__ bool cornerNode(Node* matrix, int idx, int idy, int row, int colum){

	if (matrix[idx*colum + idy].x_index_global == 0 && matrix[idx*colum + idy].y_index_global == 0){
		return true;
	}
	/*else if (matrix[idx*colum + idy].x_index_global == ((row -1)/2) && matrix[idx*colum + idy].y_index_global == 0){
		return true;
	}*/
	else if (matrix[idx*colum + idy].x_index_global == (row -1) && matrix[idx*colum + idy].y_index_global == 0){
		return true;
	}
	/*else if (matrix[idx*colum + idy].x_index_global == 0 && matrix[idx*colum + idy].y_index_global == ((colum -1)/2)){
		return true;
	}*/
	else if (matrix[idx*colum + idy].x_index_global == 0 && matrix[idx*colum + idy].y_index_global == (colum -1)){
		return true;
	}
	/*else if (matrix[idx*colum + idy].x_index_global == (row -1) && matrix[idx*colum + idy].y_index_global == ((colum - 1)/2)){
		return true;
	}*/
	/*else if (matrix[idx*colum + idy].x_index_global == ((row -1)/2) && matrix[idx*colum + idy].y_index_global == (colum - 1)){
		return true;
	}*/
	else if (matrix[idx*colum + idy].x_index_global == (row -1) && matrix[idx*colum + idy].y_index_global == (colum - 1)){
		return true;
	}
	else if (matrix[idx*colum + idy].x_index_global == (row/2) && matrix[idx*colum + idy].y_index_global == (colum/2)){
		return true;
	}
	else{return false;}
}

 __device__ bool interpolDot(float p1, float p2, float p3, float tol){

 	float p;

 	p = (p1 + p2)/2.0;
 	p = p3 - p;
 	
 	if (p<0){

 		p = p*(-1);
 	}	

 	if (p > tol){

 		return true;
 	} 
 	else{return false;}
 }

/*__global__ void wavelet_kernal(Node* matrix, int row, int colum, float tol, int step, int layers, int countTrue){


	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;

	

	if (cornerNode(matrix, idx, idy, row, colum) == true){
		atomicExch(&matrix[idx*colum + idy].isPicked, 1);
		atomicExch(&matrix[idx*colum + idy].layer, layers);
		//matrix[idx*colum + idy].isPicked = 1;
		//matrix[idx*colum + idy].layer = layers;

	}
	__syncthreads();
	for (int i=1; i<=layers; i++){

		int isPicked_x = 0;
		int isPicked_y = 0;
		int layer_x = 0;
		int layer_y = 0; 

		bool x = false;
		bool y = false;

		float p1 = matrix[idx*colum + idy].vort;
		float p2 = matrix[idx*colum + idy + step].vort;
		float p3 = matrix[(idx + step)*colum + idy].vort;
		float p4 = matrix[idx*colum + idy + step/2].vort;
		float p5 = matrix[(idx + step/2)*colum + idy].vort;

		__syncthreads();

		if(matrix[idx*colum + idy].x_index_global < row -step){ 
		
			if ((matrix[idx*colum + idy].x_index_global % step == 0)){ 		
			
				if (interpolDot(p1, p3, p5, tol) == true){

					isPicked_x = 1;
					layer_x = i;
					x = true;
					
					
					//int atomicExch(int* matrix[(idx + step/2)*colum + idy].layer, int 1);
					//matrix[(idx + step/2)*colum + idy].isPicked = true;
					//matrix[(idx + step/2)*colum + idy].layer = i;	
					
				}
			}
		}

		__syncthreads();				
		
		if(matrix[idx*colum + idy].y_index_global < colum -step){

			if(matrix[idx*colum + idy].y_index_global % step == 0){			

				if (interpolDot(p1, p2, p4, tol) == true){

					isPicked_y = 1;
					layer_y = i;
					y = true;

					//atomicExch(&matrix[idx*colum + idy + step/2].isPicked, 1);					

					//matrix[idx*colum + idy + step/2].isPicked = true;
					//matrix[idx*colum + idy + step/2].layer = i;					
					//atomicExch(&matrix[idx*c__syncthreads();olum + idy + step/2].layer, i);
					
				}
			}
		}

		__syncthreads();	
		
		
		step = step*2;

	if(x == true){
	atomicExch(&matrix[(idx + step/2)*colum + idy].isPicked, isPicked_x);
	atomicExch(&matrix[(idx + step/2)*colum + idy].layer, layer_x);
	}
	__syncthreads();
	if(y == true){
	atomicExch(&matrix[idx*colum + idy + step/2].isPicked, isPicked_y);
	atomicExch(&matrix[idx*colum + idy + step/2].layer, layer_y);
	}	
	__syncthreads();
	}
	
					
	

}*/








__global__ void wavelet_kernal(Node* matrix, int row, int colum, float tol, int step, int layers, int countTrue){


	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;

	/*int isPicked = 0;
	int layer= 0;
	bool x = false;
	bool y = false;

	Node* temp;*/

	//cudaMalloc(&temp, row*colum);

	if (cornerNode(matrix, idx, idy, row, colum) == true){
		matrix[idx*colum + idy].isPicked = true;
		matrix[idx*colum + idy].layer = layers;
		//isPicked = 1;
		//layer = layers;
		//atomicExch(&matrix[idx*colum + idy].isPicked, 1);
		//atomicExch(&matrix[idx*colum + idy].layer, layers);
		//x = true;
	}

	

	for (int i=1; i<=layers; i++){			

		float p1 = matrix[idx*colum + idy].vort;
		float p2 = matrix[idx*colum + idy + step].vort;

		float p3 = matrix[(idx + step/2)*colum + idy].vort;
		float p4 = matrix[(idx - step/2)*colum + idy].vort;

		float p5 = matrix[idx*colum + idy + step/2].vort;
		float p6 = matrix[idx*colum + idy - step/2].vort;	

	

		if(matrix[idx*colum + idy].x_index_global < (row -step/2)){ 
		
			if (((matrix[idx*colum + idy].x_index_global - step/2) % step) == 0){ 		
			
				if (interpolDot(p3, p4, p1, tol) == true){

					matrix[(idx)*colum + idy].isPicked = true;
					
					if(matrix[(idx)*colum + idy].layer == 0){
						matrix[(idx)*colum + idy].layer = i;
					}
					/*isPicked = 1;
					layer = i;
					x = true;*/
					
				}
			}
		}

		
		
					
		if(matrix[idx*colum + idy].y_index_global < colum -step/2){

			if(((matrix[idx*colum + idy].y_index_global - step/2) % step) == 0){			

				if (interpolDot(p5, p6, p1, tol) == true){

					matrix[(idx)*colum + idy].isPicked = true;
					
					if(matrix[(idx)*colum + idy].layer == 0){
						matrix[(idx)*colum + idy].layer = i;
					}
					/*isPicked = 1;
					layer = i;
					x = true;*/
						
				}
			}
		}
		
		
		/*if(x == true){
			break;
		}*/
		step = step*2;
	
	}

	//if(x == true){
	//	countTrue = 5;
	//}

	//matrix[(idx)*colum + idy].isPicked = isPicked;
	//matrix[(idx)*colum + idy].layer = layer;

	//atomicExch(&matrix[(idx)*colum + idy].isPicked, isPicked);
	//atomicExch(&matrix[(idx)*colum + idy].layer, layer);						
	
	//__syncthreads();



}

//int main(){
void wavelet_compression(Node* matrix, int* countTrue){				///*, Node* ordedNodelist, int* origoArray*/

	/*const int row = 17;
	const int colum = 17;
	int step = 2;
	int layers = 4;	
	int size = row*colum* sizeof(Node);
	float tol = 0.2;
			
	Node* matrix;*/
	Node *d_matrix;	
	
	//int *d_countTrue;

	const int size = row*colum* sizeof(Node);

	//matrix = (Node*) calloc(row*colum,sizeof(Node));	

	/*int x = 0;	
	int y = 0;
	
	for(int i=0; i< (row*colum); i++){
		
			matrix[i].x_index_global = x;
			matrix[i].y_index_global = y;
			matrix[i].x = rand()%5;
			matrix[i].y = rand()%5;
			
			if (x*x/row > y){
			matrix[i].vort = 2;
			}		

			if (y<colum - 1){

				y++;
			}
			else{

				y = 0;
				x++;
			}

			//std::cout<< /*"x = "<<matrix[i].x<<std::endl<< "y = "<<matrix[i].y<<std::endl<< "vort = "<<matrix[i].vort<<std::endl;/*<< "x_index = "<<matrix[i].x_index<<std::endl<< "y_index = "<<matrix[i].y_index<<std::endl<< "isPicked = "<<matrix[i].isPicked<<std::endl;*/
	//	}
	

	if (cudaMalloc(&d_matrix, size) != cudaSuccess){

		std::cout<< "Can't allocate memory 1!"<<std::endl;
		assert(0);
	}

	/*if (cudaMalloc(&d_countTrue, sizeof(int)) != cudaSuccess){

		std::cout<< "Can't allocate memory 2!"<<std::endl;
		assert(0);
	}*/

	if(cudaMemcpy(d_matrix, matrix, size, cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 1!"<<std::endl;
		cudaFree(d_matrix);
		assert(0);
	}

	/*if(cudaMemcpy(d_countTrue, countTrue, sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){

		std::cout<< "Could not copy to GPU 2!"<<std::endl;
		cudaFree(d_countTrue);
		assert(0);
	}*/

	
	
	dim3 blockDim(row, colum);
	dim3 gridDim(5, 5);

	wavelet_kernal<<<gridDim, blockDim>>>(d_matrix, row, colum, tol, step, layers, *countTrue); //storlek på de som ska komma tillbaka, vaktor med sparade värden

	cudaError_t err = cudaThreadSynchronize();
	
	std::cout<<"Run kernel: \n" << cudaGetErrorString(err)<<std::endl;

	if(cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost) != cudaSuccess){
		err = cudaMemcpy(matrix, d_matrix, size, cudaMemcpyDeviceToHost);
		std::cout<<"Copy to CPU: \n" << cudaGetErrorString(err)<<std::endl;
		delete[] matrix;
		cudaFree(d_matrix);
		std::cout<< "Can't copy back to CPU 1!"<<std::endl;
		assert(0);

	}

	cudaDeviceReset();
	cudaThreadExit();

	/*if(cudaMemcpy(countTrue, d_countTrue, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess){
		//err = cudaMemcpy(countTrue, d_countTrue, sizeof(int), cudaMemcpyDeviceToHost);
		std::cout<<"Copy to CPU: \n" << cudaGetErrorString(err)<<std::endl;
		delete[] countTrue;
		cudaFree(d_countTrue);
		std::cout<< "Can't copy back to CPU 1!"<<std::endl;
		assert(0);

	}*/


	//for(int i=0; i<row*colum; i++){

		//std::cout<<"tja"<< /*matrix[i].x_index<<std::endl<< matrix[i].y_index<<std::endl<< matrix[i].x<<std::endl<< matrix[i].y<<std::endl<< matrix[i].vort<<std::endl<<matrix[i].isPicked<<std::endl;
	
	//}

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

   //Beräkna antal valda noder

  *countTrue = 0;

   for (int i=0; i<row; i++){

	   	for (int j=0; j<colum; j++){

	   		if (matrix[i*colum + j].isPicked == true){

	   			*countTrue += 1;
	   		}	    	
	    }
	}


	std::cout<<"countTrue: "<<*countTrue<<std::endl;

	
	//Skapa lista med valda noder
	//Node* ordedNodelist;
	/*//free(ordedNodelist);
	ordedNodelist = (Node*) calloc(*countTrue , sizeof(Node));

	//*origoArray[layers*2];

	for(int i = 0; i<layers*2; i++){

		if(i%2 == 0){

			origoArray[i] = row;
		}

		else{

			origoArray[i] = colum;
		}
	}	
	
    
    int orderPlace = *countTrue -1;

  

	for (int m = 1; m <= layers; m++) {

		for (int i=0; i<row; i++){

	   		for (int j=0; j<colum; j++){
	   			
		    	if (matrix[i*colum + j].isPicked == true && matrix[i*colum + j].layer == m){

		    		ordedNodelist[orderPlace] = matrix[i*colum + j];

		    		if(matrix[i*colum + j].x_index_global<origoArray[(m-1)*2]){
		    		origoArray[(m-1)*2] = matrix[i*colum + j].x_index_global;
		    		}

		    		if(matrix[i*colum + j].y_index_global<origoArray[(m-1)*2 + 1]){
		    		origoArray[(m-1)*2 + 1] = matrix[i*colum + j].y_index_global;
		    		}

		    		orderPlace --;
		    	}
		    }
		}
	}

	


	for(int i=0; i<*countTrue; i++){

		std::cout<<ordedNodelist[i].x_index_global<<" "<< ordedNodelist[i].y_index_global<<" "<<ordedNodelist[i].layer<<std::endl<<std::endl;
	}

	/*for(int i = 0; i<layers*2; i++){

		std::cout<<origoArray[i]<<std::endl;
	}*/
	//*count_in = countTrue;

	cudaFree(d_matrix);
	//cudaFree(d_countTrue);

	//assert(0);


	//delete[] matrix;	
	//return ordedNodelist;

}