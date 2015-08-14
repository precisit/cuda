#include <assert.h>

#include "node.h"
#include "adaptive_grid.h"
#include "adaptive_multigrid_cuda.h"
#include "wavelet.h"
#include "wavelet_d_cuda.h"
#include "RK4.h"

int main(int argc, char const *argv[])
{
	Node* matrix;
	matrix = (Node*) calloc(row*colum,sizeof(Node));
	Node* nodeArray;
	nodeArray = (Node*) malloc(sizeof(Node));
	int* origoArray;
	//origoArray = (int*) malloc(layers*2*sizeof(int));	
	int* countTrue;
	countTrue = (int*) malloc(1*sizeof(int));

	datatype dt = 0.000001;

	float x = 0;	
	float y = 0;
	
	for(int i=0; i< (row*colum); i++){
		
		matrix[i].x_index_global = x;
		matrix[i].y_index_global = y;
		matrix[i].x = x/(row-1);
		matrix[i].y = y/(colum-1);
			
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

		std::cout<<"x = "<<matrix[i].x<<std::endl<< "y = "<<matrix[i].y<<std::endl<< "vort = "<<matrix[i].vort<<std::endl;/*<< "x_index = "<<matrix[i].x_index<<std::endl<< "y_index = "<<matrix[i].y_index<<std::endl<< "isPicked = "<<matrix[i].isPicked<<std::endl;*/
	}
	

	

	//indata();
	//make_node();

	std::cout<<"fgdfg: "<<nodeArray<<std::endl;

	nodeArray = (Node*) calloc(*countTrue , sizeof(Node));
	origoArray = (int*) malloc(layers*2 * sizeof(int));

	wavelet_compression(matrix, nodeArray, origoArray, countTrue);
	
	

	//*origoArray[layers*2];

	std::cout<<"CountTrue main: "<<*countTrue<<std::endl;

	for(int i = 0; i<layers*2; i++){

		if(i%2 == 0){

			origoArray[i] = row;
		}

		else{

			origoArray[i] = colum;
		}
	}	
	
    std::cout<<"CountTrue main: "<<*countTrue<<std::endl;

    int orderPlace = *countTrue -1;

  std::cout<<"CountTrue main: "<<*countTrue<<std::endl;

	for (int m = 1; m <= layers; m++) {

		for (int i=0; i<row; i++){

	   		for (int j=0; j<colum; j++){
	   			
		    	if (matrix[i*colum + j].isPicked == true && matrix[i*colum + j].layer == m){

		    		nodeArray[orderPlace] = matrix[i*colum + j];

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


	//assert(0);

	std::cout<<"CountTrue main: "<<countTrue<<std::endl;

	for(int i=0; i<*countTrue; i++){

		std::cout<<nodeArray[i].x_index_global<<" "<< nodeArray[i].y_index_global<<" "<<nodeArray[i].layer<<std::endl<<std::endl;
	}

	delete[] matrix;

	std::cout<<"fgdfg: "<<nodeArray<<std::endl;

	RK4(dt, nodeArray, origoArray, *countTrue);

	//adaptive_multigrid(nodeArray, origoArray, *countTrue);
	//assert(0);

	matrix = (Node*) calloc(row*colum,sizeof(Node));

	int x_int = 0;	
	int y_int = 0;
	
	for(int i=0; i< (row*colum); i++){
		
			matrix[i].x_index_global = x_int;
			matrix[i].y_index_global = y_int;				
		
			
		if (y_int<colum - 1){

			y_int++;
		}
		else{

			y_int = 0;
			x_int++;
		}

			//std::cout<< /*"x = "<<matrix[i].x<<std::endl<< "y = "<<matrix[i].y<<std::endl<< "vort = "<<matrix[i].vort<<std::endl;/*<< "x_index = "<<matrix[i].x_index<<std::endl<< "y_index = "<<matrix[i].y_index<<std::endl<< "isPicked = "<<matrix[i].isPicked<<std::endl;*/
	}

	wavelet_decompression(nodeArray, matrix, countTrue);
	//visualize(matrix);

	free (countTrue);
	free (matrix);
	free (origoArray);
	return 0;
};
