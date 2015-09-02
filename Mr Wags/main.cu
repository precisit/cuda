#include <assert.h>

#include <iostream>
#include <fstream>

#include <iomanip>
#include <sstream>

#include "node.h"
#include "adaptive_grid.h"
#include "adaptive_multigrid_cuda_new.h"
#include "wavelet_compression.h"
#include "wavelet_decompression.h"
#include "RK4.h"







void BubbleSort(Node* num, const int numLength)
{
      int i, j, flag = 1;    // set flag to 1 to start first pass
      Node temp;             // holding variable
      //int numLength = num.length( ); 
      for(i = 1; (i <= numLength) && flag; i++)
     {
          flag = 0;
          for (j=0; j < (numLength -1); j++)
         {
               if (num[j+1].layer > num[j].layer)      // ascending order simply changes to <
              { 
                    temp = num[j];             // swap elements
                    num[j] = num[j+1];
                    num[j+1] = temp;
                    flag = 1;               // indicates that a swap occurred.
               }
          }
     }
     return;   //arrays are passed to functions by address; nothing is returned
}

void ERRORCHECK2(){
	cudaThreadSynchronize();
 	cudaError_t err = cudaGetLastError();
  	if (err != cudaSuccess) {
    	printf("Error (!): %s\n", cudaGetErrorString(err));
    	exit(-1);
  	}
}




int main(int argc, char const *argv[])
{
	Node* matrix;
	matrix = (Node*) calloc(LEN_OF_MATRIX*LEN_OF_MATRIX,sizeof(Node));
	datatype t = 0.0f;

	float x = 0;	
	float y = 0;
	
	for(int i=0; i < (LEN_OF_MATRIX * LEN_OF_MATRIX); i++){
		
		matrix[i].x_index_global = x;
		matrix[i].y_index_global = y;
		matrix[i].x = x/(LEN_OF_MATRIX-1);
		matrix[i].y = y/(LEN_OF_MATRIX-1);
			
		/*if (x*x/LEN_OF_MATRIX> y){
			matrix[i].vort = 1.0f;				
		}*/

		if (y < 0.000001f)
		{
			matrix[i].vort = 1.0f;
		}
			
			if (y<LEN_OF_MATRIX - 1){
				y++;
			}
			else{
				y = 0;
				x++;
			}

		//std::cout<<"x = "<<matrix[i].x<<std::endl<< "y = "<<matrix[i].y<<std::endl<< "vort = "<<matrix[i].vort<<std::endl;/*<< "x_index = "<<matrix[i].x_index<<std::endl<< "y_index = "<<matrix[i].y_index<<std::endl<< "isPicked = "<<matrix[i].isPicked<<std::endl;*/
	}

	Node* nodeArray;
	int* origoArray;
	int* countTrue;


	std::string STRING;



	while(t < END_TIME){

	
	countTrue = (int*) malloc(1*sizeof(int));
	std::cout<<"hej?"<<std::endl;
	ERRORCHECK2();


	wavelet_compression(matrix, countTrue);
	
	nodeArray = (Node*) calloc(*countTrue , sizeof(Node));
	origoArray = (int*) malloc(LAYERS*2 * sizeof(int));

	//*origoArray[layers*2];
	ERRORCHECK2();
	//std::cout<<"CountTrue main: "<<*countTrue<<std::endl;

	for(int i = 0; i<LAYERS*2; i++){
		origoArray[i] = LEN_OF_MATRIX;
	}	
	
    //std::cout<<"CountTrue main: "<<*countTrue<<std::endl;
	ERRORCHECK2();
    int orderPlace = *countTrue -1;

	//std::cout<<"CountTrue main: "<<*countTrue<<std::endl;

	for (int m = 1; m <= LAYERS; m++) {

		for (int i=0; i<LEN_OF_MATRIX; i++){

	   		for (int j=0; j<LEN_OF_MATRIX; j++){
	   			
		    	if (matrix[i*LEN_OF_MATRIX + j].isPicked == true && matrix[i*LEN_OF_MATRIX + j].layer == m){

		    		nodeArray[orderPlace] = matrix[i*LEN_OF_MATRIX + j];

		    		if(matrix[i*LEN_OF_MATRIX + j].x_index_global<origoArray[(m-1)*2]){
		    		origoArray[(m-1)*2] = matrix[i*LEN_OF_MATRIX + j].x_index_global;
		    		}

		    		if(matrix[i*LEN_OF_MATRIX + j].y_index_global<origoArray[(m-1)*2 + 1]){
		    		origoArray[(m-1)*2 + 1] = matrix[i*LEN_OF_MATRIX + j].y_index_global;
		    		}

		    		orderPlace --;
		    	}
		    }
		}
	}

	ERRORCHECK2();

	int stepSize;

	for (int lay = 1; lay <= LAYERS; ++lay)
	{
		//Round origo array down.
		//origoArray[2*i] = 2*(origoArray[2*i]/2);
		//origoArray[2*i+1] = 2*(origoArray[2*i+1]/2);

		stepSize = 1<<(lay);

		if (origoArray[2*(lay-1)] % stepSize != 0)
		{
			origoArray[2*(lay-1)] -= stepSize/2;
		}

		if (origoArray[2*(lay-1)+1] % stepSize != 0)
		{
			origoArray[2*(lay-1)+1] -= stepSize/2;
		}


		//int stepSize = 2;
		//int layerCount = 1;
		//while(nodeArray[i].x_index_global % stepSize == 0 && nodeArray[i].y_index_global % stepSize == 0){
		//	stepSize *= 2;
		//	layerCount++;
		//	if(layerCount == LAYERS){
		//		break;
		//	}
		//}
		//nodeArray[i].layer = layerCount;
	}

	//BubbleSort(nodeArray, *countTrue);




	ERRORCHECK2();



	//assert(0);
	
	std::cout<<"CountTrue main: "<<countTrue<<std::endl;

	for(int i=0; i<*countTrue; i++){

		std::cout<<nodeArray[i].x_index_global<<" "<< nodeArray[i].y_index_global<<" "<<nodeArray[i].layer<<std::endl<<std::endl;
	}

	for(int i=0; i<LAYERS; i++){

		std::cout<<"origo: "<<origoArray[2*i]<<" "<< origoArray[2*i+1]<<std::endl;
	}

	//delete[] matrix;

	ERRORCHECK2();
	for (int i = 0; i < ITERATIONS_UNTIL_GRID_UPDATE; ++i)
	{
		std::cout<<"T=========================================================: "<<t<<std::endl;
		RK4(DELTA_T, nodeArray, origoArray, *countTrue);
	}
	

	//adaptive_multigrid(nodeArray, origoArray, *countTrue);
	//assert(0);

	//matrix = (Node*) calloc(LEN_OF_MATRIX*LEN_OF_MATRIX,sizeof(Node));

	int x_int = 0;	
	int y_int = 0;
	
	for(int i=0; i< (LEN_OF_MATRIX*LEN_OF_MATRIX); i++){
		
			matrix[i] = Node();

			matrix[i].x_index_global = x_int;
			matrix[i].y_index_global = y_int;				
		
			
		if (y_int<LEN_OF_MATRIX - 1){

			y_int++;
		}
		else{

			y_int = 0;
			x_int++;
		}

			//std::cout<< /*"x = "<<matrix[i].x<<std::endl<< "y = "<<matrix[i].y<<std::endl<< "vort = "<<matrix[i].vort<<std::endl;/*<< "x_index = "<<matrix[i].x_index<<std::endl<< "y_index = "<<matrix[i].y_index<<std::endl<< "isPicked = "<<matrix[i].isPicked<<std::endl;*/
	}

	std::cout<<"haj1"<<std::endl;
	ERRORCHECK2();
	std::cout<<"haj2"<<std::endl;
	wavelet_decompression(nodeArray, matrix, countTrue);
	std::cout<<"haj3"<<std::endl;
	ERRORCHECK2();
	std::cout<<"haj4"<<std::endl;
	//visualize(matrix);

	free (countTrue);
	
	free (origoArray);


	t = t + DELTA_T;

	/*
	for (int i =0; i < LEN_OF_MATRIX*LEN_OF_MATRIX; ++i)
	{
		std::ostringstream ss;
		ss << matrix[i].vort;
		STRING += (ss.str())+" ";
	}
	*/
	for (int i =0; i < LEN_OF_MATRIX*LEN_OF_MATRIX; ++i)
	{
		std::ostringstream ss;
		ss << matrix[i].stream;
		STRING += (ss.str())+" ";
	}
	

	/*
	STRING += "A = [";
	for (int x = 0; x < LEN_OF_MATRIX; ++x)
	{
		for (int y = 0; y < LEN_OF_MATRIX; ++y)
		{
			std::ostringstream ss;
			ss << matrix[x*LEN_OF_MATRIX + y].vort;
			STRING += (ss.str())+" ";
		}
			STRING += "\n";
	}
	STRING += "]; \n subplot(1,2,1); surf(A); \n";

	STRING += "B = [";
	for (int x = 0; x < LEN_OF_MATRIX; ++x)
	{
		for (int y = 0; y < LEN_OF_MATRIX; ++y)
		{
			std::ostringstream ss;
			ss << matrix[x*LEN_OF_MATRIX + y].stream;
			STRING += (ss.str())+" ";
		}
			STRING += "\n";
	}
	STRING += "]; \n subplot(1,2,2); surf(B); drawnow; pause(0.5); \n";
	*/


}




    std::ofstream myfile;
  	myfile.open ("MrWagsOutputStream.txt");
  	myfile << STRING;
  	myfile.close();

	free (matrix);
	return 0;
};
