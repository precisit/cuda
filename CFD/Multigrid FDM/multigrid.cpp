#include "grid.h"
#include <cmath>

void multigridCycleLaplacian(int k, Grid *grid, int preSmoothingIterations, int postSmoothingIterations, int solvingSmoothingIterations){
			/*
				This function does a multigrid V-cycle for the Laplacian part of the 
				Stream/Vorticies equation.
			*/

			//Pre-smoothing.
			for (int iter = 0; iter < preSmoothingIterations; iter++){
				grid->jacobiSmootherLaplacian();
			}

			//Calculate error
			//d = f-L*u;
			grid->calculateErrorLaplacian();

			if(k==2){
				datatype sum=0.0f;
				for(int i=0; i<grid->len*grid->len;i++){
					sum += std::abs(grid->getB(i));
				}
				std::cout<<"error: "<<sum<<std::endl;
			}

			std::cout<<std::endl;
			std::cout<<"at first:"<<std::endl;
			grid->print();
			std::cout<<std::endl;
			Grid coarseGrid(grid->lengthOfCoarserGrid());
			grid->restrictDtoB(&coarseGrid);
			//Go to a coarser grid
			//grid->restrict(grid->coarserGrid);
			std::cout<<"restricted:"<<std::endl;
			coarseGrid.printB();

			if(k==1){
				//"Solves" the coarse system. This should probably be something
				//better than a Jacobi smoother. OPT!
			    for (int iter = 0; iter < solvingSmoothingIterations; iter++){
					coarseGrid.jacobiSmootherLaplacian();
				}
			}
			else{

				//Node* tmp;
				//tmp = coarseGrid.B;
				//coarseGrid.B = coarseGrid.ptr;
				//coarseGrid.ptr = tmp;

				//Den här ska ändra saker i d. Gör den det?
			    multigridCycleLaplacian(k-1, &coarseGrid, preSmoothingIterations, postSmoothingIterations, solvingSmoothingIterations);

			}
			
			//HÄR ÄR FELET!
			//Vi kör en liten quick-fix o ser vad som händer. OPT! FIX!
			Grid tmpGrid(coarseGrid.lengthOfFinerGrid());
			coarseGrid.interpolate(&tmpGrid, coarseGrid.lengthOfFinerGrid());
			
			for(int i=0; i<grid->len * grid-> len; i++){
				grid->ptr[i] = grid->ptr[i] + tmpGrid.ptr[i];
			}

			/*
			std::cout<<std::endl;
			std::cout<<"coarse again:"<<std::endl;
			coarseGrid.print();
			std::cout<<std::endl;

			std::cout<<"d:"<<std::endl;
			grid->printB();
			std::cout<<std::endl;
			
			//Add the stuff together.
			//u = u + d;
			grid->addSomeStuffTogether();
			std::cout<<"in the end:"<<std::endl;
			grid->print();
			*/

			//post-smoothing
			for(int iter = 0; iter < postSmoothingIterations; iter++){
			   grid->jacobiSmootherLaplacian();
			}
			
		}


void testingSmoother(const int n){

	Grid grid(n);
	for(int i=0; i<n*n; i++){
		grid.B[i].vort = 1.0f;
	}

	grid.print();

	grid.jacobiSmootherLaplacian();


	grid.print();

}


int main(){

	int n = 5;

	Node* allPointers;
	allPointers = (Node*) malloc(n*n*sizeof(Node));
	
	for(int i=0; i<n*n; i++){
		allPointers[i].vort = 1.0f;//(rand()%100)/100.0f;
		allPointers[i].stream_x = 2.0f;
		allPointers[i].stream_y = 3.0f;
	}

	Grid finestGrid(allPointers, n);
	for(int i=0; i<n*n; i++){
		finestGrid.B[i].vort = 1.0f;
		finestGrid.B[i].stream_x = 1.0f;
		finestGrid.B[i].stream_y = 1.0f;
	}

	for(int i=0; i<n*n; i++){
		for(int j=0; j<n*n; j++){
			std::cout<<finestGrid.getLaplacian(i,j)<<" ";
		}
		std::cout<<std::endl;
	}


	//testingSmoother(n);


	
	for(int i=0; i<10; i++){
		std::cout<<"i: "<<i<<std::endl;
		multigridCycleLaplacian(2, &finestGrid, 10, 10, 20);
		std::cout<<std::endl;
		finestGrid.print();
	}


	/*
	std::cout<<"hoppsan!"<<std::endl<<std::endl;
			for(int y=0; y<n*n; y++){
				for(int x=0; x<n*n; x++){
					if(finestGrid.getVortTranspDisc(x,y) == 0.0f){
						printf("0         ");
					}
					else if(finestGrid.getVortTranspDisc(x,y)>0.0f){
				 		printf ("%3f  ",finestGrid.getVortTranspDisc(x,y));
					}
					else{
						printf ("%3f ", finestGrid.getVortTranspDisc(x,y));
					}
				}
				std::cout<<std::endl;
			}
	*/

	return 0;
};