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

/*
			if(k==2){
				datatype sum=0.0f;
				for(int i=0; i<grid->len*grid->len;i++){
					sum += std::abs(grid->getB(i));
				}
				std::cout<<"error: "<<sum<<std::endl;
			}*/

			//std::cout<<std::endl;
			//std::cout<<"at first:"<<std::endl;
			//grid->print();
			//std::cout<<std::endl;
			Grid coarseGrid(grid->lengthOfCoarserGrid());
			grid->restrictDtoB(&coarseGrid);
			//Go to a coarser grid
			//grid->restrict(grid->coarserGrid);
			//std::cout<<"restricted:"<<std::endl;
			//coarseGrid.printB();

			if(k==1){
				//"Solves" the coarse system. This should probably be something
				//better than a Jacobi smoother. OPT!
			    for (int iter = 0; iter < solvingSmoothingIterations; iter++){
					coarseGrid.jacobiSmootherLaplacian();
				}
			}
			else{
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

			//post-smoothing
			for(int iter = 0; iter < postSmoothingIterations; iter++){
			   grid->jacobiSmootherLaplacian();
			}
			
		};

datatype stream_dx(Grid* grid, int x, int y){

	if(grid->isInternal(x,y) || y==0 || y==grid->len-1){
		return (grid->ptr[(x+1)*grid->len +y].stream_x-grid->ptr[(x-1)*grid->len +y].stream_x)/(2.0f*grid->h);
	}
	else{
		if(x==0){
			return (grid->ptr[(x+1)*grid->len +y].stream_x-grid->ptr[(x)*grid->len +y].stream_x)/(grid->h);
		}
		else if(x== grid->len-1){
			return (grid->ptr[(x)*grid->len +y].stream_x - grid->ptr[(x-1)*grid->len +y].stream_x)/(grid->h);
		}
		else{
			assert(0);
		}
	}

}

datatype vort_dx(Grid* grid, int x, int y){

	if(grid->isInternal(x,y) || y==0 || y==grid->len-1){
		return (grid->ptr[(x+1)*grid->len +y].vort-grid->ptr[(x-1)*grid->len +y].vort)/(2.0f*grid->h);
	}
	else{
		if(x==0){
			return (grid->ptr[(x+1)*grid->len +y].vort-grid->ptr[(x)*grid->len +y].vort)/(grid->h);
		}
		else if(x== grid->len-1){
			return (grid->ptr[(x)*grid->len +y].vort - grid->ptr[(x-1)*grid->len +y].vort)/(grid->h);
		}
		else{
			assert(0);
		}
	}

}

datatype stream_dy(Grid* grid, int x, int y){

	if(grid->isInternal(x,y) || x==0 || x==grid->len-1){
		return (grid->ptr[(x+1)*grid->len +y].stream_x-grid->ptr[(x-1)*grid->len +y].stream_x)/(2.0f*grid->h);
	}
	else{
		if(y==0){
			return (grid->ptr[(x)*grid->len +y+1].stream_x-grid->ptr[(x)*grid->len +y].stream_x)/(grid->h);
		}
		else if(y== grid->len-1){
			return (grid->ptr[(x)*grid->len +y].stream_x - grid->ptr[x*grid->len +y-1].stream_x)/(grid->h);
		}
		else{
			assert(0);
		}
	}

}

datatype vort_dy(Grid* grid, int x, int y){

	if(grid->isInternal(x,y) || x==0 || x==grid->len-1){
		return (grid->ptr[x*grid->len +y+1].vort-grid->ptr[x*grid->len +y-1].vort)/(2.0f*grid->h);
	}
	else{
		if(y==0){
			return (grid->ptr[x*grid->len +y+1].vort-grid->ptr[(x)*grid->len +y].vort)/(grid->h);
		}
		else if(y== grid->len-1){
			return (grid->ptr[(x)*grid->len +y].vort - grid->ptr[x*grid->len +y-1].vort)/(grid->h);
		}
		else{
			assert(0);
		}
	}

}

datatype vort_dLaplace(Grid* grid, int x, int y){

	const int n = grid->len;
	if(grid->isInternal(x,y)){
		return (-4.0f*grid->ptr[x*n + y].vort+grid->ptr[x*n + y+1].vort+grid->ptr[x*n + y-1].vort+grid->ptr[(x+1)*n + y].vort+grid->ptr[(x-1)*n + y].vort)/(grid->h*grid->h);
	}
	else{
		if(x==0){
			assert(0);
		}
		else if(y==0){
			assert(0);
		}
		else if(x==grid->len-1){
			assert(0);
		}
		else if(y==grid->len-1){
			assert(0);
		}
		else{
			assert(0);
		}
	}
}

void calculateVortInterior(Grid* grid){

	const int n = grid->len;
	//Start with the internal points
	for(int x=1; x < n-1; x++){
		for(int y=1; y < n-1; y++){

			grid->D[x*n + y].vort = -(grid->ptr[x*n + y+1].stream_x-grid->ptr[x*n + y-1].stream_x)*
				(grid->ptr[(x+1)*n + y].vort-grid->ptr[(x-1)*n + y].vort)/(2.0f*2.0f*grid->h*grid->h)+
				(grid->ptr[(x+1)*n + y].stream_x-grid->ptr[(x-1)*n + y].stream_x)*
				(grid->ptr[x*n + y+1].vort-grid->ptr[x*n + y-1].vort)/(2.0f*2.0f*grid->h*grid->h)+
				1.0f/(Re*grid->h*grid->h) * (-4.0f*grid->ptr[x*n + y].vort+grid->ptr[x*n + y+1].vort+grid->ptr[x*n + y-1].vort+
				grid->ptr[(x+1)*n + y].vort+grid->ptr[(x-1)*n + y].vort);

		}
	}

}

void calculateVortExterior(Grid* grid){

	const int n = grid->len;
	int x,y;
	const datatype h = grid->h;

	x=0;
	for(y=0; y < n; y++){
		grid->D[x*n + y].vort = -(grid->ptr[0*n + y].stream_x - grid->ptr[1*n + y].stream_x)*2.0f/(h*h);
	}
	x=n-1;
	for(y=0; y < n; y++){
		grid->D[x*n + y].vort = (grid->ptr[(n-2)*n + y].stream_x - grid->ptr[(n-1)*n + y].stream_x)*2.0f/(h*h);
	}
	y=0;
	for(x=0; x < n; x++){
		grid->D[x*n + y].vort = -(grid->ptr[x*n + 0].stream_x - grid->ptr[x*n + 1].stream_x)*2.0f/(h*h);
	}
	y=n-1;
	for(x=0; x < n; x++){
		grid->D[x*n + y].vort = (grid->ptr[x*n + n-2].stream_x - grid->ptr[x*n + n-1].stream_x)*2.0f/(h*h)+1.0*2.0f/h;
	}

}




void ALLOFTHESHITTOGETHER(Grid* grid, datatype dt, datatype t_end){

	datatype t=0.0f;
	while(t < t_end){

		for(int i=0; i<grid->len * grid->len; i++){
			grid->B[i].vort = grid->ptr[i].vort;
		}

		for(int iter = 0; iter < 20; iter++){
			multigridCycleLaplacian(2, grid, 10, 10, 100);
		}

		for(int i=0; i<grid->len * grid->len; i++){
			grid->ptr[i].stream_x = grid->ptr[i].vort;
			grid->ptr[i].vort = grid->B[i].vort;
		}

		//Update boundary vort here
		calculateVortExterior(grid);

		calculateVortInterior(grid);

		for(int i=0; i<grid->len * grid->len; i++){
			grid->ptr[i].vort += dt*grid->D[i].vort;
		}

		std::cout<<"A=[";
		grid->print();
		std::cout<<std::endl<<"];\n mesh(A); drawnow; pause(0.1);";
		std::cout<<std::endl;

		t += dt;

	}

}

int main(){

	int n = 9;

	Node* allPointers;
	allPointers = (Node*) malloc(n*n*sizeof(Node));
	
	for(int i=0; i<n*n; i++){
		allPointers[i].vort = 0.0f;
	}

	Grid finestGrid(allPointers, n);
	/*for(int i=0; i<n*n; i++){
		finestGrid.B[i].vort = 0.0f;
		finestGrid.B[i].stream_x = 1.0f;
		finestGrid.B[i].stream_y = 1.0f;
	}*/


	ALLOFTHESHITTOGETHER(&finestGrid, 0.02f, 1.0f);

	/*
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
	*/


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