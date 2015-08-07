#include "adaptive_grid.cpp"

void multigrid(int k, AdaptiveGrid* grid, int pre, int sol, int post, func_t externalFunc){

	//std::cout<<"k: "<<k<<std::endl;

	//Pre-smoothing
	for (int i = 0; i < pre; ++i)
	{
		grid->jacobiSmootherLaplacianStream();
	}

	//Calculate error
	//d = f-L*u;
	grid->calculateErrorLaplacian();

	grid->restrictDtoD(grid->coarserGrid);
	grid->restrictUtoU(grid->coarserGrid);

	grid->coarserGrid->calculateRHS();
	for (int i = 0; i < grid->coarserGrid->len; ++i)
	{
		grid->coarserGrid->w[i].stream = grid->coarserGrid->b[i].stream;
	}

	if(k==1){
		//"Solves" the coarse system. This should probably be something
		//better than a Jacobi smoother. OPT!
	    for (int iter = 0; iter < sol; iter++){
			grid->coarserGrid->jacobiSmootherLaplacianStream();
		}
	}
	else{
	    multigrid(k-1, grid->coarserGrid, pre,sol, post, externalFunc);
	}

	//for (int i = 0; i < grid->len; ++i)
	//{
	//	std::cout<<"d[i]: "<<grid->d[i].x_index<<" ; "<<grid->d[i].y_index<<" , "<<grid->d[i].stream<<std::endl;
	//}
	//std::cout<<std::endl;

	//std::cout<<"k (again): "<<k<<std::endl;

	for (int i = 0; i < grid->coarserGrid->len; ++i)
	{
		grid->coarserGrid->d[i].stream = grid->coarserGrid->u[i].stream - grid->coarserGrid->w[i].stream;
	}

	grid->coarserGrid->interpolateD(grid);

	for (int i = 0; i < grid->len; ++i)
	{
		grid->u[i].stream += grid->d[i].stream;
	}

	/*
	//This part is slow and dumb. OPT!
	grid->resetBoundaryLength();
	grid->setBoundaryLength();
	free(grid->boundaryIndex);
	free(grid->boundaryVals);
	//free(grid->b);
	grid->setBoundary();
	//grid->b = (Node*) calloc(grid->len,sizeof(Node));
	for (int i = 0; i < grid->len; ++i)
	{
		grid->b[i].stream = 0.0f;
	}
	grid->updateBFromBoundary();

	//grid->updateBFromFunction(externalFunc);
	*/

	//post-smoothing
	for(int iter = 0; iter < post; iter++){
	   grid->jacobiSmootherLaplacianStream();
	}
}

void fullMultigridInterpolation(AdaptiveGrid* grid, AdaptiveGrid* fineGrid){
	assert(fineGrid != NULL);
	grid->interpolateU(fineGrid);

}

void fullMultigrid(int k_max, AdaptiveGrid* grid, int pre, int sol, int post, func_t externalFunc){
	//grid->updateBFromFunction(externalFunc);
	for (int iter = 0; iter < sol; iter++){
		grid->jacobiSmootherLaplacianStream();
	}

	for (int k = 1; k <= k_max; ++k)
	{
		//std::cout<<"FMG"<<std::endl;
		AdaptiveGrid *tmpGrid = grid;

		for (int i = 0; i < k; ++i)
		{
			tmpGrid = tmpGrid->finerGrid;
		}
		fullMultigridInterpolation(tmpGrid->coarserGrid, tmpGrid);

		multigrid(k, tmpGrid, pre, sol, post, externalFunc);
		//grid->updateBFromFunction(externalFunc);
	}

}

#ifndef UNITTESTING
int main(int argc, char const *argv[])
{
	/* code */
	return 0;
}
#endif