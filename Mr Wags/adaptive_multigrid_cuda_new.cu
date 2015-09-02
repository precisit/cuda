#include "adaptive_multigrid_cuda_new.h"

void multigrid_gpu(int k, AdaptiveGrid* grid, int pre, int sol, int post, const int maxGlobIndex, const int maxLayerNr){
	/*
		This is the actual multigrid part of the code. It finds the solution to Poissions eq with -vorticity as the RHS.
	*/

	const int N = grid->len;
	const int N_coarse = grid->coarserGrid->len;

	//Define the size of the current grid in a CUDA format
	dim3 block_size_1d(N);
	dim3 grid_size_1d((N+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	//And the same for the coarser grid.
	dim3 block_size_1d_coarse(N_coarse);
	dim3 grid_size_1d_coarse((N_coarse+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	//Pre-smoothing
	for (int i = 0; i < pre/2 + 1; ++i)
	{
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->u, grid->d, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->d, grid->u, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
	}

	ERRORCHECK();

	//Calculate error
	//d = f-L*u;
	calculateErrorLaplacian<<<grid_size_1d, block_size_1d>>>(grid->b, grid->u, grid->d, grid->len, grid->h, maxGlobIndex, grid->layerNr, maxLayerNr);

	ERRORCHECK();

	//Restrict d
	restrictMat<<<grid_size_1d, block_size_1d>>>(grid->coarserGrid->d, grid->coarserGrid->len, grid->d , grid->len ); 
	
	//Restrict u
	restrictMat<<<grid_size_1d, block_size_1d>>>(grid->coarserGrid->u, grid->coarserGrid->len, grid->u , grid->len ); 

	//Calculate the right hand side b for the next layer.
	//updateBNew<<<grid_size_1d, block_size_1d>>>(grid->u, grid->b, grid->coarserGrid->u, grid->len, grid->coarserGrid->len_coarse, maxGlobIndex);
	//Copy u to w.
	copy<<<grid_size_1d, block_size_1d>>>(grid->coarserGrid->w, grid->coarserGrid->u, grid->coarserGrid->len);

	ERRORCHECK();

	//Solve the system either using jacobi, or multigrid.
	if(k == 1){
		//"Solves" the coarse system. This should probably be something
		//better than a Jacobi smoother. OPT!
	    for (int iter = 0; iter < sol/2 + 1 ; iter++){
			jacobiSmootherLaplacianStream<<<grid_size_1d_coarse, block_size_1d_coarse>>>( grid->coarserGrid->u, grid->coarserGrid->d, grid->coarserGrid->b, grid->h, grid->coarserGrid->len, maxGlobIndex, grid->coarserGrid->layerNr, maxLayerNr );
			jacobiSmootherLaplacianStream<<<grid_size_1d_coarse, block_size_1d_coarse>>>( grid->coarserGrid->d, grid->coarserGrid->u, grid->coarserGrid->b, grid->h, grid->coarserGrid->len, maxGlobIndex, grid->coarserGrid->layerNr, maxLayerNr );
		}
	}
	else{
		multigrid_gpu(k-1, grid->coarserGrid, pre, sol, post, maxGlobIndex, maxLayerNr);
	}

	ERRORCHECK();
	subtract_gpu<<<grid_size_1d_coarse, block_size_1d_coarse>>>(grid->coarserGrid->d, grid->coarserGrid->u, grid->coarserGrid->w, grid->coarserGrid->len);

	ERRORCHECK();
	interpolate<<<grid_size_1d, block_size_1d>>>( grid->d, grid->coarserGrid->d, grid->len, grid->coarserGrid->len, grid->layerNr, maxLayerNr, maxGlobIndex, grid->h, grid->b);

	ERRORCHECK();

	add_gpu<<<grid_size_1d, block_size_1d>>>(grid->w, grid->d, grid->len);

	//Copy u from w
	copy<<<grid_size_1d, block_size_1d>>>(grid->u, grid->w, grid->len);

	//update b
	updateBNew<<<grid_size_1d, block_size_1d>>>(grid->u, grid->b, grid->coarserGrid->u, grid->len, grid->coarserGrid->len, maxGlobIndex, grid->layerNr);

	//Post-smoothing
	for (int i = 0; i < post/2 + 1; ++i)
	{
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->w, grid->u, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->u, grid->w, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
	}
	jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->w, grid->u, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
}


bool isLeaf(Node* u, const int first, const int last, const int x_min, const int x_max, const int y_min, const int y_max){
	for (int i = first; i <= last; ++i)
	{
		if (u[i].x_index_global <= x_max && u[i].x_index_global >= x_min && u[i].y_index_global <= y_max && u[i].y_index_global >= y_min)
		{
			return false;
		}
	}
	return true;
}

void recursiveGridFill(Node* u, int layer, AdaptiveGrid** gridList, int* firstPtr, const int x_min, const int x_max, const int y_min, const int y_max, const float h){
	int diff = (x_max-x_min)/2;

	int *x_loc, *y_loc;
	x_loc = (int*) calloc(1, sizeof(int));
	y_loc = (int*) calloc(1, sizeof(int));

	if (isLeaf(u, firstPtr[LAYERS - layer + 1], firstPtr[LAYERS], x_min, x_min+diff, y_min, y_min+diff) == false)
	{
		gridList[LAYERS - layer + 1]->global2local((x_min+diff/2), (y_min+diff/2),x_loc, y_loc);

		const datatype tmpVal = findClosestNodeValVortGlob(u, (x_min+diff/2), (y_min+diff/2), firstPtr[LAYERS]);

		Node tmp = Node(h*(x_min+diff)/2.0f, h*(y_min+diff)/2.0f, *x_loc, *y_loc,(x_min+diff/2), (y_min+diff/2), tmpVal, 0.0f, 0.0f );
		gridList[LAYERS - layer + 1]->vec.push_back(tmp);
		if (layer != 2)
		{
			recursiveGridFill(u, layer-1, gridList, firstPtr, x_min, x_min+diff, y_min, y_min+diff, h);
		}
	}

	if (isLeaf(u, firstPtr[LAYERS - layer + 1], firstPtr[LAYERS], x_min+diff, x_min+2*diff, y_min, y_min+diff) == false)
	{
		gridList[LAYERS - layer + 1]->global2local((x_min+diff+diff/2), (y_min+diff/2),x_loc, y_loc);

		const datatype tmpVal = findClosestNodeValVortGlob(u, (x_min+diff+diff/2), (y_min+diff/2), firstPtr[LAYERS]);

		Node tmp = Node(h*(x_min+diff)/2.0f, h*(y_min+diff)/2.0f, *x_loc, *y_loc,(x_min+diff+diff/2), (y_min+diff/2), tmpVal, 0.0f, 0.0f );
		gridList[LAYERS - layer + 1]->vec.push_back(tmp);
		if (layer != 2)
		{
			recursiveGridFill(u, layer-1, gridList, firstPtr,x_min+diff, x_min+2*diff, y_min, y_min+diff, h);
		}
	}

	if (isLeaf(u, firstPtr[LAYERS - layer + 1], firstPtr[LAYERS], x_min+diff, x_min+2*diff, y_min+diff, y_min+2*diff) == false)
	{
		gridList[LAYERS - layer + 1]->global2local((x_min+diff+diff/2), (y_min+diff +diff/2),x_loc, y_loc);

		const datatype tmpVal = findClosestNodeValVortGlob(u, (x_min+diff+diff/2), (y_min+diff +diff/2), firstPtr[LAYERS]);

		Node tmp = Node(h*(x_min+diff)/2.0f, h*(y_min+diff)/2.0f, *x_loc, *y_loc,(x_min+diff+diff/2), (y_min+diff +diff/2), tmpVal, 0.0f, 0.0f );
		gridList[LAYERS - layer + 1]->vec.push_back(tmp);
		if (layer != 2)
		{
			recursiveGridFill(u, layer-1, gridList, firstPtr,  x_min+diff, x_min+2*diff, y_min+diff, y_min+2*diff, h);
		}
	}

	if (isLeaf(u, firstPtr[LAYERS - layer + 1], firstPtr[LAYERS], x_min, x_min+diff, y_min+diff, y_min+2*diff) == false)
	{
		gridList[LAYERS - layer + 1]->global2local((x_min+diff/2), (y_min+diff+diff/2),x_loc, y_loc);

		const datatype tmpVal = findClosestNodeValVortGlob(u, (x_min+diff/2), (y_min+diff+diff/2), firstPtr[LAYERS]);

		Node tmp = Node(h*(x_min+diff)/2.0f, h*(y_min+diff)/2.0f, *x_loc, *y_loc,(x_min+diff/2), (y_min+diff+diff/2), tmpVal, 0.0f, 0.0f );
		gridList[LAYERS - layer + 1]->vec.push_back(tmp);
		if (layer != 2)
		{
			recursiveGridFill(u, layer-1, gridList, firstPtr, x_min, x_min+diff, y_min+diff, y_min+2*diff, h);
		}
	}
	free(y_loc);
	free(x_loc);
}

__device__ bool isInBoundaryNew(const Node* u, const int ind){
	return u[ind].nodeRight == NULL || u[ind].nodeLeft == NULL || u[ind].nodeBelow == NULL || u[ind].nodeAbove == NULL;
}

__device__ datatype __BC_Function(const int x, const int y, const int maxGlobIndex){
	//This function doesn't really care if it's a von Neumann or Dirichlet
	//BC. But make sure to avoid mixed BCs.
	if (y == maxGlobIndex)
	{
		return 0.0f; 
	}
	else if(y == 0){
		return 0.0f;
	}
	else{
		return 0.0f;//-1.0f+(2.0f*y)/16.0f;
	}
}

__global__ void updateCoarsestB(const Node* u, Node* b, const int len, const int maxGlobIndex){
	const int id = threadIdx.x + blockDim.x*blockIdx.x;

	if (id < len)
	{
		b[id].stream = 0.0f;
		b[id].stream += -u[id].vort;
		if (isInBoundaryNew(u, id))
		{
			if (u[id].x_index_global == 0 || u[id].x_index_global == maxGlobIndex || u[id].y_index_global == 0 || u[id].y_index_global == maxGlobIndex)
			{
				b[id].stream += KAPPA*__BC_Function(u[id].x_index_global, u[id].y_index_global, maxGlobIndex);
			}
			else{
				printf("updateCoarsestB gives strange results\n");
			}
		}
	}
}

void move2gpu(AdaptiveGrid * grid){

	const int N = grid->len;

	//Define the size of the current grid in a CUDA format
	dim3 block_size_1d(N);
	dim3 grid_size_1d((N+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	Node *u_dev, *b_dev, *w_dev, *d_dev;

	//dim3 block_size_1d(grid->len);

	//std::cout<<"grid len: "<<grid->len<<" \nlen*sizeofNode: "<<grid->len*sizeof(Node)<<std::endl;

	if( cudaMalloc( (void**) &u_dev, grid->len * sizeof(Node)) != cudaSuccess){
		std::cout<< "Can't allocate memory 1!"<<std::endl;
		assert(0);
	}

	//std::cout<<" test 1 \n";
	cudaMalloc( (void**) &b_dev, grid->len*sizeof(Node));
	//std::cout<<" test 2 \n";
	cudaMalloc( (void**) &w_dev, grid->len*sizeof(Node));
	//std::cout<<" test 3 \n";
	cudaMalloc( (void**) &d_dev, grid->len*sizeof(Node));
	//std::cout<<"1\n";

	if(cudaMemcpy( u_dev, grid->u, sizeof(Node)*grid->len, cudaMemcpyHostToDevice) != cudaSuccess){
		std::cout<< "Can't copy memory 1!"<<std::endl;
		assert(0);
	}
	cudaMemcpy( b_dev, grid->b, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	cudaMemcpy( w_dev, grid->w, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	cudaMemcpy( d_dev, grid->d, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	//std::cout<<"2\n";

	free(grid->u);
	free(grid->b);
	free(grid->w);
	free(grid->d);
	//std::cout<<"3\n";

	grid->u = u_dev;
	grid->w = w_dev;
	grid->d = d_dev;
	grid->b = b_dev;


	//std::cout<<"testar saker!\n";
	//u_test<<< 1, 1 >>>(grid->u, 0);



	//std::cout<<"findNeighbours\n";
	ERRORCHECK();

	//u_test<<< 1, 1 >>>( grid->u , 0 );
	findNeighbours<<<grid_size_1d, block_size_1d>>>( grid->u, grid->len );
	findNeighbours<<<grid_size_1d, block_size_1d>>>( grid->b, grid->len );
	findNeighbours<<<grid_size_1d, block_size_1d>>>( grid->w, grid->len );
	findNeighbours<<<grid_size_1d, block_size_1d>>>( grid->d, grid->len );

	ERRORCHECK();
	std::cout<<"Moved data to the gpu."<<std::endl;
}

void move2host(AdaptiveGrid * grid){

	ERRORCHECK();

	Node *u_host, *b_host, *w_host, *d_host;

	dim3 block_size_1d(grid->len);
	//gpuPrint<<<1,block_size_1d>>>(grid->d,grid->len);

	u_host = (Node*) malloc( grid->len*sizeof(Node));
	b_host = (Node*) malloc( grid->len*sizeof(Node));
	w_host = (Node*) malloc( grid->len*sizeof(Node));
	d_host = (Node*) malloc( grid->len*sizeof(Node));

	cudaMemcpy( u_host, grid->u, sizeof(Node)*grid->len, cudaMemcpyDeviceToHost);
	cudaMemcpy( b_host, grid->b, sizeof(Node)*grid->len, cudaMemcpyDeviceToHost);
	cudaMemcpy( w_host, grid->w, sizeof(Node)*grid->len, cudaMemcpyDeviceToHost);
	cudaMemcpy( d_host, grid->d, sizeof(Node)*grid->len, cudaMemcpyDeviceToHost);

	cudaFree(grid->u);
	cudaFree(grid->b);
	cudaFree(grid->w);
	cudaFree(grid->d);

	grid->u = u_host;
	grid->w = w_host;
	grid->d = d_host;
	grid->b = b_host;
	ERRORCHECK();

	grid->findNeighbours(grid->u);
	grid->findNeighbours(grid->w);
	grid->findNeighbours(grid->d);
	grid->findNeighbours(grid->b);
	std::cout<<"Moved data from the gpu."<<std::endl;
}

__global__ void printNorm(Node* u, const int len){
	const int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id == 0)
	{
		datatype sum = 0.0f;
		for (int i = 0; i < len; ++i)
		{
			sum += u[i].stream*u[i].stream;
		}
		printf("Error: %f\n",sum);
	}
}

void adaptive_multigrid_new(Node* array, int* origoArray, int countTrue){

	
	int * pointsInLayer;
	pointsInLayer = (int*) calloc(LAYERS+1, sizeof(int));

	int bigCount = 1;

	int tmpLayer = LAYERS;
	int counter = 0;

	for (int i = 0; i < countTrue; ++i)
	{
		if (array[i].layer == tmpLayer)
		{
			counter++;
		}
		else{
			pointsInLayer[bigCount] = counter;
			counter++;
			tmpLayer = array[i].layer;
			bigCount++;
		}
	}

	pointsInLayer[bigCount] = counter;

	AdaptiveGrid** gridList;
	gridList = (AdaptiveGrid**) malloc(LAYERS*sizeof(AdaptiveGrid*));

	std::cout<<"adaptive_multigrid\n";

	int pointCounter = 0;
	for (int i = 1; i <= LAYERS; ++i)
	{
		AdaptiveGrid* grid = new AdaptiveGrid(i,LAYERS, origoArray[2*(LAYERS-i)],origoArray[2*(LAYERS-i)+1], NULL, 0, 1.0f/(1<<i));
		pointCounter += pointsInLayer[i-1];
		gridList[i-1] = grid;
	}

	std::cout<<"Test1\n";

	for (int i = 0; i < pointsInLayer[1]; ++i)
	{
		gridList[0]->vec.push_back(array[i]);
	}

	std::cout<<"Test2\n";

	recursiveGridFill(array, LAYERS, gridList, pointsInLayer, 0, LEN_OF_MATRIX-1, 0, LEN_OF_MATRIX-1, 1.0f/(1<<LAYERS));

	std::cout<<"Test3\n";

	for (int i = 0; i < LAYERS; ++i)
	{
		std::cout<<"Testigen"<<i<<std::endl;
		gridList[i]->setupFromVector();
	}

	std::cout<<"Test4\n";

	//std::cout<<"move2gpu\n";
	for (int i = 0; i < LAYERS; ++i)
	{
		move2gpu(gridList[i]);
		if(i>0)
			gridList[i]->coarserGrid = gridList[i-1];
		if(i<LAYERS-1)
			gridList[i]->finerGrid = gridList[i+1];


		//const int N = gridList[i]->len;
		//Define the size of the current grid in a CUDA format
		//dim3 block_size_1d(N);
		//dim3 grid_size_1d((N+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

		//updateFromDirichletBC<<< grid_size_1d , block_size_1d >>>( gridList[i]->u, gridList[i]->w,gridList[i]->d,gridList[i]->len, 1<<LAYERS);
	}

	//assert(0);

	updateCoarsestB<<<1,9>>>(gridList[0]->u, gridList[0]->b, gridList[0]->len, LEN_OF_MATRIX-1);



	for (int i = 0; i < 6; ++i)
	{
		std::cout<<"muuuuuultigrid!\n";
		//setVectorsToZero<<<1, grid1.len>>>(grid1.b, grid1.len);
		//dev_updateBFromBoundary<<<1, gridList[0]->len>>>(gridList[0]->b, gridList[0]->u, gridList[0]->len, gridList[0]->layerNr, 4, LEN_OF_MATRIX-1, gridList[0]->h);
		multigrid_gpu(LAYERS-1, gridList[LAYERS-1], 200, 400, 200, LEN_OF_MATRIX-1, LAYERS);
	}
	

	/*
	const int N = gridList[3]->len;

	//Define the size of the current grid in a CUDA format
	dim3 block_size_1d(N);
	dim3 grid_size_1d((N+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);
	updateBNew<<<grid_size_1d, block_size_1d>>>(gridList[3]->u, gridList[3]->b, NULL, gridList[3]->len, 0, LEN_OF_MATRIX -1 , 3);

	for (int i = 0; i < 400; ++i)
	{
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( gridList[3]->u, gridList[3]->d, gridList[3]->b, gridList[3]->h, gridList[3]->len, LEN_OF_MATRIX-1, gridList[3]->layerNr, LAYERS );
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( gridList[3]->d, gridList[3]->u, gridList[3]->b, gridList[3]->h, gridList[3]->len, LEN_OF_MATRIX-1, gridList[3]->layerNr, LAYERS );
	}
	*/

	for (int i = 0; i < LAYERS; ++i)
	{
		calcVortFromStream(gridList[i]);
	}
	//calcVortFromStream(gridList[LAYERS - 1]);

	for (int i = 0; i < LAYERS; ++i)
	{
		move2host(gridList[i]);
	}

	counter = 0;

	Node* tmp;
	for (int j = 0; j < countTrue; ++j)
	{
		for(int lay=LAYERS-1; lay>=0; lay--)
		//for(int lay=0; lay<LAYERS; lay++)
		{
			tmp = gridList[lay]->findGlobNodeGeneral(gridList[lay]->u, 
				array[j].x_index_global, 
				array[j].y_index_global, 
				gridList[lay]->len);

			if (tmp != NULL){
				array[j].vort = tmp->vort;
				//std::cout<<"tmp_vort "<<tmp->vort<<std::endl;
				array[j].stream = tmp->stream;
				break;
			}
			//assert(lay != LAYERS-1);
			assert(lay != 0);
		}
	}
	for (int i = 0; i < LAYERS; ++i)
	{
		free(gridList[i]->u);
		free(gridList[i]->d);
		free(gridList[i]->w);
		free(gridList[i]->b);

		gridList[i]->vec.clear();

		delete gridList[i];
	}

	free(gridList);
	free(pointsInLayer);
	
}

__global__ void gaussSeidelSmootherStream(Node* u, const Node* b, const int len, const int maxGlobIndex, const int layerNr, const datatype h){
	const int id = threadIdx.x + blockDim.x*blockIdx.x;
	if (id==0)
	{
		for (int i = 0; i < len; ++i)
		{
			datatype sum = 0.0f;
			for (int j = 0; j < len; ++j)
			{
				if (i!=j)
				{
					sum +=  getLaplacianStreamNew(u, i, j, h, maxGlobIndex, layerNr, LAYERS) * u[j].stream;
				}
			}
			u[i].stream = (b[i].stream-sum)/getLaplacianStreamNew(u, i, i, h, maxGlobIndex, layerNr, LAYERS);
		}
	}
}

void calcVortFromStream(AdaptiveGrid* grid){
	const int N = grid->len;

	//Define the size of the current grid in a CUDA format
	dim3 block_size_1d(N);
	dim3 grid_size_1d((N+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	//Updates the vorticity and puts the new values in w.

	
	
	vortExterior<<<grid_size_1d, block_size_1d>>>(grid->w, grid->u, grid->len, grid->h, LEN_OF_MATRIX-1);
	copyVortExterior<<<grid_size_1d, block_size_1d>>>(grid->u, grid->w, grid->len);
	vortInterior<<<grid_size_1d, block_size_1d>>>(grid->w, grid->u, grid->len, grid->h);
	
	copyVort<<<grid_size_1d, block_size_1d>>>(grid->u, grid->w, grid->len);
}

__global__ void vortInterior(Node* to, Node* from, const int len, const datatype h){
	const int id = threadIdx.x + blockDim.x*blockIdx.x;

	if (id < len)
	{
		if (isInterior(from, id))
		{
			//This is slow and dumb and shared memory won't help until the memory gets a proper structure. OPT!
			const datatype vort_dx2 = (findNodeValVort(from, from[id].x_index+1, from[id].y_index, len) + findNodeValVort(from, from[id].x_index-1, from[id].y_index, len) - 2.0f * from[id].vort)/(h*h);
			const datatype vort_dy2 = (findNodeValVort(from, from[id].x_index, from[id].y_index+1, len) + findNodeValVort(from, from[id].x_index, from[id].y_index-1, len) - 2.0f * from[id].vort)/(h*h);

			const datatype stream_dx = (findNodeVal(from, from[id].x_index+1, from[id].y_index, len) - findNodeVal(from, from[id].x_index-1, from[id].y_index, len))/(2.0f*h);
			const datatype stream_dy = (findNodeVal(from, from[id].x_index, from[id].y_index+1, len) - findNodeVal(from, from[id].x_index, from[id].y_index-1, len))/(2.0f*h);

			const datatype vort_dx = (findNodeValVort(from, from[id].x_index+1, from[id].y_index, len) - findNodeValVort(from, from[id].x_index-1, from[id].y_index, len))/(2.0f*h);
			const datatype vort_dy = (findNodeValVort(from, from[id].x_index, from[id].y_index+1, len) - findNodeValVort(from, from[id].x_index, from[id].y_index-1, len))/(2.0f*h);

			//__syncthreads();

			to[id].vort = (vort_dy2+vort_dx2)/Re - stream_dy*vort_dx + stream_dx*vort_dy;
			//if(to[id].vort != 0.0f){
			//	printf("HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEJ: %f %d\n", to[id].vort, id);
			//}
			//if (from[id].x_index_global==15 && from[id].y_index_global == 15)
			//{
			//	printf("jovars: %d, %f, d2: %f, %f, h: %f \n",id, to[id].vort, vort_dx2, vort_dy2, h);
			//}

			//printf("vort is %f and h is %f\n",to[id].vort, h);
			//printf("stream is %f and h is %f\n", from[id].stream, h);
		}
	}
}

__global__ void vortExterior(Node* to, Node* from, const int len, const datatype h, const int maxGlobIndex){
	const int id = threadIdx.x + blockDim.x*blockIdx.x;

	if (id < len)
	{
		if (isInterior(from, id) == false)
		{
			if (is1StepFromBoundary(from, id, maxGlobIndex))
			{
				//printf(" ------------- ");
				if (from[id].y_index_global == maxGlobIndex)
				{
					to[id].vort =  2.0f/(h*h)*(from[id].stream - findNodeVal(from, from[id].x_index, from[id].y_index-1, len))-2.0f/h;
					//to[id].vort =  -2.0f/(h*h)*(findNodeVal(from, from[id].x_index, from[id].y_index-1, len))-2.0f/h;
					//printf("DOES THIS EVER HAPPPPPPEN? %f %d %f\n", h, id, to[id].vort);
				}
				else if(from[id].y_index_global == 0){
					to[id].vort =  2.0f/(h*h)*(from[id].stream - findNodeVal(from, from[id].x_index, from[id].y_index+1, len));
					//to[id].vort =  -2.0f/(h*h)*(findNodeVal(from, from[id].x_index, from[id].y_index+1, len));
				}
				else if(from[id].x_index_global == maxGlobIndex){
					to[id].vort =  2.0f/(h*h)*(from[id].stream - findNodeVal(from, from[id].x_index-1, from[id].y_index, len));
					//to[id].vort =  -2.0f/(h*h)*(findNodeVal(from, from[id].x_index-1, from[id].y_index, len));
				}
				else if(from[id].x_index_global == 0){
					to[id].vort =  2.0f/(h*h)*(from[id].stream - findNodeVal(from, from[id].x_index+1, from[id].y_index, len));
					//to[id].vort =  -2.0f/(h*h)*(findNodeVal(from, from[id].x_index+1, from[id].y_index, len));
				}
				else{
					printf("ooooooooooooooookay. Vort exterior is a little bit funky.");
				}
			}
			else{
				//Update according to a closest neigbour.
				//printf("thiiiiiiiiiiiiiiis shouldn't happen. I think. \n");
				if (from[id].nodeRight == NULL)
				{
					to[id].vort = findNodeValVort(from, from[id].x_index-1, from[id].y_index, len);
				}
				else if(from[id].nodeBelow == NULL)
				{
					to[id].vort = findNodeValVort(from, from[id].x_index, from[id].y_index+1, len);
				}
				else if(from[id].nodeLeft == NULL)
				{
					to[id].vort = findNodeValVort(from, from[id].x_index+1, from[id].y_index, len);
				}
				else if (from[id].nodeAbove == NULL)
				{
					to[id].vort = findNodeValVort(from, from[id].x_index, from[id].y_index-1, len);
				}
				else{
					printf("Something has gone wrong in vortExterior!\n");
				}
			}
		}
	}
}

__global__ void copyVort(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		to[i].vort = from[i].vort;
		//printf("From copyVort(), vort: %f, ind: %d \n",to[i].vort, i);
	}
}

__global__ void copyVortExterior(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		if (isInterior(from, i) == false)
		{
			to[i].vort = from[i].vort;
			//printf("From copyVort(), vort: %f, ind: %d \n",to[i].vort, i);
		}
	}
}

__device__ datatype findNodeVal(const Node* arr, const int x, const int y, const int n){
	for (int i = 0; i < n; ++i)
	{
		if (arr[i].x_index == x && arr[i].y_index == y)
		{
			return arr[i].stream;
		}
	}
	printf("			SOMETHING HAS FUCKED UP! (stream edition)!");
	return 0.0f;
}

__device__ datatype findNodeValVort(const Node* arr, const int x, const int y, const int n){
	for (int i = 0; i < n; ++i)
	{
		if (arr[i].x_index == x && arr[i].y_index == y)
		{
			return arr[i].vort;
		}
	}
	printf("			SOMETHING HAS FUCKED UP (vort edition)!");
	return 0.0f;
}

datatype findClosestNodeValVortGlob(const Node* arr, const int x, const int y, const int n){
	int dist = 1000000;
	int dist_best = dist;
	int bestIndex = -1;
	for (int i = 0; i < n; ++i)
	{
		dist = (arr[i].x_index_global-x)*(arr[i].x_index_global-x) + (arr[i].y_index_global-y)*(arr[i].y_index_global-y);
		if (dist < dist_best)
		{
			dist_best = dist;
			bestIndex = i;
			if (dist == 0)
			{
				return arr[i].vort;
			}
		}
	}
	assert(bestIndex != -1);
	return arr[bestIndex].vort;
}

__device__ bool isInterior(const Node* arr, const int i){
	return (arr[i].nodeRight != NULL && arr[i].nodeLeft != NULL && arr[i].nodeBelow != NULL && arr[i].nodeAbove != NULL);
}

__device__ bool is1StepFromBoundary(const Node * u, const int ind, const int maxGlobIndex){
  return u[ind].x_index_global == 0 || u[ind].x_index_global == maxGlobIndex || u[ind].y_index_global == 0 ||
    u[ind].y_index_global == maxGlobIndex;
}

__global__ void jacobiSmootherLaplacianStream(Node *from, Node * to, Node * b, const datatype h, const int len, const int maxGlobIndex, const int layerNr, const int maxLayerNr){
	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < len){

		datatype sum = 0.0f;
		for (int j = 0; j < len; ++j)
		{
			if(j != i){
				sum += getLaplacianStreamNew(from, i, j, h, maxGlobIndex, layerNr, maxLayerNr) * from[j].stream;
			}
		}
		const datatype w = 2.0f/3.0f;
		to[i].stream = w*(b[i].stream - sum) / getLaplacianStreamNew(from, i,i, h , maxGlobIndex, layerNr, maxLayerNr) + from[i].stream*(1.0f-w);
	}
}


__device__ datatype getLaplacianStreamNew(const Node* u, const int index1, const int index2, const datatype h,
	const int maxGlobIndex, const int layerNr, const int maxLayerNr){
	//This function works for the new data structure in which the boundary points 
	//are included in the grid itself. The nodes in the boundary are then heavily weighted 
	//to make sure that they get the correct value.

	const int dist = abs_dev(u[index1].x_index - u[index2].x_index) +abs_dev(u[index1].y_index - u[index2].y_index);

	if( dist > 1){
	 	return 0.0f;
	}
	else{
		if(isInBoundaryNew(u, index1)){
			if (u[index1].x_index_global == 0 || u[index1].x_index_global == maxGlobIndex || u[index1].y_index_global == 0 || u[index1].y_index_global == maxGlobIndex)
			{
				if(isVonNeumannBC(u[index1].x_index_global, u[index1].y_index_global )){
					if (dist == 0)
					{
						return KAPPA/h;
					}
					else if( isOneStepIn(u, index1, index2) )
					{
						return -(u[index1].stream-u[index2].stream)/h;
					}
					else
					{
						return 0.0f;
					}
				}
				else{
					//Only Dirichlev BC should matter here.
					if(	dist == 0){
						return KAPPA;
					}
					else{
						return 0.0f;
					}
				}
			}
			//Only Dirichlev BC should matter here.
			if(	dist == 0){
				return KAPPA;
			}
			else{
				return 0.0f;
			}
		}
		else{
			if( dist == 0 ){
		 		return -4.0f/(h*h);
		 	}
			else{
			 	return 1.0f/(h*h);
			}
		}
	}
}
 
void ERRORCHECK(){
	//cudaThreadSynchronize();
 	//cudaError_t err = cudaGetLastError();
  	//if (err != cudaSuccess) {
    //	printf("Error (!): %s\n", cudaGetErrorString(err));
    //	exit(-1);
  	//}
}

__global__ void calculateErrorLaplacian(const Node* b, const Node* u, Node* d, const int len, const datatype h, const int maxGlobIndex, const int layerNr, const int maxLayerNr){
	int id = threadIdx.x + blockIdx.x*blockDim.x;

	if(id < len){
		datatype sum = 0.0f;
		for(int j=0; j<len; j++){
			sum += getLaplacianStreamNew(u, id, j, h, maxGlobIndex, layerNr, maxLayerNr)*u[j].stream;
		}
		d[id].stream = b[id].stream - sum;
	}

}

__global__ void copy(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		to[i].stream = from[i].stream;
	}
}

__global__ void subtract_gpu(Node* a, const Node* b, const Node* c, const int len){
	//a = b - c
	const int i = threadIdx.x + blockIdx.x*blockDim.x;
	if( i < len ){
		a[i].stream = b[i].stream - c[i].stream;
	}
}

__global__ void findNeighbours(Node* array, const int len ){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	//printf("testar!\n");
	if(i < len){
		for (int j = 0; j < len; ++j)
		{
			//check if below
			if(array[i].x_index_global == array[j].x_index_global && array[i].y_index_global-1 == array[j].y_index_global){
				array[i].nodeBelow = &array[j];
			}
			//check if it's above
			if(array[i].x_index_global == array[j].x_index_global && array[i].y_index_global +1 == array[j].y_index_global){
				array[i].nodeAbove = &array[j];
			}
			//to the right
			if(array[i].x_index_global+1 == array[j].x_index_global && array[i].y_index_global == array[j].y_index_global){
				array[i].nodeRight = &array[j];
			}
			//and to the left
			if(array[i].x_index_global-1 == array[j].x_index_global && array[i].y_index_global == array[j].y_index_global){
				array[i].nodeLeft = &array[j];
			}
		}
	}
}

__device__ int abs_dev(const int x){
	if(x < 0){
		return -x;
	}
	else{
		return x;
	}
}

__device__ bool isVonNeumannBC(const int x, const int y){
	return false; //FIX!
}

__global__ void add_gpu(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		to[i].stream += from[i].stream;
	}
}

__device__ bool isOneStepIn(const Node* u, const int index1, const int index2){
	const int x_diff = u[index1].x_index - u[index2].x_index;
	const int y_diff = u[index1].y_index - u[index2].y_index;

	if (x_diff == -1)
	{
		return u[index1].nodeLeft == NULL;
	}
	else if(x_diff == 1){
		return u[index1].nodeRight == NULL;
	}
	else if(y_diff == -1){
		return u[index1].nodeBelow == NULL;
	}
	else if(y_diff == 1){
		return u[index1].nodeAbove == NULL;
	}
	else{
		printf("Something's wrong. Check isOneStepIn-func!\n");
		return false;
	}
}

__global__ void interpolate(Node* u_fine, const Node* u_coarse, const int len_fine, const int len_coarse, const int layerNr, const int maxLayerNr, const int maxGlobIndex, const datatype h, Node* b){

	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	//const int thr = threadIdx.x;

	//extern __shared__ datatype valdata[];
	//extern __shared__ int idxdata[];

	if (i < len_fine)
	{
		bool isOnCoarserGrid = false;

		//Take the values from the coarser grid.
		for (int j = 0; j < len_coarse; ++j)
			{
			if( u_coarse[j].x_index_global == u_fine[i].x_index_global && u_coarse[j].y_index_global == u_fine[i].y_index_global )
			{
				u_fine[i].stream = u_coarse[j].stream;
				isOnCoarserGrid = true;
			}
		}

		if (u_fine[i].nodeRight != NULL && isOnCoarserGrid == true)
		{
			datatype tmp_val = findNodeVal(u_fine, u_fine[i].x+2, u_fine[i].y, len_fine);
			int tmp_idx = findNodeIdx(u_fine, u_fine[i].x+1, u_fine[i].y, len_fine);
			u_fine[tmp_idx].stream = (tmp_val + u_fine[i].stream)/2.0f;
		}

		__syncthreads();

		int tmp_idx = -1;
		tmp_idx = findNodeIdx(u_fine, u_fine[i].x, u_fine[i].y-2, len_fine);
		if(tmp_idx != -1 && u_fine[i].y_index % 2 == 0){
			datatype tmp_val = u_fine[tmp_idx].stream;

			tmp_idx = findNodeIdx(u_fine, u_fine[i].x, u_fine[i].y-1, len_fine);

			u_fine[tmp_idx].stream = (u_fine[i].stream+tmp_val)/2.0f;
		}
		//updateBFromInterpolation(i, b, u_fine, len_fine, u_coarse, len_coarse, layerNr, maxLayerNr, maxGlobIndex, h);
	}
}

__device__ int findNodeIdx(const Node* arr, const int x, const int y, const int n){
	for (int i = 0; i < n; ++i)
	{
		if (arr[i].x_index == x && arr[i].y_index == y)
		{
			return i;
		}
	}
	return -1;
}

__global__ void restrictMat(Node* u_coarse, const int len_coarse, Node* u_fine, const int len_fine ){
	int j = threadIdx.x + blockIdx.x*blockDim.x;

	if( j < len_fine ){

		for (int i = 0; i < len_coarse; ++i)
		{
			if( u_fine[j].x_index_global == u_coarse[i].x_index_global 
					&& u_fine[j].y_index_global == u_coarse[i].y_index_global )
			{
				u_coarse[i].stream = u_fine[j].stream;
			}
		}
	}
}

__device__ datatype __findNodeGlobStream(const Node* u, const int x, const int y, const int len){

	for (int i = 0; i < len; ++i)
	{
		if (u[i].x_index_global == x && u[i].y_index_global == y)
		{
			return u[i].stream;
		}
	}
	printf("x: %d, y: %d len: %d | findNodeGlobStream gives dumb values. \n", x, y, len );
	return -1;
}

__global__ void updateBNew(const Node* u, Node* b, const Node* u_coarse, const int len, const int len_coarse, const int maxGlobIndex, const int layerNr){
	const int id = threadIdx.x + blockDim.x*blockIdx.x;

	if (id < len)
	{
		b[id].stream = 0.0f;
		b[id].stream += -u[id].vort;
		if (isInBoundaryNew(u, id))
		{
			if (u[id].x_index_global == 0 || u[id].x_index_global == maxGlobIndex || u[id].y_index_global == 0 || u[id].y_index_global == maxGlobIndex)
			{
				b[id].stream += KAPPA*__BC_Function(u[id].x_index_global, u[id].y_index_global, maxGlobIndex);
			}
			else
			{
				int tmp_index = __findNodeGlobIdx(u_coarse, u[id].x_index_global, u[id].y_index_global, len_coarse);
				if (tmp_index == -1)
				{
					if (u[id].nodeRight == NULL || u[id].nodeLeft == NULL)
					{
						int step = (1<<(LAYERS - layerNr));
						b[id].stream += KAPPA*(__findNodeGlobStream(u_coarse, u[id].x_index_global, u[id].y_index_global+step, len_coarse)+
							__findNodeGlobStream(u_coarse, u[id].x_index_global, u[id].y_index_global-step, len_coarse))/2.0f;
					}
					else if (u[id].nodeAbove == NULL || u[id].nodeBelow == NULL)
					{
						int step = (1<<(LAYERS - layerNr));
						b[id].stream += KAPPA*(__findNodeGlobStream(u_coarse, u[id].x_index_global+step, u[id].y_index_global, len_coarse)+
							__findNodeGlobStream(u_coarse, u[id].x_index_global-step, u[id].y_index_global, len_coarse))/2.0f;
					}
					else{
						printf("updateBNew gives weird values.\n");
					}

				}
				else{
					b[id].stream += KAPPA*u_coarse[tmp_index].stream;
				}
			}
		}
	}
}

__device__ int __findNodeGlobIdx(const Node* u, const int x, const int y, const int len){

	for (int i = 0; i < len; ++i)
	{
		if (u[i].x_index_global == x && u[i].y_index_global == y)
		{
			return i;
		}
	}
	return -1;

}

__global__ void updateFromDirichletBC(Node* u, Node* d, Node* w, const int len, const int maxGlobIndex){
	const int id = threadIdx.x + blockDim.x*blockIdx.x;

	if (u[id].x_index_global == 0 || u[id].x_index_global == maxGlobIndex || u[id].y_index_global == 0 || u[id].y_index_global == maxGlobIndex)
	{
		const datatype tmpVal = __BC_Function( u[id].x_index_global, u[id].y_index_global, maxGlobIndex );
		u[id].stream = tmpVal;
	}
}
