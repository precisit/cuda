#include "adaptive_multigrid_cuda.h"

typedef datatype (*func_def)(int, int);

__device__ int abs_dev(const int x){
	if(x < 0){
		return -x;
	}
	else{
		return x;
	}
}

__device__ bool is1StepFromBoundary(const Node * u, const int ind, const int maxGlobIndex){
  return u[ind].x_index_global == 0 || u[ind].x_index_global == maxGlobIndex || u[ind].y_index_global == 0 ||
    u[ind].y_index_global == maxGlobIndex;
}

__device__ bool isInCorner(const Node* u, const int ind, const int maxGlobIndex){

  return (u[ind].x_index_global == 0 || u[ind].x_index_global == maxGlobIndex) && (u[ind].y_index_global == 0 ||
    u[ind].y_index_global == maxGlobIndex);

}

__device__ bool is2StepsFromBoundary(const Node * u, const int ind, const int maxGlobIndex){
  return u[ind].x_index_global == 1 || u[ind].x_index_global == maxGlobIndex-1 || u[ind].y_index_global == 1 ||
    u[ind].y_index_global == maxGlobIndex-1;
}

__device__ datatype dev_invPow(datatype x, const int n){
  //x^(-n)
  if (n==0) {
    return 1.0f;
  }
  else{
    datatype out = x;
    for(int i=0; i<n-1; ++i){
      out = out * x;
    }
    return 1.0f/out;
  }
}

__device__ datatype getLaplacianStream(const Node * u, const int index1, const int index2, const datatype h,
  const int maxGlobIndex, const int layerNr, const int maxLayerNr){
    //We can probably get some massive speed-up by reordering the if clauses
    //in this function. OPT!

	const int dist = abs_dev(u[index1].x_index - u[index2].x_index) +abs_dev(u[index1].y_index - u[index2].y_index);

	if( dist > 1){
	 	return 0.0f;
	}
    else{
	    if(is1StepFromBoundary(u, index1, maxGlobIndex) || is1StepFromBoundary(u, index2, maxGlobIndex) ){
	        if (dist == 0) {
	    		if (isInCorner(u, index1, maxGlobIndex)) {
		            //beta+beta
		            int n = maxLayerNr-layerNr;
		            datatype tmp = dev_invPow(2.0f, n);

		            return -(1.0f+tmp)*4.0f/( (h*h) * (tmp+tmp*tmp) );
	          	}
	        	else{
		            //-2/h^2 + beta
		            int n = maxLayerNr-layerNr;
		            datatype tmp = dev_invPow(2.0f, n);
		            return -2.0f/(h*h) + -(1.0f+tmp)*2.0f/( (h*h) * (tmp+tmp*tmp) );
		       	}
	        }
	        else if ( is2StepsFromBoundary(u, index1, maxGlobIndex) || is2StepsFromBoundary(u, index2, maxGlobIndex)  ) {
	        	//alfa
	        	int n = maxLayerNr-layerNr;
	        	return 2.0f/( (h*h) * (1+dev_invPow(2.0f, n)) );
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

__device__ bool isInBoundaryNew(const Node* u, const int ind){
	return u[ind].nodeRight == NULL || u[ind].nodeLeft == NULL || u[ind].nodeLeft == NULL || u[ind].nodeAbove == NULL;
}

__device bool isVonNeumannBC(const int x, const int y){
	return false; //FIX!
}

__device__ bool isOneStepIn(u, index1, index2){
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
				if(isVonNeumannBC(u[index1].x_index_global, u[index1].y_index_global ){
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
			else{
				int tmp_index = __findNodeGlobIdx(u_coarse, u[id].x_index_global, u[id].y_index_global, len_coarse);
				if (tmp_index == -1)
				{
					if (u[id].nodeRight == NULL || u[id].nodeLeft == NULL)
					{
						int step = 1<<(LAYERS - layerNr);
						b[id].stream += KAPPA*(__findNodeGlobStream(u_coarse, u[id].x_index_global, u[id].y_index_global+step, len_coarse)+
							__findNodeGlobStream(u_coarse, u[id].x_index_global, u[id].y_index_global-step, len_coarse))/2.0f;
					}
					else if (u[id].nodeAbove == NULL || u[id].nodeBelow == NULL)
					{
						int step = 1<<(LAYERS - layerNr);
						b[id].stream += KAPPA*(__findNodeGlobStream(u_coarse, u[id].x_index_global+step, u[id].y_index_global, len_coarse)+
							__findNodeGlobStream(u_coarse, u[id].x_index_global-step, u[id].y_index_global, len_coarse))/2.0f;
					}
					else{
						printf("updateBNew gives weird values.\n");
					}

				}
				else{
					b[id].stream += KAPPA*u_coarse[tmp_index];
				}
			}
		}
	}
}

__device__ datatype __BC_Function(const int x, const int y, const int maxGlobIndex){

	//This function doesn't really care if it's a von Neumann or Dirichlet
	//BC. But make sure to avoid mixed BCs.
	if (y == maxGlobIndex)
	{
		return 1.0f;
	}
	else if(y == 0){
		return -1.0f;
	}
	else{
		return 0.0f;
	}
}

__device__ datatype __BC_FunctionVort(const int x, const int y, const int maxLayer, const int maxGlobIndex){

	//This needs to be fixed for the corners and stuff. FIX!
	if (y == maxGlobIndex)
	{
		return 1.0f;
	}
	else if(y == 0){
		return -1.0f;
	}
	else{
		return 0.0f;
	}
}

__device__ int __findNodeGlobIdx(const Node* u, const int x, const int y, const int len){

	for (int i = 0; i < len; ++i)
	{
		if (u[i].x_index_global == x && u[i].y_index_global)
		{
			return i;
		}
	}
	return -1;

}

__device__ datatype __findNodeGlobStream(const Node* u, const int x, const int y, const int len){
	for (int i = 0; i < len; ++i)
	{
		if (u[i].x_index_global == x && u[i].y_index_global)
		{
			return u[i].stream;
		}
	}
	return -1;
}

__device__ int __pow_2(const datatype x){
	return x*x;
}

__device__ datatype __findValOfClosestPoint(const Node* u, const int x, const int y, const int len){
	int idx;
	int dist = 10000000;
	datatype tmp;
	for (int i = 0; i < len; ++i)
	{
		tmp = __pow_2(u[i].x_index_global-x)+__pow_2(u[i].x_index_global-x);
		if( tmp < dist){
			dist = tmp;
			idx = i;
		}
	}
	return u[idx].stream;
}

__device__ datatype __interpolateGhostPoint(const Node* u, const Node* u_coarse, const int len_coarse, const int x, const int y, const int layerNr, const int maxLayer, const int maxGlobIndex, const datatype h){

	//So this is stupid. But it should be quick. I hope.
	//Instead of interpolating the ghost points it just takes the value of the closest point
	//from the coarser grid.
	//FIX!

	return __findValOfClosestPoint(u_coarse, x, y, len_coarse);


	/*
	if ( x%2 == 0)
	{
		return ( u[__findNodeGlobIdx()].stream + u[__findNodeGlobIdx()].stream ) / 2.0f;
	}
	else if( y%2 == 0){

	}
	else{

	}
	*/
}

__global__ void updateBFromInterpolation(Node* b, const Node* u, const int len, const Node* u_coarse, const int len_coarse, const int layerNr, const int maxLayer, const int maxGlobIndex, const datatype h){
	const int id = threadIdx.x + blockIdx.x*blockDim.x;

	if (id < len)
	{
		const int 	   n 		= 	maxLayer-layerNr;
		const int 	   step 	=	1<<(n);
		const datatype tmp_val 	= 	dev_invPow(2.0f, n);
		const datatype gamma 	= 	2.0f/( (h*h) * (tmp_val+tmp_val*tmp_val) );
		

		if(u[id].nodeRight == NULL){
			if (u[id].x_index_global == maxGlobIndex)
			{
				b[id].stream += -gamma*__BC_Function(u[id].x_index_global, u[id].y_index_global, maxLayer, maxGlobIndex);
			}
			else{
				datatype tmp = __interpolateGhostPoint(u, u_coarse, len_coarse, u[id].x_index_global+step, u[id].y_index_global, layerNr, maxLayer, maxGlobIndex, h);
				b[id].stream += -(tmp)/(h*h);
			}
		}
		if (u[id].nodeLeft == NULL){
			if(u[id].x_index_global == 0){
				b[id].stream += -gamma*__BC_Function(u[id].x_index_global, u[id].y_index_global, maxLayer, maxGlobIndex);
			}
			else{
				datatype tmp = __interpolateGhostPoint(u, u_coarse, len_coarse, u[id].x_index_global-step, u[id].y_index_global, layerNr, maxLayer, maxGlobIndex, h);
				b[id].stream += -(tmp)/(h*h);
			}
			
		}
		if(u[id].nodeAbove == NULL){
			if (u[id].y_index_global == maxGlobIndex)
			{
				b[id].stream += -gamma*__BC_Function(u[id].x_index_global, u[id].y_index_global, maxLayer, maxGlobIndex);
			}
			else{
				datatype tmp = __interpolateGhostPoint(u, u_coarse, len_coarse, u[id].x_index_global, u[id].y_index_global+step, layerNr, maxLayer, maxGlobIndex, h);
				b[id].stream += -(tmp)/(h*h);
			}
		}
		if (u[id].nodeBelow == NULL){
			if (u[id].y_index_global == 0)
			{
				b[id].stream += -gamma*__BC_Function(u[id].x_index_global, u[id].y_index_global, maxLayer, maxGlobIndex);
			}
			else{
				datatype tmp = __interpolateGhostPoint(u, u_coarse, len_coarse, u[id].x_index_global, u[id].y_index_global-step, layerNr, maxLayer, maxGlobIndex, h);
				b[id].stream += -(tmp)/(h*h);
			}
		}
	}
}

__global__ void dev_updateBFromBoundary(Node* b, const Node* u, const int len, const int layerNr, const int maxLayer, const int maxGlobIndex, const datatype h){

	const int id = threadIdx.x + blockIdx.x*blockDim.x;

	if (id < len)
	{
		if (is1StepFromBoundary(u, id, maxGlobIndex))
		{
			const int n = maxLayer-layerNr;
			const datatype tmp = dev_invPow(2.0f, n);
			const datatype gamma = 2.0f/( (h*h) * (tmp+tmp*tmp) );

			if (isInCorner(u, id, maxGlobIndex))
			{
				b[id].stream += -2.0f*gamma*__BC_Function(u[id].x_index_global, u[id].y_index_global, maxLayer, maxGlobIndex);
			}
			else{
				b[id].stream += -gamma*__BC_Function(u[id].x_index_global, u[id].y_index_global, maxLayer, maxGlobIndex);
			}
		}
	}
}

void setupBoundaryOnCoarsestGrid(Node* bc, const int len , func_def BC_Function){

	int counter = 0;

	int x = -1;
	int y;
	for (int y = -1; y < len+1; ++y)
	{
		bc[counter].stream = BC_Function(x,y);
		counter++;
	}

	x = len;
	for (int y = -1; y < len+1; ++y)
	{
		bc[counter].stream = BC_Function(x,y);
		counter++;
	}

	y = -1;
	for (int x = -1; x < len+1; ++x)
	{
		bc[counter].stream = BC_Function(x,y);
		counter++;
	}

	y = len;
	for (int x = -1; x < len+1; ++x)
	{
		bc[counter].stream = BC_Function(x,y);
		counter++;
	}
}

#ifdef OLDCODE
__device__ datatype getLaplacianStream(const Node * u, const int index1, const int index2, const datatype h){
	const int dist = abs_dev(u[index1].x_index - u[index2].x_index) +abs_dev(u[index1].y_index - u[index2].y_index);

	if( dist > 1){
		 	return 0.0f;
		 }
		 else{
		 	if( dist == 0 ){
		 		return -4.0f/(h*h);
		 	}
			else{
			 	//if(isInBoundary(u[index2].x_index, u[index2].y_index)){
			 	//	return getFromBoundary(u[index2].x_index, u[index2].y_index)/h;
			 	//}
			 	//else{
			 		return 1.0f/(h*h);
			 	//}
			}
		}
}
#endif

__global__ void calculateErrorLaplacian(const Node* b, const Node* u, Node* d, const int len, const datatype h, const int maxGlobIndex, const int layerNr, const int maxLayerNr){
	int id = threadIdx.x + blockIdx.x*blockDim.x;

	if(id < len){
		datatype sum = 0.0f;
		for(int j=0; j<len; j++){
			sum += getLaplacianStream(u, j, id, h, maxGlobIndex, layerNr, maxLayerNr)*u[j].stream;
		}
		d[id].stream = b[id].stream - sum;
	}

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

__global__ void calculateRHS(const Node* u, const Node* d, Node* b, const int len, const datatype h, const int maxGlobIndex, const int layerNr, const int maxLayerNr){
	const int id = threadIdx.x + blockIdx.x*blockDim.x;
	if( id < len){

		datatype sum = 0.0f;
		for (int j = 0; j < len; ++j)
		{
			sum += getLaplacianStream(u, j, id, h , maxGlobIndex, layerNr, maxLayerNr )*u[j].stream;
		}
		b[id].stream = d[id].stream + sum;
	}
}

__global__ void copy(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		to[i].stream = from[i].stream;
	}
}

__global__ void copyVort(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		to[i].vort = from[i].vort;
	}
}

__global__ void subtract_gpu(Node* a, const Node* b, const Node* c, const int len){
	//a = b - c
	const int i = threadIdx.x + blockIdx.x*blockDim.x;
	if( i < len ){
		a[i].stream = b[i].stream - c[i].stream;
	}
}

__global__ void add_gpu(Node* to, const Node* from, const int len){
	const int i = threadIdx.x + blockIdx.x*blockDim.x;

	if( i < len ){
		to[i].stream += from[i].stream;
	}
}

__device__ datatype findNodeIdx(const Node* arr, const int x, const int y, const int n){
	for (int i = 0; i < n; ++i)
	{
		if (arr[i].x_index == x && arr[i].y_index == y)
		{
			return i;
		}
	}
	return -1;
}

__device__ datatype findNodeVal(const Node* arr, const int x, const int y, const int n){
	for (int i = 0; i < n; ++i)
	{
		if (arr[i].x_index == x && arr[i].y_index == y)
		{
			return arr[i].stream;
		}
	}
	printf("			SOMETHING HAS FUCKED UP! (stream edition)");
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

#ifdef OLDCODE
__global__ void interpolate_old( Node* u_fine, Node* u_coarse, const int len_fine, const int len_coarse){
	//This isn't done yet. FIX!
	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < len_fine){


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

		syncthreads();

		printf("threadIdx: %d, isCoarse: %d, x: %d, y: %d, stream: %f \n", i, isOnCoarserGrid, u_fine[i].x_index, u_fine[i].y_index, u_fine[i].stream);

		if( isOnCoarserGrid == true ){
			Node *middleNode;
			middleNode = &u_fine[i];

			datatype tmp = -7.0f;




			//Interpolate right
			if(middleNode->nodeRight != NULL){
				printf("node right: %p \n", middleNode->nodeRight);
				printf("node right right: %p \n", middleNode->nodeRight->nodeRight);
				//return;
				if(middleNode->nodeRight->nodeRight != NULL){
					//return;
					//middleNode->nodeRight->stream =  6;//middleNode->stream;// + middleNode->nodeRight->nodeRight->stream ) / 2.0f;
					tmp = 0.2;//middleNode->nodeRight->stream;
				}
			}

			printf("tmp = %f \n", tmp);

			//return;

			//Interpolate left
			if(middleNode->nodeLeft != NULL){
				if(middleNode->nodeLeft->nodeLeft != NULL){
					middleNode->nodeLeft->stream =  ( middleNode->stream + middleNode->nodeLeft->nodeLeft->stream ) / 2.0f;
				}
			}

			//Interpolate up
			if(middleNode->nodeAbove != NULL){
				if(middleNode->nodeAbove->nodeAbove != NULL){
					middleNode->nodeAbove->stream =  ( middleNode->stream + middleNode->nodeAbove->nodeAbove->stream ) / 2.0f;
				}
			}

			//Interpolate down
			if(middleNode->nodeBelow != NULL){
				if(middleNode->nodeBelow->nodeBelow != NULL){
					middleNode->nodeBelow->stream =  ( middleNode->stream + middleNode->nodeBelow->nodeBelow->stream ) / 2.0f;
				}
			}

			//And at last, the middle point.
			if(middleNode->nodeBelow != NULL){
				if(middleNode->nodeBelow->nodeRight != NULL){
					middleNode->nodeBelow->nodeRight->stream =  ( middleNode->nodeBelow->stream + middleNode->nodeBelow->nodeRight->nodeRight->stream ) / 2.0f;
				}
			}

		}

		
	}
}
#endif

__global__ void findNeighbours(Node* array, const int len ){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
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

__global__ void jacobiSmootherLaplacianStream(Node *from, Node * to, Node * b, const datatype h, const int len, const int maxGlobIndex, const int layerNr, const int maxLayerNr){
	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < len){

		datatype sum = 0.0f;
		for (int j = 0; j < len; ++j)
		{
			if(j != i){
				sum += getLaplacianStream(from, i, j, h, maxGlobIndex, layerNr, maxLayerNr) * from[j].stream;
			}
		}
		to[i].stream = (b[i].stream - sum) / getLaplacianStream(from, i,i, h , maxGlobIndex, layerNr, maxLayerNr);
	}
}

void ERRORCHECK(){
	cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("Error (!): %s\n", cudaGetErrorString(err));
    exit(-1);
  }
}

__global__ void setVectorsToZero(Node * arr, const int len){
	int i = threadIdx.x + blockIdx.x*blockDim.x;

	if(i < len){
		arr[i].stream = 0.0f;
		//grid->d[i].stream = 0.0f;
		//grid->w[i].stream = 0.0f;
	}
}

__global__ void dev_updateBFromVort(Node* b, Node* u, const int len){
	const int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id < len)
	{
		b[id].stream += -u[id].vort;
	}
}


#define THREADS_PER_BLOCK 17

void multigrid_gpu(int k, AdaptiveGrid* grid, int pre, int sol, int post, const int maxGlobIndex, const int maxLayerNr){

	const int N = grid->len;
	const int N_coarse = grid->coarserGrid->len;

	//Define the size of the current grid in a CUDA format
	dim3 block_size_1d(N);
	dim3 grid_size_1d((N+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	//And the same for the coarser grid.
	dim3 block_size_1d_coarse(N_coarse);
	dim3 grid_size_1d_coarse((N_coarse+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	//Pre-smoothing
	for (int i = 0; i < k*pre/2 + 1; ++i)
	{
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->u, grid->d, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->d, grid->u, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
	}

	//std::cout<<"  1  "<<std::endl;
	ERRORCHECK();

	//Calculate error
	//d = f-L*u;
	calculateErrorLaplacian<<<grid_size_1d, block_size_1d>>>(grid->b, grid->u, grid->d, grid->len, grid->h, maxGlobIndex, grid->layerNr, maxLayerNr);

	//std::cout<<"  2  "<<std::endl;
	ERRORCHECK();

	//Restrict d
	restrictMat<<<grid_size_1d, block_size_1d>>>(grid->coarserGrid->d, grid->coarserGrid->len, grid->d , grid->len ); 
	
	//Restrict u
	restrictMat<<<grid_size_1d, block_size_1d>>>(grid->coarserGrid->u, grid->coarserGrid->len, grid->u , grid->len ); 

	//Calculate the right hand side b for the next layer.
	//calculateRHS<<<1, block_size_1d_coarse>>>(grid->coarserGrid->u, grid->coarserGrid->d, grid->coarserGrid->b, grid->coarserGrid->len, grid->coarserGrid->h, maxGlobIndex, grid->coarserGrid->layerNr, maxLayerNr);

	//Copy u to w.
	copy<<<grid_size_1d, block_size_1d>>>(grid->coarserGrid->w, grid->coarserGrid->u, grid->coarserGrid->len);


	//std::cout<<"  3  "<<std::endl;
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

	//std::cout<<"  4  "<<std::endl;
	ERRORCHECK();

	subtract_gpu<<<grid_size_1d_coarse, block_size_1d_coarse>>>(grid->coarserGrid->d, grid->coarserGrid->u, grid->coarserGrid->w, grid->coarserGrid->len);

	//std::cout<<"  5  "<<std::endl;
	ERRORCHECK();
	interpolate<<<grid_size_1d, block_size_1d>>>( grid->d, grid->coarserGrid->d, grid->len, grid->coarserGrid->len, grid->layerNr, maxLayerNr, maxGlobIndex, grid->h, grid->b);

	//std::cout<<"  6  "<<std::endl;
	ERRORCHECK();
	//std::cout<<"  7  "<<std::endl;

	add_gpu<<<grid_size_1d, block_size_1d>>>(grid->w, grid->d, grid->len);

	//Copy u from w
	copy<<<grid_size_1d, block_size_1d>>>(grid->coarserGrid->u, grid->coarserGrid->w, grid->coarserGrid->len);

	//update b
	updateBNew<<<grid_size_1d, block_size_1d>>>(grid->u, grid->b, grid->coarserGrid->u, grid->len, grid->coarserGrid->len_coarse, maxGlobIndex);
	//setVectorsToZero<<<grid_size_1d, block_size_1d>>>(grid->b, grid->len);
	//dev_updateBFromBoundary<<<grid_size_1d, block_size_1d>>>(grid->b, grid->u, grid->len, grid->layerNr, maxLayerNr, maxGlobIndex, grid->h);
	//dev_updateBFromVort<<<grid_size_1d, block_size_1d>>>(grid->b, grid->u, grid->len);
	//updateBFromInterpolation<<<grid_size_1d, block_size_1d>>>(grid->b, grid->u, grid->len, grid->coarserGrid->u, grid->coarserGrid->len, grid->layerNr, maxLayerNr, maxGlobIndex, grid->h);

	//Post-smoothing
	for (int i = 0; i < k*post/2 + 1; ++i)
	{
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->w, grid->u, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
		jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->u, grid->w, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
	}
	jacobiSmootherLaplacianStream<<<grid_size_1d, block_size_1d>>>( grid->w, grid->u, grid->b, grid->h, grid->len, maxGlobIndex, grid->layerNr, maxLayerNr);
}

__global__ void printAll(Node * arr, const int len){
	int i = threadIdx.x + blockDim.x*blockIdx.x;

	if(i<len){
		printf("threadIdx: %d, x: %d, y: %d, stream: %f , ptr: %p \n", i, arr[i].x_index, arr[i].y_index, arr[i].stream, &arr[i]);
	}
}

__global__ void gpuPrint(Node* arr, const int len){
	int i = threadIdx.x + blockDim.x*blockIdx.x;

	if(i<len){
		printf("threadIdx: %d, x: %d, y: %d, stream: %f \n", i, arr[i].x_index, arr[i].y_index, arr[i].stream);
	}
}

void move2gpu(AdaptiveGrid * grid){

	Node *u_dev, *b_dev, *w_dev, *d_dev;

	dim3 block_size_1d(grid->len);

	/*Node* tmp, *tmp2;
	std::cout<<"before \n";
	std::cout<<"node size: "<<sizeof(Node)<<"\n";
	cudaMalloc(&tmp,0);
	std::cout<<"middle \n";

	//cudaMalloc(&tmp2,1);
	//assert(0);
	*/
	std::cout<<"grid len: "<<grid->len<<" \nlen*sizeofNode: "<<grid->len*sizeof(Node)<<std::endl;

	cudaMalloc( (void**) &u_dev, grid->len*sizeof(Node));
	std::cout<<" test 1 \n";
	cudaMalloc( (void**) &b_dev, grid->len*sizeof(Node));
	std::cout<<" test 2 \n";
	cudaMalloc( (void**) &w_dev, grid->len*sizeof(Node));
	std::cout<<" test 3 \n";
	cudaMalloc( (void**) &d_dev, grid->len*sizeof(Node));

	std::cout<<"1\n";

	cudaMemcpy( u_dev, grid->u, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	cudaMemcpy( b_dev, grid->b, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	cudaMemcpy( w_dev, grid->w, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);
	cudaMemcpy( d_dev, grid->d, sizeof(Node)*grid->len, cudaMemcpyHostToDevice);

	std::cout<<"2\n";

	free(grid->u);
	free(grid->b);
	free(grid->w);
	free(grid->d);

	std::cout<<"3\n";

	grid->u = u_dev;
	grid->w = w_dev;
	grid->d = d_dev;
	grid->b = b_dev;

	std::cout<<"findNeighbours\n";

	findNeighbours<<<1, block_size_1d>>>( grid->u, grid->len );
	findNeighbours<<<1, block_size_1d>>>( grid->b, grid->len );
	findNeighbours<<<1, block_size_1d>>>( grid->w, grid->len );
	findNeighbours<<<1, block_size_1d>>>( grid->d, grid->len );

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

	std::cout<<"Moved data from the gpu."<<std::endl;
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

		Node tmp = Node(h*(x_min+diff)/2.0f, h*(y_min+diff)/2.0f, *x_loc, *y_loc,(x_min+diff/2), (y_min+diff/2), 0.0f, 0.0f, 0.0f );
		gridList[LAYERS - layer + 1]->vec.push_back(tmp);
		if (layer != 2)
		{
			recursiveGridFill(u, layer-1, gridList, firstPtr, x_min, x_min+diff, y_min, y_min+diff, h);
		}
	}

	if (isLeaf(u, firstPtr[LAYERS - layer + 1], firstPtr[LAYERS], x_min+diff, x_min+2*diff, y_min, y_min+diff) == false)
	{
		gridList[LAYERS - layer + 1]->global2local((x_min+diff+diff/2), (y_min+diff/2),x_loc, y_loc);

		Node tmp = Node(h*(x_min+diff)/2.0f, h*(y_min+diff)/2.0f, *x_loc, *y_loc,(x_min+diff+diff/2), (y_min+diff/2), 0.0f, 0.0f, 0.0f );
		gridList[LAYERS - layer + 1]->vec.push_back(tmp);
		if (layer != 2)
		{
			recursiveGridFill(u, layer-1, gridList, firstPtr,x_min+diff, x_min+2*diff, y_min, y_min+diff, h);
		}
	}

	if (isLeaf(u, firstPtr[LAYERS - layer + 1], firstPtr[LAYERS], x_min+diff, x_min+2*diff, y_min+diff, y_min+2*diff) == false)
	{
		gridList[LAYERS - layer + 1]->global2local((x_min+diff+diff/2), (y_min+diff +diff/2),x_loc, y_loc);

		Node tmp = Node(h*(x_min+diff)/2.0f, h*(y_min+diff)/2.0f, *x_loc, *y_loc,(x_min+diff+diff/2), (y_min+diff +diff/2), 0.0f, 0.0f, 0.0f );
		gridList[LAYERS - layer + 1]->vec.push_back(tmp);
		if (layer != 2)
		{
			recursiveGridFill(u, layer-1, gridList, firstPtr,  x_min+diff, x_min+2*diff, y_min+diff, y_min+2*diff, h);
		}
	}

	if (isLeaf(u, firstPtr[LAYERS - layer + 1], firstPtr[LAYERS], x_min, x_min+diff, y_min+diff, y_min+2*diff) == false)
	{
		gridList[LAYERS - layer + 1]->global2local((x_min+diff/2), (y_min+diff+diff/2),x_loc, y_loc);

		Node tmp = Node(h*(x_min+diff)/2.0f, h*(y_min+diff)/2.0f, *x_loc, *y_loc,(x_min+diff/2), (y_min+diff+diff/2), 0.0f, 0.0f, 0.0f );
		gridList[LAYERS - layer + 1]->vec.push_back(tmp);
		if (layer != 2)
		{
			recursiveGridFill(u, layer-1, gridList, firstPtr, x_min, x_min+diff, y_min+diff, y_min+2*diff, h);
		}
	}
	free(y_loc);
	free(x_loc);
}

__device__ bool isInterior(Node* arr, const int i){
	return (arr[i].nodeRight != NULL && arr[i].nodeLeft != NULL && arr[i].nodeBelow != NULL && arr[i].nodeAbove != NULL);
}

__global__ void vortInterior(Node* to, Node* from, const int len, const datatype h){
	const int id = threadIdx.x + blockDim.x*blockIdx.x;

	if (id < len)
	{
		if (isInterior(from, id))
		{
			//This is slow and dumb and shared memory won't help until the memory gets a proper structure. OPT!
			const datatype vort_dx2 = (findNodeValVort(from, from[id].x_index+1, from[id].y_index, len) + findNodeValVort(from, from[id].x_index-1, from[id].y_index, len) - 2.0f* from[id].vort)/(h*h);
			const datatype vort_dy2 = (findNodeValVort(from, from[id].x_index, from[id].y_index+1, len) + findNodeValVort(from, from[id].x_index, from[id].y_index-1, len) -2.0f* from[id].vort)/(h*h);
			//printf("Has anything fucked up? ");
			const datatype stream_dx = (findNodeVal(from, from[id].x_index+1, from[id].y_index, len) - findNodeVal(from, from[id].x_index-1, from[id].y_index, len))/(2.0f*h);
			//printf("Has anything fucked up? ");
			const datatype stream_dy = (findNodeVal(from, from[id].x_index, from[id].y_index+1, len) - findNodeVal(from, from[id].x_index, from[id].y_index-1, len))/(2.0f*h);

			const datatype vort_dx = (findNodeValVort(from, from[id].x_index+1, from[id].y_index, len) - findNodeValVort(from, from[id].x_index-1, from[id].y_index, len))/(2.0f*h);
			const datatype vort_dy = (findNodeValVort(from, from[id].x_index, from[id].y_index+1, len) - findNodeValVort(from, from[id].x_index, from[id].y_index-1, len))/(2.0f*h);

			to[id].vort = (vort_dy2+vort_dx2)/Re-stream_dy*vort_dx+stream_dx*vort_dy;
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
					to[id].vort =  2.0f/(h*h)*(from[id].stream - findNodeVal(from, from[id].x_index, from[id].y_index-1, len))+1.0f;;
				}
				else if(from[id].y_index_global == 0){
					to[id].vort =  2.0f/(h*h)*(from[id].stream - findNodeVal(from, from[id].x_index, from[id].y_index+1, len));
				}
				else if(from[id].x_index_global == maxGlobIndex){
					to[id].vort =  2.0f/(h*h)*(from[id].stream - findNodeVal(from, from[id].x_index-1, from[id].y_index, len));
				}
				else if(from[id].x_index_global == 0){
					to[id].vort =  2.0f/(h*h)*(from[id].stream - findNodeVal(from, from[id].x_index+1, from[id].y_index, len));
				}
				else{
					printf("ooooooooooooooookej. Vort exterior is a little bit funky.");
				}
			}
			else{
				//Update according to a closest neigbour.
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

void calcVortFromStream(AdaptiveGrid* grid){
	const int N = grid->len;

	//Define the size of the current grid in a CUDA format
	dim3 block_size_1d(N);
	dim3 grid_size_1d((N+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);

	//Updates the vorticity and puts the new values in w.
	vortInterior<<<grid_size_1d, block_size_1d>>>(grid->w, grid->u, grid->len, grid->h);
	vortExterior<<<grid_size_1d, block_size_1d>>>(grid->w, grid->u, grid->len, grid->h, 1<<LAYERS);

	//
	copyVort<<<grid_size_1d, block_size_1d>>>(grid->u, grid->w, grid->len);
}

void adaptive_multigrid(Node* array, int* origoArray, int countTrue, int layers){

	std::cout<<"hej\n";

	int * pointsInLayer;
	pointsInLayer = (int*) calloc(LAYERS+1, sizeof(int));

	int bigCount = 1;

	int tmpLayer = layers;
	int counter = 0;



	for(int i = 0; i<tmpLayer*2; i++){

		std::cout<<origoArray[i]<<std::endl;

	}

/*
	//global2local
	int idx = tmpLayer*2 -2;
	int idy = idx +1;

	for(int j = tmpLayer; j>0; j--){

		for(int i = 0; i<countTrue; i++){

			if(array[i].layer == j){

				array[i].x_index = (array[i].x_index_global - origoArray[idx])/ (1<<(j-1));
				array[i].y_index = (array[i].y_index_global - origoArray[idy])/ (1<<(j-1));

				//std::cout<<"origo x: "<<origoArray[idx]<<" globalx: "<<array[i].x_index_global<<" Local: "<<array[i].x_index<<std::endl<<"origo y: "<<origoArray[idy]<<" globaly: "<<array[i].y_index_global<<" Local: "<<array[i].y_index<<std::endl;
			}
		}

		idx = idx -2;
		idy = idx +1;
	}	

*/
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


	for (int i = 0; i < 4; ++i)
	{
		//std::cout<<"points: "<<pointsInLayer[i]<<std::endl;
	}

	AdaptiveGrid** gridList;
	gridList = (AdaptiveGrid**) malloc(layers*sizeof(AdaptiveGrid*));

	std::cout<<"adaptive_multigrid\n";

	//assert(0);
	int pointCounter = 0;
	for (int i = 1; i <= layers; ++i)
	{
		AdaptiveGrid* grid = new AdaptiveGrid(i,LAYERS, origoArray[2*(LAYERS-i)],origoArray[2*(LAYERS-i)+1], NULL, 0,1.0f/(1<<LAYERS));
		pointCounter += pointsInLayer[i-1];
		gridList[i-1] = grid;
	}

	for (int i = 0; i < 5; ++i)
	{
		gridList[0]->vec.push_back(array[i]);
	}

	assert(row == colum);
	recursiveGridFill(array, LAYERS, gridList, pointsInLayer, 0, row-1, 0, row-1, 1.0f/(1<<LAYERS));

	for (int i = 0; i < LAYERS; ++i)
	{
		gridList[i]->setupFromVector();
	}






















	

	std::cout<<"move2gpu\n";
	for (int i = 0; i < LAYERS; ++i)
	{
		move2gpu(gridList[i]);
		if(i>0)
			gridList[i]->coarserGrid = gridList[i-1];
		if(i<LAYERS-1)
			gridList[i]->finerGrid = gridList[i+1];
	}


	for (int i = 0; i < 80; ++i)
	{
		//std::cout<<"muuuuuultigrid!\n";
		//setVectorsToZero<<<1, grid1.len>>>(grid1.b, grid1.len);
		dev_updateBFromBoundary<<<1, gridList[0]->len>>>(gridList[0]->b, gridList[0]->u, gridList[0]->len, gridList[0]->layerNr, 4, row-1, gridList[0]->h);
		multigrid_gpu(layers-1, gridList[layers-1], 20, 40, 20, row-1, layers);

	}

	for (int i = 0; i < LAYERS; ++i)
	{
		calcVortFromStream(gridList[i]);	
	}

	for (int i = 0; i < LAYERS; ++i)
	{
		move2host(gridList[i]);
	}

	counter = 0;

	for (int j = 0; j < countTrue; ++j)
	{
		//std::cout<<"j "<<j<<"\n";
		for(int lay=0; lay<LAYERS; lay++){
			Node* tmp;

			tmp = gridList[lay]->findGlobNodeGeneral(gridList[lay]->u, 
				array[j].x_index_global, 
				array[j].y_index_global, 
				gridList[lay]->len);

			//std::cout<<"lay "<<lay<<"\n";

			if (tmp != NULL){
				array[j].vort = tmp->vort;
				break;
			}
			assert(lay != LAYERS-1);
		}
	}

	std::cout<<"free stuff in adaptive_multigrid_cuda.cu"<<std::endl;

	for (int i = 0; i < LAYERS; ++i)
	{
		delete gridList[i];
	}

	free(gridList);
	free(pointsInLayer);

}
