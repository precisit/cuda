#include <assert.h>
#include "node.cpp"
class AdaptiveGrid{
public:

	Node* u;
    int len;
    Node* b;

    Node* boundaryVals;
    int* boundaryIndex;
    int lenOfBoundary;

    Node** pointsChosenByWavelet; //An array of Node pointers
    int numOfPointsChosen;

    int layerNr;
    int numOfLayers;
    int origo_x;
    int origo_y;

    AdaptiveGrid* finerGrid = NULL;
    AdaptiveGrid* coarserGrid = NULL;
    datatype h;

	bool isInBoundary(const int x, const int y){
		//std::cout<<"hello there!: "<<lenOfBoundary<<std::endl;
		for (int i = 0; i < lenOfBoundary; ++i)
		{
			if(boundaryIndex[2*i]==x && boundaryIndex[2*i+1]==y){
				//*val = boundaryVals[i];
				return true;
			}
		}
		return false;
	}

	void setB(const datatype val, const int i){
		assert(i<len && i>=0);
		b[i].stream = val;
	}

	datatype getB(const int i){
		assert(i<len && i>=0);
		return b[i].stream;
	}

	Node getFromBoundary(const int x, const int y){
		std::cout<<"hi!"<<std::endl;
		for (int i = 0; i < lenOfBoundary; ++i)
		{
			if(boundaryIndex[2*i] == x && boundaryIndex[2*i+1] == y){
				return boundaryVals[i];
			}
		}
		std::cout<<x<<" "<<y<<std::endl;
		assert(0);
		return boundaryVals[0];
	}

	datatype getLaplacianStream(const int index1, const int index2){

		assert(index1 >= 0 && index1 < len);
		assert(index2 >= 0);
		assert(index2 < len);

		int dist = abs(u[index1].x_index - u[index2].x_index) +abs(u[index1].y_index - u[index2].y_index);

		if( dist > 1){
		 	return 0.0f;
		 }
		 else{
		 	if( dist == 0 ){
		 		return -4.0f/(h*h);
		 	}
			else{
			 	if(isInBoundary(u[index2].x_index, u[index2].y_index)){
			 		return getFromBoundary(u[index2].x_index, u[index2].y_index).stream/h;
			 	}
			 	else{
			 		return 1.0f/(h*h);
			 	}
			}
		}
	}

	void resetBoundaryLength(){
		lenOfBoundary = 0;
	}

	bool isInGrid(const int x, const int y){
		for (int i = 0; i < len; ++i)
		{
			if(u[i].x_index==x && u[i].y_index==y){
				return true;
			}
		}
		//std::cout<<x<<" ; "<<y<<std::endl;
		return false;
	}

	void setBoundaryLength(){
		//std::cout<<"                    "<<std::endl;
		int counter = 0;
		int x, y;
		for (int i = 0; i < len; ++i)
		{
			x = u[i].x_index + 1;
			y = u[i].y_index;

			if(isInGrid(x,y)==false){
				++counter;
			}
			x = u[i].x_index - 1;
			y = u[i].y_index;

			if(isInGrid(x,y)==false){
				++counter;
			}
			x = u[i].x_index;
			y = u[i].y_index + 1;

			if(isInGrid(x,y)==false){
				++counter;
			}
			x = u[i].x_index;
			y = u[i].y_index - 1;

			if(isInGrid(x,y)==false){
				++counter;
			}
		}
		lenOfBoundary = counter;
	}

	//void setB(const int i, const datatype val){
	//	assert(i>=0 && i<len);
//
//		b[i].stream = val;
//	}


	void setU(const datatype val, const int i){
		assert(i<len && i>=0);
		u[i].stream = val;
	}

	datatype getU(const int i){
		assert(i<len && i>=0);
		return u[i].stream;
	}

	void jacobiSmootherLaplacianStream(){
		/*
			Smoothes the error to D u = b.
			u = (p, v_x, v_y)^T
		*/

			datatype tmpSum;
			for(int i=0; i< len; i++){
				tmpSum = 0.0f;
				for(int j=0; j< len ; j++){
					if(j != i){

						tmpSum += this->getLaplacianStream(i,j)*(this->getU(j));

					}
				}
				setU( (this->getB(i)-tmpSum) / this->getLaplacianStream(i,i), i );
			}
		}

	void restrictU(AdaptiveGrid coarse){
		for (int i = 0; i < coarse.len; ++i)
		{
			for (int j = 0; j < this->len; ++j)
			{
				if( this->u[j].x_index_global == coarse.u[i].x_index_global 
					&& this->u[j].y_index_global == coarse.u[i].y_index_global )
				{
					coarse.u[i].stream = this->u[j].stream;
				}
			}
		}
	}


	Node * findNode(const int ind_x, const int ind_y){

		for (int i = 0; i < len; ++i)
		{
			if(u[i].x_index == ind_x && u[i].y_index == ind_y){
				return &u[i];
			}
		}
		//std::cout<<ind_x<<"; "<<ind_y<<std::endl;
		//assert(0);
		return NULL;

	}



	void interpolateU(AdaptiveGrid fine){
		/*for (int i = 0; i < fine.len; ++i)
		{
			for (int j = 0; j < this->len; ++j)
			{
				if( this->u[j].x_index_global == fine.u[i].x_index_global 
					&& this->u[j].y_index_global == fine.u[i].y_index_global )
				{
					fine.u[i].stream = this->u[j].stream;
				}
			}
		}*/

		for (int i = 0; i < numOfPointsChosen; ++i)
		{

			//std::cout<<"i: "<<i<<std::endl;
			Node *middleNode = findNode(pointsChosenByWavelet[i]->x_index, pointsChosenByWavelet[i]->y_index);
			Node *tmpNodePtr;
			//Interpolate to the left
			tmpNodePtr = middleNode->nodeLeft;

			//assert(tmpNodePtr != NULL);
			//assert(tmpNodePtr->nodeBelow != NULL);
			//assert(tmpNodePtr->nodeAbove != NULL);

			tmpNodePtr->stream = (tmpNodePtr->nodeAbove->stream + tmpNodePtr->nodeBelow->stream )/2.0f;

			//to the right
			tmpNodePtr = middleNode->nodeRight;
			tmpNodePtr->stream = (tmpNodePtr->nodeAbove->stream + tmpNodePtr->nodeBelow->stream )/2.0f;

			//up
			tmpNodePtr = middleNode->nodeAbove;
			tmpNodePtr->stream = (tmpNodePtr->nodeLeft->stream + tmpNodePtr->nodeRight->stream )/2.0f;

			//down
			tmpNodePtr = middleNode->nodeBelow;
			tmpNodePtr->stream = (tmpNodePtr->nodeLeft->stream + tmpNodePtr->nodeRight->stream )/2.0f;

			//and don't forget the point itself.
			middleNode->stream = ( middleNode->nodeLeft->stream + 
				middleNode->nodeRight->stream )/2.0f;
		}
	}



	void setNeighbours(){

		for (int i = 0; i < numOfPointsChosen; ++i)
		{
			//Up
			(pointsChosenByWavelet[i]->nodeAbove) = findNode(pointsChosenByWavelet[i]->x_index, pointsChosenByWavelet[i]->y_index + 1);

			//down
			(pointsChosenByWavelet[i]->nodeBelow) = findNode(pointsChosenByWavelet[i]->x_index, pointsChosenByWavelet[i]->y_index - 1);
			
			//Left
			(pointsChosenByWavelet[i]->nodeLeft) = findNode(pointsChosenByWavelet[i]->x_index-1, pointsChosenByWavelet[i]->y_index);

			//Right
			(pointsChosenByWavelet[i]->nodeRight) = findNode(pointsChosenByWavelet[i]->x_index+1, pointsChosenByWavelet[i]->y_index);
			

		}
	}

	void local2global(const int* x_loc, const int* y_loc, int* x_glo, int* y_glo ){

		assert( this->layerNr > 0 );
		assert( this->layerNr <= this->numOfLayers );

		*y_glo = *y_loc*(1<<(numOfLayers-layerNr))+this->origo_y;
		*x_glo = *x_loc*(1<<(numOfLayers-layerNr))+this->origo_x;

	}


	void setupGrid(Node* savedNodes, const int numberOfPoints){
		//free(pointsChosenByWavelet);
		pointsChosenByWavelet = (Node**) malloc(numberOfPoints*sizeof(Node*));
		this->numOfPointsChosen = numberOfPoints;

		//free(u);
		//9 is the maximum. It's probably waaaay too large. OPT!
		u = (Node*) malloc(9*numberOfPoints*sizeof(Node));
		len = 9*numberOfPoints;

		int counter = numberOfPoints;

		//Take the points given by 
		for (int i = 0; i < numberOfPoints; ++i)
		{
			pointsChosenByWavelet[i] = &savedNodes[i];
			u[i] = savedNodes[i];
			
			//Left
			if(isInGrid(savedNodes[i].x_index-1, savedNodes[i].y_index) == false){
				//assert(i != 0);
				Node tmpNode;
				tmpNode = u[i];
				u[i].nodeLeft = &tmpNode;
				u[counter] = tmpNode;
				u[counter].x_index--;
				
				//Make sure to do something about the global index here! FIX!
				counter++;
			}
			//Right
			if(isInGrid(savedNodes[i].x_index+1, savedNodes[i].y_index) == false){
				Node tmpNode;
				tmpNode = u[i];
				u[i].nodeRight = &tmpNode;
				u[counter] = *u[i].nodeRight;
				u[counter].x_index++;
				//Make sure to do something about the global index here! FIX!
				counter++;
			}

			//Up
			if(isInGrid(savedNodes[i].x_index, savedNodes[i].y_index+1) == false){
				Node tmpNode;
				tmpNode = u[i];
				u[i].nodeAbove = &tmpNode;
				u[counter] = *u[i].nodeAbove;
				u[counter].y_index++;
				//Make sure to do something about the global index here! FIX!
				counter++;
			}

			//down
			if(isInGrid(savedNodes[i].x_index, savedNodes[i].y_index-1) == false){
				Node tmpNode;
				tmpNode = u[i];
				u[i].nodeBelow = &tmpNode;
				u[counter] = *u[i].nodeBelow;
				u[counter].y_index--;
				//Make sure to do something about the global index here! FIX!
				counter++;
			}

			//Don't forget to create the corners as well.
			for (int x_mod = -1; x_mod <= 1; x_mod += 2)
			{
				for (int y_mod = -1; y_mod <= 1; y_mod += 2)
				{
					if(isInGrid(savedNodes[i].x_index+x_mod, savedNodes[i].y_index+y_mod) == false){
						u[counter] = u[i];
						u[counter].x_index += x_mod;
						u[counter].y_index += y_mod;
						//std::cout<<counter<<std::endl;
						counter++;
					}
				}
			}
		}

		//Make sure everyone finds thier neigbours.
		for (int i = 0; i < counter+1; ++i)
		{
			if(u[i].nodeLeft == NULL)
			{
				this->u[i].nodeLeft = findNode(u[i].x_index-1, u[i].y_index);
			}
			if (u[i].nodeRight == NULL)
			{
				this->u[i].nodeRight = findNode(u[i].x_index+1, u[i].y_index);
			}
			if(u[i].nodeAbove == NULL)
			{
				this->u[i].nodeAbove = findNode(u[i].x_index, u[i].y_index+1);
			}
			if (u[i].nodeBelow == NULL)
			{
				this->u[i].nodeBelow = findNode(u[i].x_index, u[i].y_index-1);
			}
		}

		this->len = counter+1;
		int *x_loc, *y_loc, *y_glo, *x_glo;
		x_loc =(int*) malloc(sizeof(int));
		y_loc =(int*) malloc(sizeof(int));
		x_glo =(int*) malloc(sizeof(int));
		y_glo =(int*) malloc(sizeof(int));
		for (int i = 0; i < len; ++i)
		{
			std::cout<<i<<std::endl;
			*x_loc = this->u[i].x_index;
			*y_loc = u[i].y_index;

			local2global(x_loc, y_loc, x_glo, y_glo);

			u[i].x_index_global = *x_glo;
			u[i].y_index_global = *y_glo;
		}
		free(x_loc);
		free(y_loc);
		free(x_glo);
		free(y_glo);
		//std::cout<<"fill rate: "<<(float)counter/(9.0f*numberOfPoints)<<std::endl;

	}

	AdaptiveGrid(){}

	AdaptiveGrid(const int layerIn ,const int numLayersIn, const int origo_x_in, const int origo_y_in, 
		Node* savedNodesIn , const int numberOfPointsIn, const datatype h_in){

		this->h = h_in;

		this->layerNr = layerIn;
		this->numOfLayers = numLayersIn;
		this->origo_y = origo_y_in;
		this->origo_x = origo_x_in;
		this->setupGrid(savedNodesIn, numberOfPointsIn);
		boundaryIndex = (int*) malloc(1);
		boundaryVals = (Node*) malloc(1);
		b = (Node*) malloc(1);
	}

	Node interpolateGhostPointFromGlobal(int x_glo, int y_glo){
		std::cout<<"GLLLLLLLLLLLLLLLLLLLLLLOBAL!!!!!"<<std::endl;
		std::cout<<"Llllllllllllllllllayer: "<<this->layerNr<<std::endl;
		assert(this->layerNr != 1);

		assert( this->coarserGrid != NULL );
		assert( this->coarserGrid->numOfLayers == this->numOfLayers );

		Node outNode;

		//int x_loc_mod = x_glo / (1<<(numOfLayers-coarserGrid->layerNr));
		//int y_loc_mod = y_glo / (1<<(numOfLayers-coarserGrid->layerNr));

		//std::cout<<"x_loc_mod: "<<x_loc_mod<<" y_loc_mod: "<<y_loc_mod<<std::endl;

		if(x_glo % (1<<(numOfLayers-coarserGrid->layerNr)) == 0){
			int y1 = y_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int y2 = y_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			std::cout<<"x"<<std::endl;
			outNode =  (coarserGrid->findNodeFromGlobalIndex(x_glo,y1)+coarserGrid->findNodeFromGlobalIndex(x_glo,y2))/2.0f;
		}
		else if(y_glo % (1<<(numOfLayers-coarserGrid->layerNr)) == 0){
			int x1 = x_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int x2 = x_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			std::cout<<"y"<<std::endl;
			outNode = (coarserGrid->findNodeFromGlobalIndex(x1,y_glo)+coarserGrid->findNodeFromGlobalIndex(x2,y_glo))/2.0f;
		}
		else{
			int x1 = x_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int x2 = x_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int y1 = y_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int y2 = y_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			std::cout<<"xy"<<std::endl;
			//std::cout<<x_loc_mod<<" "<<y_loc_mod<<std::endl;
			std::cout<<x_glo<<" "<<y_glo<<std::endl;
			std::cout<<(1<<(numOfLayers-coarserGrid->layerNr))<<std::endl;
			std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
			outNode = (coarserGrid->findNodeFromGlobalIndex(x1,y1)+coarserGrid->findNodeFromGlobalIndex(x1,y2)+coarserGrid->findNodeFromGlobalIndex(x2,y2)+coarserGrid->findNodeFromGlobalIndex(x2,y1))/4.0f;
		}

		return outNode;
	}

	Node findNodeFromGlobalIndex(const int x, const int y){
		for (int i = 0; i < len; ++i)
		{
			if(u[i].x_index_global == x && u[i].y_index_global == y){
				return u[i];
			}
		}
		std::cout<<"x: "<<x<<" ; y: "<<y<<std::endl;
		//The node doesn't exist. Oh, well.
		return this->interpolateGhostPointFromGlobal(x, y);
	}

	Node interpolateGhostPoint(int x_loc, int y_loc){
		assert( this->coarserGrid != NULL );
		assert( this->coarserGrid->numOfLayers == this->numOfLayers );

		int *x_glo, *y_glo;
		x_glo =(int*) malloc(sizeof(int));
		y_glo =(int*) malloc(sizeof(int));
		local2global(&x_loc, &y_loc, x_glo, y_glo);

		Node outNode;

		if(x_loc % 2 == 0){
			int y1 = *y_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int y2 = *y_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			std::cout<<"x"<<std::endl;
			outNode =  (coarserGrid->findNodeFromGlobalIndex(*x_glo,y1)+coarserGrid->findNodeFromGlobalIndex(*x_glo,y2))/2.0f;
		}
		else if(y_loc %2 == 0){
			int x1 = *x_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int x2 = *x_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			std::cout<<"y"<<std::endl;
			outNode = (coarserGrid->findNodeFromGlobalIndex(x1,*y_glo)+coarserGrid->findNodeFromGlobalIndex(x2,*y_glo))/2.0f;
		}
		else{
			int x1 = *x_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int x2 = *x_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int y1 = *y_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int y2 = *y_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			std::cout<<"xy"<<std::endl;
			std::cout<<x_loc<<" "<<y_loc<<std::endl;
			std::cout<<*x_glo<<" "<<*y_glo<<std::endl;
			std::cout<<(1<<(numOfLayers-coarserGrid->layerNr))<<std::endl;
			std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
			outNode = (coarserGrid->findNodeFromGlobalIndex(x1,y1)+coarserGrid->findNodeFromGlobalIndex(x1,y2)+coarserGrid->findNodeFromGlobalIndex(x2,y2)+coarserGrid->findNodeFromGlobalIndex(x2,y1))/4.0f;
		}

		free(x_glo);
		free(y_glo);

		return outNode;

	}

	void setBoundary(){
		int counter = 0;
		int x, y;
		boundaryIndex = (int*) malloc(2*lenOfBoundary*sizeof(int));
		boundaryVals = (Node*) malloc(lenOfBoundary*sizeof(Node));

		for (int i = 0; i < len; ++i)
		{
			x = u[i].x_index + 1;
			y = u[i].y_index;

			if(isInGrid(x,y)==false){
				boundaryIndex[counter] = x;
				counter++;
				boundaryIndex[counter] = y;
				counter++;
			}
			x = u[i].x_index - 1;
			y = u[i].y_index;

			if(isInGrid(x,y)==false){
				boundaryIndex[counter] = x;
				counter++;
				boundaryIndex[counter] = y;
				counter++;
			}
			x = u[i].x_index;
			y = u[i].y_index + 1;

			if(isInGrid(x,y)==false){
				boundaryIndex[counter] = x;
				counter++;
				boundaryIndex[counter] = y;
				counter++;
			}
			x = u[i].x_index;
			y = u[i].y_index - 1;

			if(isInGrid(x,y)==false){
				boundaryIndex[counter] = x;
				counter++;
				boundaryIndex[counter] = y;
				counter++;
			}
		}

		for (int i = 0; i < lenOfBoundary; ++i)
		{
			boundaryVals[i] = interpolateGhostPoint(boundaryIndex[2*i], boundaryIndex[2*i+1]);
		}


	}

	void updateBFromBoundary(){
		Node boundary;
		for (int i = 0; i < len; ++i)
		{
			if (u[i].nodeAbove == NULL)
			{
				std::cout<<"above"<<std::endl;
				boundary = getFromBoundary(u[i].x_index, u[i].y_index + 1);
				b[i] = b[i] + boundary/(-h);
			}
			if (u[i].nodeBelow == NULL)
			{
				std::cout<<"below"<<std::endl;
				boundary = getFromBoundary(u[i].x_index, u[i].y_index - 1);
				b[i] = b[i] + boundary/(-h);
			}
			if (u[i].nodeRight == NULL)
			{
				std::cout<<"right"<<std::endl;
				boundary = getFromBoundary(u[i].x_index +1, u[i].y_index);
				b[i] = b[i] + boundary/(-h);
			}
			if (u[i].nodeLeft == NULL)
			{
				std::cout<<"left"<<std::endl;
				boundary = getFromBoundary(u[i].x_index -1 , u[i].y_index);
				b[i] = b[i] + boundary/(-h);
			}
		}
	}

	void updateBoundary(){
		this->resetBoundaryLength();
		this->setBoundaryLength();
		free(boundaryIndex);
		free(boundaryVals);
		free(b);
		this->setBoundary();
		b = (Node*) calloc(len,sizeof(Node));
		updateBFromBoundary();
	}

};




void multigrid(int k, AdaptiveGrid grid, int pre, int sol, int post){
	for (int i = 0; i < pre; ++i)
	{

		for (int j = 0; j < grid.len; ++j)
		{
			std::cout<<"u[i]: "<<grid.u[j].x_index<<" ; "<<grid.u[j].y_index<<" , "<<grid.u[j].stream<<std::endl;
		}
		std::cout<<std::endl;

		grid.jacobiSmootherLaplacianStream();
	}
}












int main(int argc, char const *argv[])
{
	Node* fromWavelet;
	fromWavelet = (Node*) malloc((9+3+2)*sizeof(Node));

	datatype h = 0.125f;
	for (int x = 0; x < 3; x += 1)
	{
		for (int y = 0; y < 3; y += 1)
		{
			fromWavelet[x*3+y] = Node(x*h, y*h, x,y, 0.0f,0.0f,0.0f);
		}
	}

	fromWavelet[9]  = Node(0.25,0.25,1,1, 1.0f,0.1f,0.0f);
	fromWavelet[10] = Node(0.25,0.75,1,3, 1.5f,0.2f,0.0f);
	fromWavelet[11] = Node(0.75,0.75,3,3, 2.0f,0.3f,0.0f);

	fromWavelet[12] = Node(0.125+0.25,0.5+0.125, 1,1, 0.1f,0.4f,0.0f);
	fromWavelet[13] = Node(0.125+0.50,0.5+0.125, 3,1, 0.1f,0.5f,0.0f);

	

	AdaptiveGrid fine(3,3,2,4, &fromWavelet[12], 2, 4*h);
	AdaptiveGrid middle(2,3,0,0, &fromWavelet[9], 3, 2*h);
	AdaptiveGrid coarse(1,3,0,0, &fromWavelet[0], 9, h);

	fine.coarserGrid = &middle;
	middle.coarserGrid = &coarse;
	//middle.setBoundaryLength();
	//std::cout<<middle.lenOfBoundary<<std::endl;

	std::cout<<std::endl;
	for (int i = 0; i < coarse.len; ++i)
	{
		//std::cout<<std::endl;
		//std::cout<<"u[i]: "<<coarse.u[i].x_index<<" ; "<<coarse.u[i].y_index<<std::endl;
		//std::cout<<"u[i]: "<<coarse.u[i].x_index_global<<" ; "<<coarse.u[i].y_index_global<<std::endl;
		//std::cout<<std::endl;
	}
	
	fine.updateBoundary();
	middle.updateBoundary();
	//fine.setBoundaryLength();
	//std::cout<<fine.lenOfBoundary<<std::endl;
	for (int i = 0; i < fine.len; ++i)
	{
		std::cout<<"u[i]: "<<fine.u[i].x_index<<" ; "<<fine.u[i].y_index<<" , "<<fine.u[i].stream<<std::endl;
	}
	multigrid(0,fine, 12,6,2);
	//assert(0);
	for (int i = 0; i < fine.len; ++i)
	{
		std::cout<<"u[i]: "<<fine.u[i].x_index<<" ; "<<fine.u[i].y_index<<" , "<<fine.u[i].stream<<std::endl;
	}
	std::cout<<coarse.len<<std::endl;
	/*
	middle.interpolateU(fine);

	std::cout<<"Fine: "<<std::endl;
	for (int i = 0; i < fine.len; ++i)
	{
		std::cout<<"x:"<<fine.u[i].x_index<<", y:"<<fine.u[i].y_index<<", vort:"<<middle.u[i].vort<<std::endl;
	}

	std::cout<<"Middle: "<<std::endl;
	for (int i = 0; i < middle.len; ++i)
	{
		std::cout<<"x:"<<middle.u[i].x_index<<", y:"<<middle.u[i].y_index<<", vort:"<<middle.u[i].vort<<std::endl;
	}
	*/

	free(fromWavelet);

	/*
	//AdaptiveGrid grid;

	//AdaptiveGrid coarse;
	int coarse_n = 9;

	int n = 21;

	grid.u = (Node*) malloc(n*sizeof(Node));
	grid.b = (Node*) malloc(n*sizeof(Node));
	grid.h = 0.25f;

	coarse.u = (Node*) malloc(coarse_n*sizeof(Node));
	coarse.b = (Node*) malloc(coarse_n*sizeof(Node));
	coarse.h = 0.5f;

	int counter = 0;
	for (int x = 0; x < 5; ++x)
	 {
	 	for (int y = 0; y < 5; ++y)
	 	{
	 		if(y>2 && x<2){
	 			//do nothing
	 		}
	 		else{
	 			grid.u[counter] = Node(grid.h*y, grid.h*y, x, y,0.1f,0.2f,0.3f);
	 			counter++;
	 		}
	 	}
	 } 
	 grid.len = counter;
	

	 for (int i = 0; i < grid.len; ++i)
	 {
	 	grid.setB(i, 0.0f);
	 }


	 for (int i = 0; i < grid.len; ++i)
	 {
	 	std::cout<<grid.u[i].x_index << " ; " <<grid.u[i].y_index << " ; " <<grid.u[i].stream << "\n" ;
	 }

		grid.jacobiSmootherLaplacianStream();
		grid.jacobiSmootherLaplacianStream();
		grid.jacobiSmootherLaplacianStream();
		grid.jacobiSmootherLaplacianStream();


	std::cout<<std::endl;
	for (int i = 0; i < grid.len; ++i)
	 {
	 	std::cout<<grid.u[i].x_index << " ; " <<grid.u[i].y_index << " ; " <<grid.u[i].stream << "\n" ;
	 }


	free(grid.u);
	free(grid.b);

	free(coarse.u);
	free(coarse.b);
	*/

	return 0;
}
