
#include "adaptive_grid.h"

	bool AdaptiveGrid::isInBoundary(const int x, const int y){
		for (int i = 0; i < lenOfBoundary; ++i)
		{
			if(boundaryIndex[2*i]==x && boundaryIndex[2*i+1]==y){
				return true;
			}
		}
		return false;
	}

	void AdaptiveGrid::setB(const datatype val, const int i){
		assert(i<len && i>=0);
		b[i].stream = val;
	}

	datatype AdaptiveGrid::getB(const int i){
		assert(i<len && i>=0);
		return b[i].stream;
	}

	datatype AdaptiveGrid::getFromBoundary(const int x, const int y){
		for (int i = 0; i < lenOfBoundary; ++i)
		{
			if(boundaryIndex[2*i] == x && boundaryIndex[2*i+1] == y){
				return boundaryVals[i].stream;
			}
		}
		assert(0);
		return boundaryVals[0].stream;
	}

	void AdaptiveGrid::calculateErrorLaplacian(){
		datatype sum;
		for(int i=0; i<len; i++){
			sum = 0.0f;
			for(int j=0; j<len; j++){
				sum += getLaplacianStream(j,i)*getU(j);
			}
			d[i].stream = b[i].stream-sum;
		}
	}

	datatype AdaptiveGrid::getLaplacianStream(const int index1, const int index2){

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
			 		return getFromBoundary(u[index2].x_index, u[index2].y_index)/h;
			 	}
			 	else{
			 		return 1.0f/(h*h);
			 	}
			}
		}
	}

	void AdaptiveGrid::resetBoundaryLength(){
		lenOfBoundary = 0;
	}

	bool AdaptiveGrid::isInGrid(const int x, const int y){
		for (int i = 0; i < len; ++i)
		{
			if(u[i].x_index==x && u[i].y_index==y){
				return true;
			}
		}
		//std::cout<<x<<" ; "<<y<<std::endl;
		return false;
	}

	void AdaptiveGrid::setBoundaryLength(){
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


	void AdaptiveGrid::setU(const datatype val, const int i){
		assert(i<len && i>=0);
		u[i].stream = val;
	}

	datatype AdaptiveGrid::getU(const int i){
		assert(i<len && i>=0);
		return u[i].stream;
	}

	void AdaptiveGrid::jacobiSmootherLaplacianStream(){
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

	void AdaptiveGrid::restrictUtoU(AdaptiveGrid* coarse){
		for (int i = 0; i < coarse->len; ++i)
		{
			for (int j = 0; j < this->len; ++j)
			{
				if( this->u[j].x_index_global == coarse->u[i].x_index_global 
					&& this->u[j].y_index_global == coarse->u[i].y_index_global )
				{
					coarse->u[i].stream = this->u[j].stream;
				}
			}
		}
	}

	void AdaptiveGrid::restrictDtoB(AdaptiveGrid* coarse){
		for (int i = 0; i < coarse->len; ++i)
		{
			for (int j = 0; j < this->len; ++j)
			{
				if( this->d[j].x_index_global == coarse->b[i].x_index_global 
					&& this->d[j].y_index_global == coarse->b[i].y_index_global )
				{
					coarse->b[i].stream = this->d[j].stream;
				}
			}
		}
	}

	void AdaptiveGrid::restrictDtoD(AdaptiveGrid* coarse){
		for (int i = 0; i < coarse->len; ++i)
		{
			for (int j = 0; j < this->len; ++j)
			{
				if( this->d[j].x_index_global == coarse->d[i].x_index_global 
					&& this->d[j].y_index_global == coarse->d[i].y_index_global )
				{
					coarse->d[i].stream = this->d[j].stream;
				}
			}
		}
	}

	Node* AdaptiveGrid::findNode(const int ind_x, const int ind_y){

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

	Node * AdaptiveGrid::findNodeD(const int ind_x, const int ind_y){

		for (int i = 0; i < len; ++i)
		{
			if(d[i].x_index == ind_x && d[i].y_index == ind_y){
				return &d[i];
			}
		}
		//std::cout<<"layer: "<<layerNr<<", findNodeD; "<<ind_x<<"; "<<ind_y<<std::endl;
		//assert(0);
		return NULL;

	}

	void AdaptiveGrid::interpolateU(AdaptiveGrid *fine){
		//This doesn't work. FIX!
		/*for (int i = 0; i < fine->len; ++i)
		{
			for (int j = 0; j < this->len; ++j)
			{
				if( this->u[j].x_index_global == fine->u[i].x_index_global 
					&& this->u[j].y_index_global == fine->u[i].y_index_global )
				{
					fine->u[i].stream = this->u[j].stream;
				}
			}
		}*/

			fine->findNeighboursU();
			for (int i = 0; i < fine->numOfPointsChosen; ++i)
			{

				//std::cout<<"i: "<<i<<std::endl;
				Node *middleNode = fine->findNode(fine->u[i].x_index, fine->u[i].y_index);
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


	void AdaptiveGrid::interpolateD(AdaptiveGrid *fine){
		/*
		for (int i = 0; i < fine->len; ++i)
		{
			for (int j = 0; j < this->len; ++j)
			{
				if( this->d[j].x_index_global == fine->d[i].x_index_global 
					&& this->d[j].y_index_global == fine->d[i].y_index_global )
				{
					fine->d[i].stream = this->d[j].stream;
				}
			}
		}*/

		fine->findNeighboursD();

		for (int i = 0; i < fine->numOfPointsChosen; ++i)
		{

			//std::cout<<"i: "<<i<<std::endl;
			//std::cout<<"nu kör vi! x: "<<fine->u[i].x_index<<"; "<<fine->u[i].y_index<<std::endl;
			Node *middleNode = fine->findNodeD(fine->u[i].x_index, fine->u[i].y_index);
			Node *tmpNodePtr;

			//Interpolate to the left
			tmpNodePtr = middleNode->nodeLeft;
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

	void AdaptiveGrid::findNeighboursD(){
		//Make sure everyone finds thier neigbours.
		for (int i = 0; i < len; ++i)
		{
			this->d[i].nodeLeft = findNodeD(d[i].x_index-1, d[i].y_index);
			this->d[i].nodeRight = findNodeD(d[i].x_index+1, d[i].y_index);
			this->d[i].nodeAbove = findNodeD(d[i].x_index, d[i].y_index+1);
			this->d[i].nodeBelow = findNodeD(d[i].x_index, d[i].y_index-1);
		}
	}

	void AdaptiveGrid::findNeighboursU(){
		//Make sure everyone finds thier neigbours.
		for (int i = 0; i < len; ++i)
		{
			this->u[i].nodeLeft = findNode(u[i].x_index-1, u[i].y_index);
			this->u[i].nodeRight = findNode(u[i].x_index+1, u[i].y_index);
			this->u[i].nodeAbove = findNode(u[i].x_index, u[i].y_index+1);
			this->u[i].nodeBelow = findNode(u[i].x_index, u[i].y_index-1);
		}
	}


	void AdaptiveGrid::setNeighbours(){

		for (int i = 0; i < numOfPointsChosen; ++i)
		{
			//Up
			(pointsChosenByWavelet[i]->nodeAbove) = findNode(pointsChosenByWavelet[i]->x_index, pointsChosenByWavelet[i]->y_index + 1);

			//down
			(pointsChosenByWavelet[i]->nodeBelow) = findNode(pointsChosenByWavelet[i]->x_index, pointsChosenByWavelet[i]->y_index - 1);
			
			//Left
			(pointsChosenByWavelet[i]->nodeLeft)  = findNode(pointsChosenByWavelet[i]->x_index-1, pointsChosenByWavelet[i]->y_index);

			//Right
			(pointsChosenByWavelet[i]->nodeRight) = findNode(pointsChosenByWavelet[i]->x_index+1, pointsChosenByWavelet[i]->y_index);
			

		}
	}

	void AdaptiveGrid::local2global(const int* x_loc, const int* y_loc, int* x_glo, int* y_glo ){

		assert( this->layerNr > 0 );
		assert( this->layerNr <= this->numOfLayers );

		*y_glo = *y_loc*(1<<(numOfLayers-layerNr))+this->origo_y;
		*x_glo = *x_loc*(1<<(numOfLayers-layerNr))+this->origo_x;
	}

	void AdaptiveGrid::setupGrid(Node* savedNodes, const int numberOfPoints){
		//free(pointsChosenByWavelet);
		pointsChosenByWavelet = (Node**) malloc(numberOfPoints*sizeof(Node*));
		this->numOfPointsChosen = numberOfPoints;

		//free(u);
		//9 is the maximum. It's probably waaaay too large. OPT!
		u = (Node*) malloc(9*numberOfPoints*sizeof(Node));
		len = numberOfPoints;

		int counter = numberOfPoints;

		//Take the points given by 
		for (int i = 0; i < numberOfPoints; ++i)
		{
			pointsChosenByWavelet[i] = &savedNodes[i];
			u[i] = savedNodes[i];
		}

		for (int i = 0; i < numberOfPoints; ++i){
			//Left
			if(isInGrid(savedNodes[i].x_index-1, savedNodes[i].y_index) == false){
				//assert(i != 0);
				Node tmpNode;
				tmpNode = u[i];
				u[i].nodeLeft = &tmpNode;
				u[counter] = *u[i].nodeLeft;
				u[counter].x_index--;
				len++;
				counter++;
			}
			//Right
			if(isInGrid(savedNodes[i].x_index+1, savedNodes[i].y_index) == false){
				Node tmpNode;
				tmpNode = u[i];
				u[i].nodeRight = &tmpNode;
				u[counter] = *u[i].nodeRight;
				u[counter].x_index++;
				len++;
				counter++;
			}
			//Up
			if(isInGrid(savedNodes[i].x_index, savedNodes[i].y_index+1) == false){
				Node tmpNode;
				tmpNode = u[i];
				u[i].nodeAbove = &tmpNode;
				u[counter] = *u[i].nodeAbove;
				u[counter].y_index++;
				len++;
				counter++;
			}
			//Down
			if(isInGrid(savedNodes[i].x_index, savedNodes[i].y_index-1) == false){
				Node tmpNode;
				tmpNode = u[i];
				u[i].nodeBelow = &tmpNode;
				u[counter] = *u[i].nodeBelow;
				u[counter].y_index--;
				len++;
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
						len++;
						counter++;
					}
				}
			}
		}

		//Make sure everyone finds thier neigbours.
		for (int i = 0; i < counter; ++i)
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

		assert(len == counter);

		this->len = counter; //+1 ?
		int *x_loc, *y_loc, *y_glo, *x_glo;
		x_loc =(int*) malloc(sizeof(int));
		y_loc =(int*) malloc(sizeof(int));
		x_glo =(int*) malloc(sizeof(int));
		y_glo =(int*) malloc(sizeof(int));
		for (int i = 0; i < len; ++i)
		{
			//std::cout<<i<<std::endl;
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

	void AdaptiveGrid::setupCoarsestGrid(Node* savedNodes, const int numberOfPoints){
		/*
			The coarsest grid has five chosen points! It always has five chosen points.
			ALWAYS!
			The four corners and the middle. Okay? Okay.

			(I am aware that this function is extremely ad hoc, and all.
			But it works, doesn't it?)
		*/
		
		//std::cout<<"num points:"<<numberOfPoints<<std::endl;

		pointsChosenByWavelet = (Node**) malloc(numberOfPoints*sizeof(Node*));
		this->numOfPointsChosen = numberOfPoints;

		u = (Node*) malloc(9*sizeof(Node));
		len = 9;

		int *x_local, *y_local;
		x_local = (int*) calloc(1, sizeof(int));
		y_local = (int*) calloc(1, sizeof(int));

		for (int i = 0; i < numberOfPoints; ++i)
		{
			pointsChosenByWavelet[i] = &savedNodes[i];
			u[i] = savedNodes[i];
			global2local(u[i].x_index_global, u[i].y_index_global, x_local, y_local);
			u[i].x_index = *x_local;
			u[i].y_index = *y_local;
		}
		free(x_local);
		free(y_local);
		
		std::cout<<"numPoints (coarsest grid): "<<numberOfPoints<<std::endl;

		int counter = numberOfPoints;
		len = counter;

		//std::cout<<"lay_Nr: "<<layerNr<<std::endl;
		//std::cout<<"lay_Nr_max: "<<numOfLayers<<std::endl;
		//std::cout<<"origo_x: "<<origo_x<<std::endl;
		//std::cout<<"origo_y: "<<origo_y<<std::endl;
		//std::cout<<"counter (!): "<<counter<<std::endl;
		//for (int i = 0; i < 5; ++i)
		//{
		//	std::cout<<"x_index: "<<u[i].x_index<<"  y_index: "<<u[i].y_index<<std::endl;
		//	std::cout<<"x_index_glob: "<<u[i].x_index_global<<"  y_index_glob: "<<u[i].y_index_global<<std::endl;
		//}

		//Start with node (0,1)
		Node *tmpNodePtr;

		for (int i = 0; i < numberOfPoints; ++i)
		{
			if(u[i].x_index == 0 && u[i].y_index == 0){
				tmpNodePtr = &u[i];
				break;
			}
			assert(i != numberOfPoints -1);
		}
		Node tmp = *tmpNodePtr;
		tmp.y_index++;
		if (findNodeGeneral(u, tmp.x_index, tmp.y_index) == NULL)
		{
			//std::cout<<"test 1\n";
			u[counter] = tmp;
			counter++;
		}		

		//Then (1,0)
		tmp = *tmpNodePtr;
		tmp.x_index++;
		if (findNodeGeneral(u, tmp.x_index, tmp.y_index) == NULL)
		{
			//std::cout<<"test 2\n";
			u[counter] = tmp;
			counter++;
		}

		//Now (2,1)
		for (int i = 0; i < numberOfPoints; ++i)
		{
			if(u[i].x_index == 2 && u[i].y_index == 2){
				tmpNodePtr = &u[i];
				break;
			}
			assert(i != numberOfPoints -1);
		}
		tmp = *tmpNodePtr;
		tmp.y_index--;
		if (findNodeGeneral(u, tmp.x_index, tmp.y_index) == NULL)
		{
			//std::cout<<"test 3\n";
			u[counter] = tmp;
			counter++;
		}	

		
		//Then (1,2)
		tmp = *tmpNodePtr;
		tmp.x_index--;
		//std::cout<<"count "<<counter<<std::endl;

		if (findNodeGeneral(u, tmp.x_index, tmp.y_index) == NULL)
		{
			//std::cout<<"test 4\n";
			u[counter] = tmp;
			counter++;
		}

		//for (int i = 0; i < 5; ++i)
		//{
		//	std::cout<<"x_index: "<<u[i].x_index<<"  y_index: "<<u[i].y_index<<std::endl;
		//	std::cout<<"x_index_glob: "<<u[i].x_index_global<<"  y_index_glob: "<<u[i].y_index_global<<std::endl;
		//}

		//Make sure everyone finds thier neigbours.
		for (int i = 0; i < counter; ++i)
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
		//std::cout<<"count "<<counter<<std::endl;
		assert(9 == counter);
		len = counter;
		int *x_loc, *y_loc, *y_glo, *x_glo;
		x_loc =(int*) malloc(sizeof(int));
		y_loc =(int*) malloc(sizeof(int));
		x_glo =(int*) malloc(sizeof(int));
		y_glo =(int*) malloc(sizeof(int));
		for (int i = 0; i < len; ++i)
		{
			//std::cout<<i<<std::endl;
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

		std::cout<<"setupCoarsestGrid is done\n";

		//len = counter;
		//std::cout<<"fill rate: "<<(float)counter/(9.0f*numberOfPoints)<<std::endl;
	}

	Node* AdaptiveGrid::findNodeGeneral(Node* arr, const int ind_x, const int ind_y){


		/*
		bool tmp = false;

		for (int i = 0; i < len; ++i)
		{
			if(arr[i].x_index == ind_x && arr[i].y_index == ind_y){
				//return &arr[i];
				if(tmp == true){
					assert(0);
				}
				tmp = true;
			}
		}
		*/








		for (int i = 0; i < len; ++i)
		{
			if(arr[i].x_index == ind_x && arr[i].y_index == ind_y){
				//std::cout<<"i "<<i<<std::endl;
				return &arr[i];
				
			}
		}
		return NULL;

	}		

Node* AdaptiveGrid::findGlobNodeGeneral(Node* arr, const int ind_x, const int ind_y, const int lenIn){

		for (int i = 0; i < lenIn; ++i)
		{
			if(arr[i].x_index_global == ind_x && arr[i].y_index_global == ind_y){
				return &arr[i];
				
			}
		}
		return NULL;

	}

	void AdaptiveGrid::findNeighbours(Node * arr){
		//Make sure everyone finds thier neigbours.
		for (int i = 0; i < len; ++i) 
		{
			arr[i].nodeLeft =  findNodeGeneral(arr, arr[i].x_index-1, arr[i].y_index);
			arr[i].nodeRight = findNodeGeneral(arr, arr[i].x_index+1, arr[i].y_index);
			arr[i].nodeAbove = findNodeGeneral(arr, arr[i].x_index, arr[i].y_index+1);
			arr[i].nodeBelow = findNodeGeneral(arr, arr[i].x_index, arr[i].y_index-1);
		}
	}

	int AdaptiveGrid::removeDumbPoints(Node* savedNodesIn, Node* savedNodes, const int numberOfPointsIn){
		int counter = 0;

		int x_mod = 0;
		int y_mod = 0;
		for (int i = 0; i < numberOfPointsIn; ++i)
		{
			x_mod = 0;
			y_mod = 0;
			assert( savedNodesIn[i].x_index_global != -1 && savedNodesIn[i].y_index_global != -1 );
			if( (savedNodesIn[i].x_index_global %2) != (savedNodesIn[i].x_index_global %2 == 0) )
			{
				if ((savedNodesIn[i].x_index_global %2) == 1){
					if(savedNodesIn[i].y_index_global>0){
						y_mod = -1;
					}
					else{
						y_mod = 1;
					}
				}
				else{
					if(savedNodesIn[i].x_index_global>0){
						x_mod = -1;
					}
					else{
						x_mod = 1;
					}
				}
			}
			
			if (findGlobNodeGeneral(savedNodes, savedNodesIn[i].x_index_global+x_mod, savedNodesIn[i].y_index_global+y_mod, counter+1) == NULL)
			{
				savedNodes[counter] = savedNodesIn[i];
				savedNodes[counter].x_index_global += x_mod;
				savedNodes[counter].y_index_global += y_mod;
				counter++;
			}

		}

		return counter;
	}

	AdaptiveGrid::AdaptiveGrid(){
		coarserGrid = NULL;
		finerGrid = NULL;
	}

	AdaptiveGrid::AdaptiveGrid(const int layerIn, const int numLayersIn, const int origo_x_in, const int origo_y_in, 
		Node* savedNodesIn, const int numberOfPointsIn, const datatype h_in){

		this->h = h_in;

		coarserGrid = NULL;
		finerGrid = NULL;

		this->layerNr = layerIn;
		this->numOfLayers = numLayersIn;
		this->origo_y = origo_y_in;
		this->origo_x = origo_x_in;

		if(numberOfPointsIn > 0){

			Node* savedNodes;
			savedNodes = (Node*) malloc(numberOfPointsIn*sizeof(Node));

			int numberOfPoints = removeDumbPoints(savedNodesIn, savedNodes, numberOfPointsIn);
			//free(savedNodesIn);

			//assert(0);

			if(this->layerNr == 1){
				this->setupCoarsestGrid(savedNodes, numberOfPoints);
			}
			else{
				this->setupGrid(savedNodes, numberOfPoints);
			}

			b = (Node*) calloc(len, sizeof(Node));
			d = (Node*) calloc(len, sizeof(Node));
			w = (Node*) calloc(len, sizeof(Node));

			Node * u_tmp;
			u_tmp = (Node*) malloc(len*sizeof(Node));

			for (int i = 0; i < len; ++i)
			{
				u_tmp[i] = u[i];

				b[i] = u[i];
				b[i].stream = 0;
				b[i].vort = 0;

				d[i] = u[i];
				d[i].stream = 0;
				d[i].vort = 0;

				w[i] = u[i];
				w[i].stream = 0;
				w[i].vort = 0;
			}

			free(u);
			u = u_tmp;

			findNeighbours(b);
			findNeighbours(w);
			findNeighbours(d);
			findNeighbours(u);
		}
	}

	void AdaptiveGrid::setupFromVector(){
		int pointsInVec = this->vec.size();
		Node* savedNodes;
		savedNodes = (Node*) malloc(pointsInVec*sizeof(Node));
		for (int i = 0; i < pointsInVec; ++i)
		{
			savedNodes[i] = this->vec[i];
		}

		//int numberOfPoints = removeDumbPoints(savedNodesIn, savedNodes, numberOfPointsIn);
			//free(savedNodesIn);

		//assert(0);

		if(this->layerNr == 1){
			this->setupCoarsestGrid(savedNodes, pointsInVec);
		}
		else{
			this->setupGrid(savedNodes, pointsInVec);
		}

		b = (Node*) calloc(len, sizeof(Node));
		d = (Node*) calloc(len, sizeof(Node));
		w = (Node*) calloc(len, sizeof(Node));

		Node * u_tmp;
		u_tmp = (Node*) malloc(len*sizeof(Node));

		for (int i = 0; i < len; ++i)
		{
			u_tmp[i] = u[i];

			b[i] = u[i];
			b[i].stream = 0;
			b[i].vort = 0;

			d[i] = u[i];
			d[i].stream = 0;
			d[i].vort = 0;

			w[i] = u[i];
			w[i].stream = 0;
			w[i].vort = 0;
		}

		free(u);
		free(savedNodes);
		u = u_tmp;

		findNeighbours(b);
		findNeighbours(w);
		findNeighbours(d);
		findNeighbours(u);
	}

	Node AdaptiveGrid::interpolateGhostPointFromGlobal(int x_glo, int y_glo){
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
			//std::cout<<"x"<<std::endl;
			outNode =  (coarserGrid->findNodeFromGlobalIndex(x_glo,y1)+coarserGrid->findNodeFromGlobalIndex(x_glo,y2))/2.0f;
		}
		else if(y_glo % (1<<(numOfLayers-coarserGrid->layerNr)) == 0){
			int x1 = x_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int x2 = x_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			//std::cout<<"y"<<std::endl;
			outNode = (coarserGrid->findNodeFromGlobalIndex(x1,y_glo)+coarserGrid->findNodeFromGlobalIndex(x2,y_glo))/2.0f;
		}
		else{
			int x1 = x_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int x2 = x_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int y1 = y_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
			int y2 = y_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
			//std::cout<<"xy"<<std::endl;
			//std::cout<<x_loc_mod<<" "<<y_loc_mod<<std::endl;
			//std::cout<<x_glo<<" "<<y_glo<<std::endl;
			//std::cout<<(1<<(numOfLayers-coarserGrid->layerNr))<<std::endl;
			//std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
			outNode = (coarserGrid->findNodeFromGlobalIndex(x1,y1)+coarserGrid->findNodeFromGlobalIndex(x1,y2)+coarserGrid->findNodeFromGlobalIndex(x2,y2)+coarserGrid->findNodeFromGlobalIndex(x2,y1))/4.0f;
		}

		return outNode;
	}

	Node AdaptiveGrid::findNodeFromGlobalIndex(const int x, const int y){
		for (int i = 0; i < len; ++i)
		{
			if(u[i].x_index_global == x && u[i].y_index_global == y){
				return u[i];
			}
		}
		//std::cout<<"x: "<<x<<" ; y: "<<y<<std::endl;
		//The node doesn't exist. Oh, well.
		return this->interpolateGhostPointFromGlobal(x, y);
	}

	Node AdaptiveGrid::interpolateGhostPoint(int x_loc, int y_loc){
		assert( this->coarserGrid != NULL );
		assert( this->coarserGrid->numOfLayers == this->numOfLayers );
		int *x_glo, *y_glo;
		x_glo =(int*) malloc(sizeof(int));
		y_glo =(int*) malloc(sizeof(int));
		local2global(&x_loc, &y_loc, x_glo, y_glo);
		Node outNode;

		if(*x_glo<0 || *y_glo<0 || *y_glo> (1<<numOfLayers) || *x_glo> (1<<numOfLayers)){
			free(x_glo);
			free(y_glo);
			return outNode;
		}
		else{

			if(x_loc % 2 == 0){
				int y1 = *y_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
				int y2 = *y_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
				//std::cout<<"x"<<std::endl;
				outNode =  (coarserGrid->findNodeFromGlobalIndex(*x_glo,y1)+coarserGrid->findNodeFromGlobalIndex(*x_glo,y2))/2.0f;
			}
			else if(y_loc %2 == 0){
				int x1 = *x_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
				int x2 = *x_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
				//std::cout<<"y"<<std::endl;
				outNode = (coarserGrid->findNodeFromGlobalIndex(x1,*y_glo)+coarserGrid->findNodeFromGlobalIndex(x2,*y_glo))/2.0f;
			}
			else{
				int x1 = *x_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
				int x2 = *x_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
				int y1 = *y_glo+(1<<(numOfLayers-coarserGrid->layerNr))/2;
				int y2 = *y_glo-(1<<(numOfLayers-coarserGrid->layerNr))/2;
				//std::cout<<"xy"<<std::endl;
				//std::cout<<x_loc<<" "<<y_loc<<std::endl;
				//std::cout<<*x_glo<<" "<<*y_glo<<std::endl;
				//std::cout<<(1<<(numOfLayers-coarserGrid->layerNr))<<std::endl;
				//std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
				outNode = (coarserGrid->findNodeFromGlobalIndex(x1,y1)+coarserGrid->findNodeFromGlobalIndex(x1,y2)+coarserGrid->findNodeFromGlobalIndex(x2,y2)+coarserGrid->findNodeFromGlobalIndex(x2,y1))/4.0f;
			}

			free(x_glo);
			free(y_glo);

			return outNode;
	}

	}

	void AdaptiveGrid::setBoundary(){
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

	void AdaptiveGrid::updateBFromBoundary(){
		datatype boundary;
		for (int i = 0; i < len; ++i)
		{
			if (u[i].nodeAbove == NULL)
			{
				//std::cout<<"above"<<std::endl;
				boundary = getFromBoundary(u[i].x_index, u[i].y_index + 1);
				b[i].stream += boundary/(-h);
			}
			if (u[i].nodeBelow == NULL)
			{
				//std::cout<<"below"<<std::endl;
				boundary = getFromBoundary(u[i].x_index, u[i].y_index - 1);
				b[i].stream += boundary/(-h);
			}
			if (u[i].nodeRight == NULL)
			{
				//std::cout<<"right"<<std::endl;
				boundary = getFromBoundary(u[i].x_index +1, u[i].y_index);
				b[i].stream += boundary/(-h);
			}
			if (u[i].nodeLeft == NULL)
			{
				//std::cout<<"left"<<std::endl;
				boundary = getFromBoundary(u[i].x_index -1 , u[i].y_index);
				b[i].stream += boundary/(-h);
			}
		}
	}

	/*
	void updateBoundary(){
		this->resetBoundaryLength();
		this->setBoundaryLength();
		free(boundaryIndex);
		free(boundaryVals);
		free(b);
		this->setBoundary();
		b = (Node*) calloc(len,sizeof(Node));
		updateBFromBoundary();
		for (int i = 0; i < len; ++i)
		{
			d[i] = u[i];

		}
		findNeighboursD();

	}

	void updateBoundaryCoarsestGrid(){
		this->resetBoundaryLength();
		free(boundaryIndex);
		free(boundaryVals);
		free(b);
		b = (Node*) calloc(len,sizeof(Node));
		for (int i = 0; i < len; ++i)
		{
			d[i] = u[i];
		}
		findNeighboursD();
	}
	*/

	void AdaptiveGrid::updateBFromFunction( func_t externalForceFunction){
		datatype x,y;
		for (int i = 0; i < this->len; ++i)
		{
			x = u[i].x;
			y = u[i].y;

			this->b[i].stream += externalForceFunction(x, y);
		}
	}

	void AdaptiveGrid::calculateRHS(){
		for (int i = 0; i < len; ++i)
		{
			datatype sum = 0.0f;
			for (int j = 0; j < len; ++j)
			{
				sum += getLaplacianStream(j,i)*u[j].stream;
			}
			b[i].stream = d[i].stream + sum;
		}
	}

	void AdaptiveGrid::global2local(const int x_glo, const int y_glo, int* x_loc, int* y_loc ){
		assert(this->layerNr > 0);
		assert(this->layerNr <= this->numOfLayers);
		assert(origo_x >= 0);
		assert(origo_y >= 0);

		*x_loc = x_glo;
		*x_loc = *x_loc-origo_x;
		*x_loc = *x_loc/(1<<(numOfLayers-layerNr));

		*y_loc = y_glo;
		*y_loc = *y_loc-origo_y;
		*y_loc = *y_loc/(1<<(numOfLayers-layerNr));
	}

	/*void AdaptiveGrid::calcVortFromStream(){

		const int N = grid->len;

		//Define the size of the current grid in a CUDA format
		dim3 block_size_1d(N);
		dim3 grid_size_1d((N+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK);


	}*/

