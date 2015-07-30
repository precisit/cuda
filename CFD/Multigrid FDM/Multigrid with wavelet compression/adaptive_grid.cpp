#include <assert.h>
#include "node.cpp"
class AdaptiveGrid{
public:

	Node* u;
    int len;
    Node* b;

    datatype* boundaryVals;
    int* boundaryIndex;
    int lenOfBoundary;

    Node** pointsChosenByWavelet; //An array of Node pointers
    int numOfPointsChosen;

    //adaptiveGrid* finerGrid;
    //adaptiveGrid* coarserGrid;
    datatype h;

	bool isInBoundary(const int x, const int y){

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

	datatype getLaplacianStream(const int index1, const int index2){

		assert(index1>= 0 && index1 < len);
		assert(index2>= 0);
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
			 		return 0.0f;
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
		return false;
	}

	void setBoundaryLength(){

		int counter = 0;
		int x, y;
		for (int i = 0; i < len; ++i)
		{
			for (int i = -1; i < 2; i+=2)
			{
				for (int j = -1; j < 2; j+=2)
				{
					x = u[i].x_index + i;
					y = u[i].y_index + j;

					if( !(x-1 < 0) && !(y-1 < 0)){
						if(isInGrid(x,y)==false){
							++counter;
						}
					}
				}
			}
		}
		lenOfBoundary = counter;
	}

	void setB(const int i, const datatype val){
		assert(i>=0 && i<len);

		b[i].stream = val;
	}


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

	void interpolateU(AdaptiveGrid fine){
		for (int i = 0; i < fine.len; ++i)
		{
			for (int j = 0; j < this->len; ++j)
			{
				if( this->u[j].x_index_global == fine.u[i].x_index_global 
					&& this->u[j].y_index_global == fine.u[i].y_index_global )
				{
					fine.u[i].stream = this->u[j].stream;
				}
			}
		}

		for (int i = 0; i < numOfPointsChosen; ++i)
		{
			Node *tmpNodePtr;
			//Interpolate to the left
			tmpNodePtr = pointsChosenByWavelet[i]->nodeLeft;
			tmpNodePtr->stream = (tmpNodePtr->nodeAbove->stream + tmpNodePtr->nodeBelow->stream )/2.0f;

			//to the right
			tmpNodePtr = pointsChosenByWavelet[i]->nodeRight;
			tmpNodePtr->stream = (tmpNodePtr->nodeAbove->stream + tmpNodePtr->nodeBelow->stream )/2.0f;

			//up
			tmpNodePtr = pointsChosenByWavelet[i]->nodeAbove;
			tmpNodePtr->stream = (tmpNodePtr->nodeLeft->stream + tmpNodePtr->nodeRight->stream )/2.0f;

			//down
			tmpNodePtr = pointsChosenByWavelet[i]->nodeBelow;
			tmpNodePtr->stream = (tmpNodePtr->nodeLeft->stream + tmpNodePtr->nodeRight->stream )/2.0f;

			//and don't forget the point itself.
			pointsChosenByWavelet[i]->stream = ( pointsChosenByWavelet[i]->nodeLeft->stream + 
				pointsChosenByWavelet[i]->nodeRight->stream )/2.0f;
		}
	}

	Node * findNode(const int ind_x, const int ind_y){

		for (int i = 0; i < len; ++i)
		{
			if(u[i].x_index == ind_x && u[i].y_index == ind_y){
				return &u[i];
			}
		}
		assert(0);
		return NULL;

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

	void setupGrid(Node* savedNodes , const int numberOfPoints){
		free(pointsChosenByWavelet);
		pointsChosenByWavelet = (Node**) malloc(numberOfPoints*sizeof(Node*));
		this->numOfPointsChosen = numberOfPoints;

		free(u);
		//9 is the maximum. It's probably waaaay too large. OPT!
		u = (Node*) malloc(9*numberOfPoints*sizeof(Node));

		int counter = numberOfPoints;

		//Take the points given by 
		for (int i = 0; i < numberOfPoints; ++i)
		{
			pointsChosenByWavelet[i] = &savedNodes[i];

			u[i] = savedNodes[i];

			//Left
			if(isInGrid(savedNodes[i].x_index-1, savedNodes[i].y_index) == false){
				*u[i].nodeLeft = u[i];
				u[counter] = *u[i].nodeLeft;
				u[counter].x_index--;
				//Make sure to do something about the global index here! FIX!
				counter++;
			}

			//Right
			if(isInGrid(savedNodes[i].x_index+1, savedNodes[i].y_index) == false){
				*u[i].nodeRight = u[i];
				u[counter] = *u[i].nodeRight;
				u[counter].x_index++;
				//Make sure to do something about the global index here! FIX!
				counter++;
			}

			//Up
			if(isInGrid(savedNodes[i].x_index, savedNodes[i].y_index+1) == false){
				*u[i].nodeAbove = u[i];
				u[counter] = *u[i].nodeAbove;
				u[counter].y_index++;
				//Make sure to do something about the global index here! FIX!
				counter++;
			}

			//down
			if(isInGrid(savedNodes[i].x_index-1, savedNodes[i].y_index-1) == false){
				*u[i].nodeBelow = u[i];
				u[counter] = *u[i].nodeBelow;
				u[counter].x_index--;
				//Make sure to do something about the global index here! FIX!
				counter++;
			}

			//Don't forget to create the corners as well.
			for (int x_mod = -1; x_mod <= 1; x_mod += 2)
			{
				for (int y_mod = -1; y_mod <= 1; x_mod += 2)
				{
					if(isInGrid(savedNodes[i].x_index+x_mod, savedNodes[i].y_index+y_mod) == false){
						u[counter] = u[i];
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
				u[i].nodeLeft = findNode(u[i].x_index-1, u[i].y_index);
			}
			if (u[i].nodeRight == NULL)
			{
				u[i].nodeRight = findNode(u[i].x_index+1, u[i].y_index);
			}
			if(u[i].nodeAbove == NULL)
			{
				u[i].nodeAbove = findNode(u[i].x_index, u[i].y_index+1);
			}
			if (u[i].nodeBelow == NULL)
			{
				u[i].nodeBelow = findNode(u[i].x_index, u[i].y_index-1);
			}
		}

		std::cout<<"fill rate: "<<(float)counter/(9.0f*numberOfPoints)<<std::endl;

	}

};

int main(int argc, char const *argv[])
{

	AdaptiveGrid grid;

	AdaptiveGrid coarse;
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

	return 0;
}
