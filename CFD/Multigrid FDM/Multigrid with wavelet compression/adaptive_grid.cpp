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

    //adaptiveGrid* finerGrid;
    //adaptiveGrid* coarserGrid;
    datatype h;



	/*
	void getPoint(const int x, const int y, bool* isInBoundary, datatype* val){
		if(isInBoundary(x,y,layer, val)){
			*isInBoundary = true;
		}
		else{
			*isInBoundary = false;

		}
	}
	*/

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

						tmpSum += this->getLaplacianStream(i,j)*(this->getU(j)); //-=

					}
				}
				setU( (this->getB(i)-tmpSum) / this->getLaplacianStream(i,i), i );
			}
		}

};


int main(int argc, char const *argv[])
{

	AdaptiveGrid grid;
	int n = 21;

	grid.u = (Node*) malloc(n*sizeof(Node));
	grid.b = (Node*) malloc(n*sizeof(Node));
	grid.h = 0.25f;

	int counter = 0;
	for (int x = 0; x < 5; ++x)
	 {
	 	for (int y = 0; y < 5; ++y)
	 	{
	 		if(y>2 && x<3){
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

	//for (int i = 0; i < 1; ++i)
	//{
		grid.jacobiSmootherLaplacianStream();
		grid.jacobiSmootherLaplacianStream();
		grid.jacobiSmootherLaplacianStream();
		grid.jacobiSmootherLaplacianStream();
	//}
	

	std::cout<<std::endl;
	for (int i = 0; i < grid.len; ++i)
	 {
	 	std::cout<<grid.u[i].x_index << " ; " <<grid.u[i].y_index << " ; " <<grid.u[i].stream << "\n" ;
	 }


	free(grid.u);
	free(grid.b);

	return 0;
}
