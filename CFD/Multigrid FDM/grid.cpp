#include "grid.h"
	/*
		Contains a structured uniform 2D grid of size len-by-len.

		(Each grid object needs to copy the data (the Node objects)
		instead of having pointers to the pointers or some thing 
		clever. This should probably be optimized. (OPT!))
	*/

		void Grid::print(){
			int n = this-> len;

			for(int x=0; x<n; x++){
				for(int y=0; y<n; y++){
					if(this->ptr[x+y*n].p == 0.0f){
						printf("0         ");
					}
					else if(this->ptr[x+y*n].p>0.0f){
				 		printf ("%3f  ", this->ptr[x+y*n].p);
						//std::cout<<A[x+y*n]<<"   ";
					}
					else{
						printf ("%3f ", this->ptr[x+y*n].p);
					}
				}
				std::cout<<std::endl;
			}
		}

		void Grid::restrict(Grid* coarseGrid){

			int counter = 0;
			const int n = this->len;

			for( int iter1= 0; iter1<n; iter1 += 2){
				for( int iter2 = 0; iter2<n; iter2++ ){
					if( iter2 % 2 == 0){
						coarseGrid->ptr[counter]  =  this->ptr[iter1*n + iter2];
						counter++;
					}
				}
			}
		};

		void Grid::interpolate(Grid* fineGrid, const int fineN){

			const int n = this->len;
			//int n = fineN;

			//Copy the nodal values
			for( int x = 0; x< n; x++){
	    		for( int y = 0; y< n; y++){
	        		fineGrid->ptr[(2*x)*fineN+2*y] = this->ptr[x*n+y];
	    		}
	    	}

	    	fineGrid->print();
	    	std::cout<<std::endl;

	    	//Interpolate in the x-dir (might be y-dir. It's possible I fucked up here. But it doesn't really matter.)
	    	for( int x = 0; x < n; x++){
	    		for( int y = 1; y < n; y++){
	        		fineGrid->ptr[2*x*fineN+2*y-1] = (this->ptr[x*n+y]+this->ptr[x*n+y-1])/2.0f;
	        	}
	        }

	        fineGrid->print();
	    	std::cout<<std::endl;

	    	//Interpolate in the y-dir
	    	for( int x = 1; x< fineN; x+=2){
	    		for( int y = 0; y< fineN; y++){
	        		fineGrid->ptr[x*fineN+y] = (fineGrid->ptr[(x-1)*fineN+y]+fineGrid->ptr[(x+1)*fineN+y])/2.0f;
	        	}
	        }

		};

		void Grid::setValues(Node* ptrIn, int n){
			if(this->len != 0){
				assert(n == this->len);
			}

			for( int i=0; i<n*n; i++){
				this->ptr[i] = ptrIn[i];
			}
		};

		datatype Grid::getU(const int index){
			assert(index < 3*this->len*this->len);
			assert(index >= 0);

			if( index < this->len*this->len ){
				return this->ptr[index].p;
			}
			else if( index < 2*this-> len*this->len ){
				return ptr[index-1*this-> len*this->len].v_x;
			}
			else{
				return ptr[index-2*this-> len*this->len].v_y;
			}
		}

		void Grid::setU(const datatype u, const int index){
			assert(index < 3*this->len*this->len);
			assert(index >= 0);

			if( index < this->len*this->len ){
				ptr[index].p = u;
			}
			else if( index < 2*this->len*this->len ){
				ptr[index].v_x = u;
			}
			else{
				ptr[index].v_y = u;
			}
		}

		datatype Grid::getLaplacian(const int i, const int j){
			int dist = abs(i-j);
			if( dist == 0){
				return -4.0f/(this->h*this->h);
			}
			else if( dist == 1 ){
				return 1.0f/(this->h*this->h);
			}
			else if( dist == this->len ){
				return 1.0f/(this->h*this->h);
			}
			else{
				return 0.0f;
			}
		}

		datatype Grid::getVortTranspDisc(const int i, const int j){
			/*
				This is the discreatesation (that's not how that's
				spelled. At all.) of the Vort. Transport differential
				operator.
			*/
			//i is the x index, and j is the y.
			int dist = abs(i-j);
			if( dist == 0){
				return -4.0f/(this->h*this->h*Re);
			}
			else if( dist == this->len ){
				if( i > j){
					return 1.0f/this->h*(-this->getU(1*this->len*this->len + j )/2.0f+1.0f/(Re*this->h));
				}
				else{
					return 1.0f/this->h*(this->getU(1*this->len*this->len + j )/2.0f+1.0f/(Re*this->h));
				}
			}
			else if( dist == 1 ){
				if( i > j){
					return 1.0f/this->h*(-this->getU(2*this->len*this->len + j )/2.0f+1.0f/(Re*this->h));
				}
				else{
					return 1.0f/this->h*(this->getU(2*this->len*this->len + j )/2.0f+1.0f/(Re*this->h));
				}
			}
			else{
				return 0.0f;
			}

		}

		datatype Grid::getB(const int i){
			return 0.0f;
		} //FIX!

		int Grid::lengthOfFinerGrid(){
			return (2*(this->len-1)+1);
		}

		int Grid::lengthOfCoarserGrid(){
			return ((this->len-1)/2+1);
		}

		void Grid::jacobiSmootherLaplacian(){
			/*
				Smoothes the error to D u = b.
				u = (p, v_x, v_y)^T
			*/

			datatype tmpSum;
			for(int i=0; i<this->len*this->len; i++){
				tmpSum = 0.0f;
				for(int j=0; j<this->len*this->len; j++){
					if(j != i){

						tmpSum += -this->getLaplacian(i,j)*(this->getU(i));

					}
				}
				setU( (this->getB(i)-tmpSum) / this->getLaplacian(i,i), i );
			}
		}

		void Grid::jacobiSmootherVortTranspDisc(){
			/*
				Smoothes the error to D u = b.
				u = (p, v_x, v_y)^T
			*/

			datatype tmpSum;
			for(int i=0; i<this->len*this->len; i++){
				tmpSum = 0.0f;
				for(int j=0; j<this->len*this->len; j++){
					if(j != i){

						tmpSum += -this->getVortTranspDisc(i,j)*(this->getU(i));

					}
				}
				setU( (this->getB(i)-tmpSum) / this->getVortTranspDisc(i,i), i );
			}
		}

		//Constructors
		Grid::Grid(Node* nodeIn, const int nIn){
			this->ptr = nodeIn;
			this->len = nIn;

			assert(nIn > 1);
			this->h = 1.0f / ((float) nIn - 1.0f);

			this->ptr = (Node*) malloc(nIn*nIn*sizeof(Node));

			this->setValues(nodeIn, nIn);
		};

		Grid::Grid(const int n){
			this->len = n;
			this->ptr = (Node*) malloc(n*n*sizeof(Node));
			assert(n > 1);

			this->h = 1.0f / ((float) n - 1.0f);
		};

		Grid::Grid(){
			this->len = 0;
		};


		void Grid::calculateError(){
			/*

			*/
			
		};
