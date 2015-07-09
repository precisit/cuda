#include "mrWags.h"

MrWags::MrWags(float* U_0, const int n, const float dt, const float endTime, const float tol){
	this->n = n;
	this->dt = dt;
	this->endTime = endTime;
	this->tol = tol;
	
	this->t = 0.0f;
	
	cublasCreate(&(this->handle));
	
	this->U = U_0;
	
	this->waveletGalerkinInit();	

	//this->waveletGalerkinEnd();
};

/*
// This function is called every time GLUT refreshes the display.
void MrWags::draw(void)
{
	this->t += this->dt;

 	this->waveletGalerkinIter();
 	this->drawMatrix();

	assert(this->t < this->endTime);
}
*/

void MrWags::waveletGalerkinInit(){
	
	assert(this->endTime>0.0f);
	assert(this->n>0);
	
	const int n = this->n;
	
	//Setup the Element (should be a wavelet, but isn't yet. FIX!)
	Element phi = Element(0.0f, 0.0f, 1.0f/((float)(n-1)));

	//Setup mass and stiffness matrix;
	float *M,*A;
	M = (float*) calloc (n*n*n*n,sizeof(*M));
	A = (float*) calloc (n*n*n*n,sizeof(*A));
	
	//std::cout<<"funkar? 1"<<std::endl;
	setMassMat(M, &phi,n);
	
	//prmat(M,n*n);
	//std::cout<<"funkar? 2"<<std::endl;
	setStiffMat(A, &phi,n);
	std::cout<<"funkar? 3"<<std::endl<<std::endl;
	//prmat(A,n*n);
	//assert(0);
	cudaMalloc ((void**)&this->devM, n*n*n*n*sizeof(float));
	cudaMalloc ((void**)&this->devA, n*n*n*n*sizeof(float));
	cudaMalloc ((void**)&this->devU, n*n*sizeof(float));
	cudaMalloc ((void**)&this->devB, n*n*sizeof(float));

	//The devStorage matrix is used to store matrices that won't be used 
	//multiple times.
	cudaMalloc ((void**)&this->devStorage, n*n*n*n*sizeof(float));
	
	cublasSetMatrix(n*n,n*n, sizeof(*A), A, n*n, this->devA, n*n);
	cublasSetMatrix(n*n,n*n, sizeof(*M), M, n*n, this->devM, n*n);
	cublasSetVector(n*n, sizeof(float), this->U, 1, this->devU, 1);
	
	free(A); free(M);
};


void MrWags::prmatgpu(float *devA, int n){

	float *tmp;
	tmp = (float*) calloc (n*n,sizeof(*tmp));
	cublasGetMatrix(n,n, sizeof(float), devA, n, tmp,n);
	prmat(tmp, n);
	free(tmp);
}

void MrWags::prvecgpu(float *devA, int n){

	float * tmp;
	tmp = (float*) calloc (n,sizeof(*tmp));
	cublasGetVector(n, sizeof(float), devA, 1, tmp,1);
	for(int i=0; i<n; i++){
		std::cout<<tmp[i]<<std::endl;
	}
	free(tmp);
}

void MrWags::waveletGalerkinIter(){

	/*
	for(int i= 0; i< this->n*this->n; i++){
		assert(this->U[i] >= 0.0f);
		assert(this->U[i] <= 1000.0f);
	}
	*/

	std::cout<<"iter"<<std::endl;

	const int n = this->n;
	const float zero = 0.0f;
	const float one = 1.0f;
	this->t += this->dt;
	float k = -this->dt;
	
	/*
	//THIS IS FOR FORWARD EULER (USE BACKWARD INSTEAD!)
	//Define b := (M+k*A)U_old
	//Start by S := M
	cublasScopy(this->handle,n*n*n*n,this->devM, 1, this->devStorage, 1);
	//Then S := S+k*A. Which is the same as M-dt*A 
	cublasSaxpy(this->handle,n*n*n*n, &k, this->devA, 1, this->devStorage, 1);
	//And finally b := S*U_old, thus giving b := (M+k*A)U_old
	cublasSgemv(this->handle, CUBLAS_OP_N, n*n,n*n, &one, this->devStorage, n*n, this->devU, 1, &zero, this->devB,1);
	//Solve M*u = b
	matSolve(this->devM, this->devU, this->devB, n*n, this->tol ,this->handle);
	*/

	//THIS IS FOR BACKWARD EULER.
	//Define b := M*U_old
	cublasSgemv(this->handle, CUBLAS_OP_N, n*n,n*n, &one, this->devM, n*n, this->devU, 1, &zero, this->devB,1);

	//Start by S := M
	cublasScopy(this->handle,n*n*n*n,this->devM, 1, this->devStorage, 1);
	//k = dt;
	k = this->dt;
	//Then S := S+k*A. Which is the same as M+dt*A 
	cublasSaxpy(this->handle,n*n*n*n, &k, this->devA, 1, this->devStorage, 1);
	//Solve (M+dt*A)*u = M*u_old
	matSolve(this->devStorage, this->devU, this->devB, n*n, this->tol ,this->handle);
		
	cublasGetVector(n*n, sizeof(float), this->devU, 1, this->U,1); //Det här är segt. FIX!

	for(int i= 0; i< this->n*this->n; i++){
		//std::cout<<U[i]<<std::endl;
		//assert(this->U[i] >= 0.0f);
		assert(this->U[i] <= 1.0f);
	}


	assert(this->t < this->endTime);

};


void MrWags::waveletGalerkinEnd(){
	cudaFree(this->devM);
	cudaFree(this->devA);
	cudaFree(this->devU);
	cudaFree(this->devB);

	cudaFree(this->devStorage);
	
	cublasDestroy(this->handle);
}

int MrWags::getN(){
	return this->n;
}

float MrWags::getU(const int x,const int y){
	return this->U[x*(this->n)+y];
}

int MrWags::dist_l1(const int x1, const int y1, const int x2, const int y2){
	//Returns the distance in the l^1-norm.
	return std::abs(x1-x2)+std::abs(y1-y2);
}
















void MrWags::wait(int sec){
	time_t now, timer;
	time(&now);
	time(&timer);
	int counter;
	while(std::abs(difftime(timer, now))<sec){
		counter++;
		time(&timer);
	}
	return;
}
	

void MrWags::prvec(float *A, int n){
	for(int x=0; x<n; x++){
		std::cout<<A[x]<<"   ";
	}
}

void MrWags::prmat(float *A, int n){
	for(int x=0; x<n; x++){
		for(int y=0; y<n; y++){
			if(A[x+y*n] == 0.0f){
				printf("0         ");
			}
			else if(A[x+y*n]>0.0f){
		 		printf ("%3f  ", A[x+y*n]);
				//std::cout<<A[x+y*n]<<"   ";
			}
			else{
				printf ("%3f ", A[x+y*n]);
			}
		}
		std::cout<<std::endl;
	}
}

template<typename T>
float MrWags::numIntergStiff(T func1, T func2 , float x0, float x1, float y0, float y1, int N){
	assert(x0<x1);
	assert(y0<y1);
	assert(N>0);

	float sum = 0.0f;
	float dx = (x1-x0)/((float) N);
	float dy = (y1-y0)/((float) N);
	
	assert(dx>0);
	
	float gridSize = dx*dy;
	
	float x = x0;
	float y = y0;
	
	while(x<=x1){
		
		y = y0;
		while(y<=y1){
			sum += (func1.f_x(x,y,dx/100.0f)+func1.f_y(x,y,dy/100.0f)) * (func2.f_x(x,y,dx/100.0f)+func2.f_y(x,y,dy/100.0f));
			y += dy;
		}
		x += dx;
	}
	sum = sum * gridSize;
	return sum;
}

template<typename T>
float MrWags::numIntergMass(T func1, T func2 , float x0, float x1, float y0, float y1, int N){
	assert(x0<x1);
	assert(y0<y1);
	assert(N>0);

	float sum = 0.0f;
	float dx = (x1-x0)/((float) N-1);
	float dy = (y1-y0)/((float) N-1);
	
	assert(dx>0);
	
	float gridSize = dx*dy;
	
	float x = x0;
	float y = y0;
	
	while(x<=x1){
		
		y = y0;
		while(y<=y1){
			sum += func1.f(x,y) * func2.f(x,y);
			y += dy;
		}
		x += dx;
	}
	sum = sum * gridSize;
	return sum;
}

//Solves Ax=b for x. (A is n-by-n and also symmetric).
//Based on CG.
void MrWags::matSolve(const float *A, float *x, float *b, const int n, const float absTol, cublasHandle_t handle){
	//Stolen from here: https://en.wikipedia.org/wiki/Conjugate_gradient_method


	//float *tmp;
	/*
	tmp = (float*) calloc (n*n,sizeof(*tmp));
	cublasGetMatrix(n,n, sizeof(float), A, n, tmp,n);
	for(int i=0; i<n*n; i++){
		std::cout<<"A: "<<tmp[i]<<std::endl;
	}
	free(tmp);
	
	tmp = (float*) calloc (n,sizeof(*tmp));
	cublasGetVector(n, sizeof(float), x, 1, tmp,1);
	for(int i=0; i<n; i++){
		std::cout<<"x: "<<tmp[i]<<std::endl;
	}
	free(tmp);
	
	tmp = (float*) calloc (n,sizeof(*tmp));
	cublasGetVector(n, sizeof(float), b, 1, tmp,1);
	for(int i=0; i<n; i++){
		std::cout<<"b: "<<tmp[i]<<std::endl;
	}
	free(tmp);
	*/
	
	float *r, *p, *Ap;
	float r_norm_old, r_norm;
	
	cudaMalloc ((void**)&r, n*sizeof(*r));
	cudaMalloc ((void**)&p, n*sizeof(*p));
	cudaMalloc ((void**)&Ap, n*sizeof(*Ap));
	
	//Kolla om man får nån speed-up av att flytta dessa till GPUn. FIX!
	float alpha, beta, negative_alpha;
	
	//Kolla om man får nån speed-up av att flytta dessa till GPUn. FIX!
	float one = 1.0f;
	float negative_one = -1.0f;
	float zero = 0.0f;
	
	
	cublasStatus_t status;
	//r = b-A*x
	cublasScopy(handle,n, b, 1, r,1);
	
	/*
	tmp = (float*) malloc (n*sizeof(*tmp));
	cublasGetVector(n, sizeof(float), r, 1, tmp,1); //Det här är segt. FIX!
	for(int i=0; i<n; i++){
		std::cout<<"r (copy): "<<tmp[i]<<std::endl;
	}
	free(tmp);
	*/

	/*
	//ALLT UNDER ÄR FÖR DEBUGGING!
	float *Anew, *Bnew, *Xnew;
	Anew = (float*) calloc (n*n,sizeof(*Anew));
	Bnew = (float*) calloc (n,sizeof(*Bnew));
	Xnew = (float*) calloc (n,sizeof(*Xnew));
	
	float *devAnew, *devXnew, *devBnew;
	cudaMalloc ((void**)&devAnew, n*n*sizeof(*devAnew));
	cudaMalloc ((void**)&devBnew, n*sizeof(*devBnew));
	cudaMalloc ((void**)&devXnew, n*sizeof(*devXnew));
	
	//Anew[0] = 1;
	//Anew[5] = 1;
	//Anew[10] = 1;
	//Anew[15]= 1;
	
	Anew[0]= 0.666246;
	Anew[1]= 0.208336;
	Anew[2]= 0.208339;
	Anew[3]= 0.0834405;
	Anew[4]= 0.208336;
	Anew[5]= 0.666246;
	Anew[6]= 0.0832422;
	Anew[7]= 0.208338;
	Anew[8]= 0.208339;
	Anew[9]= 0.0832422;
	Anew[10]= 0.666318;
	Anew[11]= 0.20834;
	Anew[12]= 0.0834405;
	Anew[13]= 0.208338;
	Anew[14]= 0.20834;
	Anew[15]= 0.666315;
	
	
	//Xnew[0] = 0.5;
	//Xnew[2] = 1;
	
	Xnew[0] = 0.01;
	Xnew[1] = 0.01;
	Xnew[2] = 0.01;
	Xnew[3] = 1000;
	
	
	//Bnew[0] = 3;

	//Bnew[0] = 83.4513;
	//Bnew[1] = 208.348;
	//Bnew[2] = 208.35;
	//Bnew[3] = 666.32;
	

	status = cublasSetMatrix(n,n, sizeof(float), Anew, n, devAnew,n);
	if(status != CUBLAS_STATUS_SUCCESS){
		//printf(_cudaGetErrorEnum(status));
		std::cout<<"1"<<std::endl;
		assert(0);
		//return EXIT_FAILURE;
	}
	status = cublasSetVector(n, sizeof(float), Bnew, 1, devBnew,1);
	if(status != CUBLAS_STATUS_SUCCESS){
		//printf(_cudaGetErrorEnum(status));
		std::cout<<"2"<<std::endl;
		assert(0);
		//return EXIT_FAILURE;
	}
	status = cublasSetVector(n, sizeof(float), Xnew, 1, devXnew,1);
	if(status != CUBLAS_STATUS_SUCCESS){
		//printf(_cudaGetErrorEnum(status));
		std::cout<<"3"<<std::endl;
		assert(0);
		//return EXIT_FAILURE;
	}
	
	
	
	
	status = cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_one, A,n,x,1,&zero,devBnew, 1);
	
	tmp = (float*) malloc (n*sizeof(*tmp));
	cublasGetVector(n, sizeof(float),devBnew , 1, tmp,1); //Det här är segt. FIX!
	for(int i=0; i<n; i++){
		std::cout<<"Bnew (1): "<<tmp[i]<<std::endl;
	}
	free(tmp);
	
	cublasSaxpy(handle,n, &one, b, 1, devBnew, 1);
	
	
	tmp = (float*) malloc (n*sizeof(*tmp));
	cublasGetVector(n, sizeof(float),devBnew , 1, tmp,1); //Det här är segt. FIX!
	for(int i=0; i<n; i++){
		std::cout<<"Bnew (2): "<<tmp[i]<<std::endl;
	}
	free(tmp);
	//DEBUGGING SLUTAR HÄR!
	
	*/
	
	
	//r := r -A*x
	status = cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_one, A,n,x,1,&one,r, 1);
	
	
	//std::cout<<status<<std::endl;
	assert(status == CUBLAS_STATUS_SUCCESS);
	if(status != CUBLAS_STATUS_SUCCESS){
		assert(0);
	}
	
	
	/*
	tmp = (float*) malloc (n*sizeof(*tmp));
	tmp[0] = 0.1748792347875f;
	cublasGetVector(n, sizeof(float), r, 1, tmp,1); //Det här är segt. FIX!
	for(int i=0; i<n; i++){
		std::cout<<"r: "<<tmp[i]<<std::endl;
	}
	free(tmp);
	*/
	
	
	//define r_norm_old
	cublasSdot(handle,n,r,1,r,1,&r_norm_old);
	//std::cout<<"r_norm_old "<<r_norm_old<<std::endl;
	//cublasSnrm2(handle,n,r,1, &r_norm_old);
	
	//p = r
	cublasScopy(handle, n, r, 1, p,1);
	
	for(int i=0; i<n*n; i++){
		//std::cout<<"Iter. nr: "<<i<<std::endl<<std::endl;
		
		/*
		float * tmp;
		tmp = (float*) calloc (n,sizeof(*tmp));
		cublasGetVector(n, sizeof(float), x, 1, tmp,1); //Det här är segt. FIX!
		for(int i=0; i<n; i++){
			std::cout<<"x2: "<<tmp[i]<<std::endl;
		}
		free(tmp);
		*/
		//Ap = A*p
		cublasSgemv(handle,CUBLAS_OP_N,n,n,&one, A,n,p,1,&zero,Ap, 1);
		
		//update alpha
		cublasSdot(handle,n,p,1,Ap,1,&alpha);
		alpha = r_norm_old/alpha;
		//std::cout<<"alpha: "<<alpha<<std::endl;
		assert(alpha >= 0.0f);
		
		//x = x + alpha*p
		cublasSaxpy(handle,n, &alpha, p, 1, x, 1);
		
		//r = r - alpha*A*p (=r-alpha*Ap. FIX! cublasSaxpy() borde funka)
		negative_alpha = - alpha;
		cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_alpha, A,n,p,1,&one,r, 1);
	
		//if r small: break
		cublasSdot(handle,n,r,1,r,1,&r_norm);
		//std::cout<<"r_norm: "<<r_norm<<std::endl;
		
		
		//update beta
		beta = r_norm/r_norm_old;
		//std::cout<<"beta: "<<beta<<std::endl;
		
		
		if(r_norm < absTol){
		
			/*
			float * tmp;
			tmp = (float*) calloc (n,sizeof(*tmp));
			cublasGetVector(n, sizeof(float), x, 1, tmp,1); //Det här är segt. FIX!
			for(int i=0; i<n; i++){
				std::cout<<"x: "<<tmp[i]<<std::endl;
			}
			free(tmp);
			*/
			break;
		}
		else{
			assert(r_norm != 0.0f);
		}
		
		/*
		//update beta
		beta = r_norm/r_norm_old;
		std::cout<<"beta: "<<beta<<std::endl;
		*/
		
		//p = r + beta*p
		 cublasSscal(handle, n,&beta, p,1);
		 cublasSaxpy(handle, n, &one, r, 1, p, 1);
		 
		 r_norm_old = r_norm;
	}
	//std::cout<<"matSolve END"<<std::endl;
	cudaFree(r);
	cudaFree(p);
	cudaFree(Ap);
}

template<typename Element>
void MrWags::setStiffMat(float *A, const Element *phi, const int n){
	//const float h = phi->h; 
	int d;
	
	for(int x1=0; x1<n; x1++){
		for(int y1=0; y1<n; y1++){
			for(int x2=0; x2<n; x2++){
				for(int y2=0; y2<n; y2++){
					//Element phi1 = Element(((float)x1)/((float)(n-1)), ((float)y1)/((float)(n-1)), h);
					//Element phi2 = Element(((float)x2)/((float)(n-1)), ((float)y2)/((float)(n-1)), h);
					//A[(x1*n+y1)*n*n+(x2*n+y2)] = numIntergStiff(phi1, phi2 , 0-h, 1+h, 0-h, 1+h, n*200);
					d = dist_l1(x1,y1,x2,y2);
					if(d< 3){
						if(d == 0){
							A[(x1*n+y1)*n*n+(x2*n+y2)] = 4.0f;
						}
						else if(d == 2){
							if( y1 != y2 && x1 != x2){
								if( (y1-y2)*(x1-x2) == 1){
									A[(x1*n+y1)*n*n+(x2*n+y2)] = -1.0f; //minus här
								}
								else{
									A[(x1*n+y1)*n*n+(x2*n+y2)] = 1.0f;
								}
							}
						}
					}
				}
			}
		}
	}
}

template<typename Element>
void MrWags::setMassMat(float *M, const Element *phi, const int n){
	//This is stupid and stuff. And shouldn't be done on the CPU. Oh, well.
	//But that'll be fixed later. 
	//(And also it's done ad hoc for hat funcs)
	
	
	//Antar att alla har samma storlek. FIX!
	const float h = phi->h;
	float val;
	
	for(int x1=0; x1<n; x1++){
		for(int y1=0; y1<n; y1++){
			for(int x2=0; x2<n; x2++){
				for(int y2=0; y2<n; y2++){
					if(dist_l1(x1,y1,x2,y2) < 3){
						Element phi1 = Element(((float)x1)/((float)(n-1)), ((float)y1)/((float)(n-1)), h);
						Element phi2 = Element(((float)x2)/((float)(n-1)), ((float)y2)/((float)(n-1)), h);
						/*
						if(x2<x1){
							if(y2<y1){
								val = numIntergMass(phi1, phi2 , x2-4.0f*h, x1+4.0f*h, y2-4.0f*h, y1+4.0f*h, n*200);
								val = 1;
							}
							else{
								val = numIntergMass(phi1, phi2 , x2-4.0f*h, x1+4.0f*h, y1-4.0f*h, y2+4.0f*h, n*200);
								val = 2;
							}
						}
						else{
							if(y2<y1){
								val = numIntergMass(phi1, phi2 , x1-4.0f*h, x2+4.0f*h, y2-4.0f*h, y1+4.0f*h, n*200);
								val = 3;
							}
							else{
								val = numIntergMass(phi1, phi2 , x1-4.0f*h, x2+4.0f*h, y1-4.0f*h, y2+4.0f*h, n*200);
								
							}
						}
						*/
						val = numIntergMass(phi1, phi2 , 0, 1.0f, 0.0f, 1.0f, n*10);
						M[(x1*n+y1)*n*n+(x2*n+y2)] = val;
					}
					
				}
			}
		}
	}
}
