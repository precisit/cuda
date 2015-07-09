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
	
	std::cout<<"funkar? 1"<<std::endl;
	setMassMat(M, &phi,n);
	
	prmat(M,n*n);
	std::cout<<"funkar? 2"<<std::endl;
	setStiffMat(A, &phi,n);
	std::cout<<"funkar? 3"<<std::endl;
	prmat(A,n*n);
	
	cudaMalloc ((void**)&this->devM, n*n*n*n*sizeof(float));
	cudaMalloc ((void**)&this->devA, n*n*n*n*sizeof(float));
	cudaMalloc ((void**)&this->devU, n*n*sizeof(float));
	cudaMalloc ((void**)&this->devB, n*n*sizeof(float));
	
	cublasSetMatrix(n*n,n*n, sizeof(*A), A, n*n, this->devA, n*n);
	cublasSetMatrix(n*n,n*n, sizeof(*M), M, n*n, this->devM, n*n);
	cublasSetVector(n*n, sizeof(float), this->U, 1, this->devU, 1);
	
	free(A); free(M);
	
	float k = - dt;
	const float one = 1.0f;
	
	//Reassign A as A := M+k*A; No. Don't do it. Just don't.
	//cublasSgeam(this->handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, &k, this->devA, n, &one, this->devM, n, this->devA, n);
};


void MrWags::waveletGalerkinIter(){

	for(int i= 0; i< this->n*this->n; i++){
		assert(this->U[i] >= 0.0f);
		assert(this->U[i] <= 1000.0f);
	}

	std::cout<<"iter"<<std::endl;

	const int n = this->n;
	const float zero = 0.0f;
	const float one = 1.0f;
	this->t += this->dt;
	float k = -this->dt;
	
	//b = (M+k*A)U_old
	//Start with M := (M+k*A)
	cublasSaxpy(this->handle,n*n*n*n, &k, this->devA, 1, this->devM, 1);
	
	/*float * tmp;
	tmp = (float*) calloc (n*n*n*n,sizeof(*tmp));
	cublasGetVector(n*n*n*n, sizeof(float), devM, 1, tmp,1); //Det här är segt. FIX!
	for(int i=0; i<n*n*n*n; i++){
		std::cout<<"b (but it's really M): "<<tmp[i]<<std::endl;
	}
	free(tmp);
	*/
	
	//Then b = M*U_old
	cublasSgemv(this->handle, CUBLAS_OP_N, n*n,n*n, &one, this->devM, n*n, this->devU, 1, &zero, this->devB,1);
	
	
	//
	//b = A*U_old; //HÄR ÄR FELET. FIX!
	//	cublasSgemv(this->handle, CUBLAS_OP_N, n*n, n*n, &one, this->devA, n*n, this->devU, 1, &zero, devB, 1); //M är förmodligen bandmarix. Använd annan funk!! //FIX
	
	//Solve M*u = b
	matSolve(this->devM, this->devU, this->devB, n*n, this->tol ,this->handle);
		
	cublasGetVector(n*n, sizeof(float), this->devU, 1, this->U,1); //Det här är segt. FIX!

	for(int i= 0; i< this->n*this->n; i++){
		std::cout<<U[i]<<std::endl;
		assert(this->U[i] >= 0.0f);
		assert(this->U[i] <= 1.0f);
	}


	assert(this->t < this->endTime);

};


void MrWags::waveletGalerkinEnd(){
	cudaFree(this->devM);
	cudaFree(this->devA);
	cudaFree(this->devU);
	cudaFree(this->devB);
	
	cublasDestroy(this->handle);
}

int MrWags::getN(){
	return this->n;
}

float MrWags::getU(int x, int y){
	return this->U[x*(this->n)+y];
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
	float dx = (x1-x0)/((float) N);
	float dy = (y1-y0)/((float) N);
	
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


	
	float * tmp;
	tmp = (float*) calloc (n*n,sizeof(*tmp));
	cublasGetMatrix(n,n, sizeof(float), A, n, tmp,n); //Det här är segt. FIX!
	for(int i=0; i<n*n; i++){
		std::cout<<"A: "<<tmp[i]<<std::endl;
	}
	free(tmp);
	
	tmp = (float*) calloc (n,sizeof(*tmp));
	cublasGetVector(n, sizeof(float), x, 1, tmp,1); //Det här är segt. FIX!
	for(int i=0; i<n; i++){
		std::cout<<"x: "<<tmp[i]<<std::endl;
	}
	free(tmp);
	
	tmp = (float*) calloc (n,sizeof(*tmp));
	cublasGetVector(n, sizeof(float), b, 1, tmp,1); //Det här är segt. FIX!
	for(int i=0; i<n; i++){
		std::cout<<"b: "<<tmp[i]<<std::endl;
	}
	free(tmp);
	



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
	
	//r := r -A*x
	status = cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_one, A,n,x,1,&one,r, 1);
	
	
	std::cout<<status<<std::endl;
	assert(status == CUBLAS_STATUS_SUCCESS);
	if(status != CUBLAS_STATUS_SUCCESS){
		assert(0);
	}
	else{
		std::cout<<"sdgdfgsdfg"<<std::endl;
	}
	
	
	//float * tmp;
	tmp = (float*) malloc (n*sizeof(*tmp));
	tmp[0] = .1748792347875;
	cublasGetVector(n, sizeof(float), r, 1, tmp,1); //Det här är segt. FIX!
	for(int i=0; i<n; i++){
		std::cout<<"r: "<<tmp[i]<<std::endl;
	}
	free(tmp);
	
	
	//define r_norm_old
	cublasSdot(handle,n,r,1,r,1,&r_norm_old);
	std::cout<<"r_norm_old "<<r_norm_old<<std::endl;
	//cublasSnrm2(handle,n,r,1, &r_norm_old);
	
	//p = r
	cublasScopy(handle, n, r, 1, p,1);
	
	for(int i=0; i<n*n; i++){
		std::cout<<"Iter. nr: "<<i<<std::endl<<std::endl;
		
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
		std::cout<<"alpha: "<<alpha<<std::endl;
		assert(alpha >= 0.0f);
		
		//x = x + alpha*p
		cublasSaxpy(handle,n, &alpha, p, 1, x, 1);
		
		//r = r - alpha*A*p (=r-alpha*Ap. FIX! cublasSaxpy() borde funka)
		negative_alpha = - alpha;
		cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_alpha, A,n,p,1,&one,r, 1);
	
		//if r small: break
		cublasSdot(handle,n,r,1,r,1,&r_norm);
		std::cout<<"r_norm: "<<r_norm<<std::endl;
		
		
		//update beta
		beta = r_norm/r_norm_old;
		std::cout<<"beta: "<<beta<<std::endl;
		
		
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
	std::cout<<"matSolve END"<<std::endl;
	cudaFree(r);
	cudaFree(p);
	cudaFree(Ap);
}

template<typename Element>
void MrWags::setStiffMat(float *A, const Element *phi, const int n){
	const float h = phi->h; 
	
	for(int x1=0; x1<n; x1++){
		for(int y1=0; y1<n; y1++){
			for(int x2=0; x2<n; x2++){
				for(int y2=0; y2<n; y2++){
					Element phi1 = Element(((float)x1)/((float)(n-1)), ((float)y1)/((float)(n-1)), h);
					Element phi2 = Element(((float)x2)/((float)(n-1)), ((float)y2)/((float)(n-1)), h);
					A[(x1*n+y1)*n*n+(x2*n+y2)] = numIntergStiff(phi1, phi2 , 0-h, 1+h, 0-h, 1+h, n*200);
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
	//int d;
	
	for(int x1=0; x1<n; x1++){
		for(int y1=0; y1<n; y1++){
			for(int x2=0; x2<n; x2++){
				for(int y2=0; y2<n; y2++){
						Element phi1 = Element(((float)x1)/((float)(n-1)), ((float)y1)/((float)(n-1)), h);
						Element phi2 = Element(((float)x2)/((float)(n-1)), ((float)y2)/((float)(n-1)), h);
						M[(x1*n+y1)*n*n+(x2*n+y2)] = numIntergMass(phi1, phi2 , 0-h, 1+h, 0-h, 1+h, n*200);
					
				}
			}
		}
	}
}
