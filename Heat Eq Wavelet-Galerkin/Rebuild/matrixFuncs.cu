#include "matrixFuncs.h"

void wait(int sec){
	time_t now, timer;
	time(&now);
	time(&timer);
	int counter;
	while(abs(difftime(timer, now))<sec){
		counter++;
		time(&timer);
	}
	return;
}
	

void prvec(float *A, int n){
	for(int x=0; x<n; x++){
		std::cout<<A[x]<<"   ";
	}
}

void prmat(float *A, int n){
	for(int x=0; x<n; x++){
		for(int y=0; y<n; y++){
		 	printf ("%3f  ", A[x+y*n]);
			//std::cout<<A[x+y*n]<<"   ";
		}
		std::cout<<std::endl;
	}
}

template<typename T>
float numIntergStiff(T func1, T func2 , float x0, float x1, float y0, float y1, int N){
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
			sum += (func1.f_x(x,y,dx/10.0f)+func1.f_y(x,y,dy/10.0f)) * (func2.f_x(x,y,dx/10.0f)+func2.f_y(x,y,dy/10.0f));
			y += dy;
		}
		x += dx;
	}
	sum = sum * gridSize;
	return sum;
}

template<typename T>
float numIntergMass(T func1, T func2 , float x0, float x1, float y0, float y1, int N){
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
void matSolve(const float *A, float *x, const float *b, const int n, const float absTol, cublasHandle_t handle){
	//Stolen from here: https://en.wikipedia.org/wiki/Conjugate_gradient_method
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
	
	//r = b-A*x
	cublasScopy(handle,n, b, 1, r,1);
	cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_one, A,n,x,1,&one,r, 1);
	
	//define r_norm_old
	cublasSdot(handle,n,r,1,r,1,&r_norm_old);
	//cublasSnrm2(handle,n,r,1, &r_norm_old);
	
	//p = r
	cublasScopy(handle, n, r, 1, p,1);
	
	for(int i=0; i<n; i++){
		std::cout<<"hej?"<<std::endl;
	
		//Ap = A*p
		cublasSgemv(handle,CUBLAS_OP_N,n,n,&one, A,n,p,1,&zero,Ap, 1);
		
		//update alpha
		cublasSdot(handle,n,p,1,Ap,1,&alpha);
		alpha = r_norm_old/alpha;
		assert(alpha>=0.0f);
		
		//x = x + alpha*p
		cublasSaxpy(handle,n,&alpha, p,1, x, 1);
		
		//r = r - alpha*A*p (=r-alpha*Ap. FIX! cublasSaxpy() borde funka)
		negative_alpha = - alpha;
		cublasSgemv(handle,CUBLAS_OP_N,n,n,&negative_alpha, A,n,p,1,&one,r, 1);
	
		//if r small: break
		cublasSdot(handle,n,r,1,r,1,&r_norm);
		if(r_norm < absTol){
			break;
		}
		else{
			assert(r_norm != 0.f);
		}
		
		//update beta
		beta = r_norm/r_norm_old;
		
		//p = r + beta*p
		 cublasSscal(handle, n,&beta, p,1);
		 cublasSaxpy(handle, n, &one, r, 1, p, 1);
		 
		 r_norm_old = r_norm;
	}
		
	cudaFree(r);
	cudaFree(p);
	cudaFree(Ap);
}

int dist(const int x1, const int y1, const int x2, const int y2){
	//Returns the distance (measured in the chess king's norm) between 2 points.
	return abs(x1-x2)+abs(y1-y2);
}

template<typename T>
void setStiffMat(float *A, const T *phi, const int n){
	const float h = phi->h; 
	int d;

	for(int x1=0; x1<n; x1++){
		for(int y1=0; y1<n; y1++){
			for(int x2=0; x2<n; x2++){
				for(int y2=0; y2<n; y2++){
					Element phi1 = Element(((float)x1)/((float)(n-1)), ((float)y1)/((float)(n-1)), h);
					Element phi2 = Element(((float)x2)/((float)(n-1)), ((float)y2)/((float)(n-1)), h);
					A[(x1*n+y1)*n*n+(x2*n+y2)] = numIntergStiff(phi1, phi2 , 0, 1, 0, 1, n*10);
				}
			}
		}
	}
}

template<typename T>
void setMassMat(float *M, const T *phi, const int n){
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
						M[(x1*n+y1)*n*n+(x2*n+y2)] = numIntergMass(phi1, phi2 , 0, 1, 0, 1, n*10);
					
				}
			}
		}
	}
}

