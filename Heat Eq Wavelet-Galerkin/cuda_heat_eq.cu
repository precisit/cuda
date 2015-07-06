#include "cuda_heat_eq.h"

//extern float *GLOBAL_A(nullptr);
extern int GLOBAL_N = 0;
extern float* U_0(NULL);


void prvec(float *A, int n){
	for(int x=0; x<n; x++){
		std::cout<<A[x]<<"   ";
	}
}

void prmat(float *A, int n){
	for(int x=0; x<n; x++){
		for(int y=0; y<n; y++){
			std::cout<<A[x+y*n]<<"   ";
		}
		std::cout<<std::endl;
	}
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

void setStiffMat(float *A, const class Element *phi, const int n){
	const float h = phi->h; 
	int d;

	for(int x1=0; x1<n; x1++){
		for(int y1=0; y1<n; y1++){
			for(int x2=0; x2<n; x2++){
				for(int y2=0; y2<n; y2++){
					d = dist(x1,y1,x2,y2);
					if(d == 0){
						A[(x1*n+y1)*n*n+(x2*n+y2)] = 4.0f/(h*h); //fett ad hoc! FIX! //O dessutom stämmer det nog inte. Eller?
					}
				}
			}
		}
	}
}

void setMassMat(float *M, const class Element *phi, const int n){
	//This is stupid and stuff. And shouldn't be done on the CPU. Oh, well.
	//But that'll be fixed later. 
	//(And also it's done ad hoc for hat funcs)
	
	
	//Antar att alla har samma storlek. FIX!
	const float h = phi->h;
	int d;
	
	for(int x1=0; x1<n; x1++){
		for(int y1=0; y1<n; y1++){
			for(int x2=0; x2<n; x2++){
				for(int y2=0; y2<n; y2++){
					d = dist(x1,y1,x2,y2);
					if(d == 0){
						M[(x1*n+y1)*n*n+(x2*n+y2)] = h*h/3.0f; //wolfram alpha said so. http://www.wolframalpha.com/input/?i=4*int_0%5Eh+int_0%5E%28h-x%29+%281-%28x%2By%29%2Fh%29%5E2+dy+dx
					}
					else if(d == 1){
						M[(x1*n+y1)*n*n+(x2*n+y2)] = h*h/24.0f; //http://www.wolframalpha.com/input/?i=4*int_0%5E%28h%2F2%29+int_0%5E%28x%29+%281-%28x%2By%29%2Fh%29*%281-%28-%28x-h%29%2By%29%2Fh%29+dy+dx
					}
				}
			}
		}
	}
}

void waveGal(float *U_0, const int n, const float dt, const float endTime, const float tol){


	::U_0 = U_0;


	::GLOBAL_N = n;
	assert(endTime>0.0f);
	assert(n>0);
	
	float t=0;
	
	//Setup the Element (should be a wavelet, but isn't yet. FIX!)
	Element phi = Element(0.0f, 0.0f, 1.0f/((float)(n-1)));
	float h = phi.h;
	
	/*
	//For debugging. Remove. FIX!
	float *Q;
	Q = (float*) calloc (n*n*n*n,sizeof(*Q));
	*/
	
	//Setup mass and stiffness matrix;
	float *M,*A;
	M = (float*) calloc (n*n*n*n,sizeof(*M));
	A = (float*) calloc (n*n*n*n,sizeof(*A));
	setMassMat(M, &phi,n);
	setStiffMat(A, &phi,n);
	
	//Move everything to the GPU
	float *devM, *devA, *devU, *devB;
	
    cublasHandle_t handle;
    cublasCreate(&handle); 
	
	cudaMalloc ((void**)&devM, n*n*n*n*sizeof(*devM));
	cudaMalloc ((void**)&devA, n*n*n*n*sizeof(*devA));
	cudaMalloc ((void**)&devU, n*n*sizeof(*devU));
	cudaMalloc ((void**)&devB, n*n*sizeof(*devB));
	
	cublasSetMatrix(n*n,n*n, sizeof(*A), A, n*n, devA, n*n);
	cublasSetMatrix(n*n,n*n, sizeof(*M), M, n*n, devM, n*n);
	cublasSetVector(n*n, sizeof(*U_0), U_0, 1, devU, 1);
	
	free(A); free(M);
	
	//Integrate forward through time
	//Antar att M och A är konstanta. Vilket de är än så länge. 
	//Adaptive grids kommer dessvärre ändra det. 
	//FIX!
	
	//Reassign A as A := M+k*A; //FIX!
	cublasSgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, &h, devA, n, 0, devM, n, devA, n);

	const float one = 1.0f;
	float zero = 0.0f;
	
	display();
	while(t<endTime){
		//b = M_old*U_old;
		cublasSgemv(handle, CUBLAS_OP_N, n*n, n*n,&one, devM, n*n, devU, 1, &zero, devB, 1); //M är förmodligen bandmarix. Använd annan funk!! //FIX
		
		//Solve A*u = b
		matSolve(devA, devU, devB, n*n, 0.1f, handle);
		
		cublasGetVector(n*n, sizeof(float), devU, 1, U_0,1); //Det här är segt. FIX!

		::U_0 = U_0;
		
		//måla här.
		display();
		
		t += dt;
		return;
	}
	
	cublasGetVector(n*n, sizeof(*U_0), devU, 1, U_0, 1);
	
	cudaFree(devM);
	cudaFree(devA);
	cudaFree(devU);
	cudaFree(devB);
	
	cublasDestroy(handle);
	
}






#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <time.h>


void color(int x, int array[]){
		if (x==0){
			//std::cout<<"red"<<std::endl;
			array[0] = 1.0; array[1] = 0.0; array[2] = 0.0;	
						
		}			
		else if (x==1){
			std::cout<<"not red"<<std::endl;
			array[0] = 1.0; array[1] = 0.5; array[2] = 0.0;
			
		}
		else if (x==2){
			std::cout<<"not red"<<std::endl;
			array[0] = 1.0; array[1] = 1.0; array[2] = 0.0;
			
		}
		else if (x==3){
			std::cout<<"not red"<<std::endl;
			array[0] = 0.5; array[1] = 1.0; array[2] = 0.0;
			
		}
		else if (x==4){
			std::cout<<"not red"<<std::endl;
			array[0] = 0.5; array[1] = 1.0; array[2] = 0.5;
			
		}
		else if (x==5){
			std::cout<<"not red"<<std::endl;
			array[0] = 0.0; array[1] = 0.5; array[2] = 0.0;
			
		}
		else if (x==6){
			std::cout<<"not red"<<std::endl;
			array[0] = 0.0; array[1] = 1.0; array[2] = 1.0;
			
		}
		else {
			std::cout<<"not red"<<std::endl;
			array[0] = 0.0; array[1] = 0.0; array[2] = 1.0;
		}
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);			//clear buffers to preset values (Indicates the buffers currently enabled for color writing.)
	//int number = 100;
	int r;
	int g;
	int b;	
	int colorarray[3];
	//float m [number];
	
	int n = ::GLOBAL_N;
	
	int val;
	float val_f;
	
	glBegin(GL_QUADS);
	for (int x = 0; x < n; x ++){
		for (int y = 0; y < n; y ++){
		//std::cout<<"tjeeeeeeeeeena!"<<std::endl;
			glVertex2f(x, y);  
			glVertex2f(x+1, y);
			glVertex2f(x+1, y+1);
			glVertex2f(x, y+1);

			//for (int k = 0; k < number; k++) {
				//int random = rand()%8;
				
				val_f = ::U_0[x*n+y];
				val = (int) (val_f*7.0f);
				//m[j] = random; 				
				color(val, colorarray);
				r = colorarray[0];
				g = colorarray[1];
				b = colorarray[2];
				glColor3f(r, g, b);
				//}
			//glMultMatrixf(m);
		}
	
			
	}
	
	//std::cout<<"klaaaaar"<<std::endl;
	glEnd();
	glFlush();
			
}
		
	

	
	
	

/*glBegin(GL_QUADS);
	for (GLfloat i = 0; i<= 100; i ++){
		for (GLfloat j = 0; j<= 100; j ++){
		glVertex2f(i*2, j*2);
		glVertex2f(i, j);
		color(4, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			glColor3f(r, g, b);
	}
	}
	glEnd();*/




void init() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, ::GLOBAL_N, ::GLOBAL_N, 0, -1, 1);
}

void timer(int v) {
	glutDisplayFunc(display);
	glutPostRedisplay();
	glutTimerFunc(v, timer, v);
}


/*
void reshape() {
}*/

/*

int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); 	
	glutInitWindowPosition(200, 200);		
	glutInitWindowSize(800, 800);		
	glutCreateWindow("Test of Grid");
	init();			
	//glutDisplayFunc(display);
	//glutReshapeFunc(reshape);
	glutTimerFunc(100, timer, 0);			
	glutMainLoop();				
}

*/













/*
int main(){
	int n = 3;
	float* U_0;
	U_0 = (float*) calloc (n*n,sizeof(*U_0));
	
	U_0[n+1] = 1.0f;

	waveGal(U_0, n, 0.02f, 0.04f, 0.1f);

	for(int i=0; i<n; i++){
		std::cout<<U_0[i]<<std::endl;
	}

	std::cout<<"it works?"<<std::endl;
	
	return 0;
}

*/
