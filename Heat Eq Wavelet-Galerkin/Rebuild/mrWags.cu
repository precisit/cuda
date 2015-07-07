#include "mrWags.h"

MrWags::MrWags(float* U_0, const int n, const float dt, const float endTime, const float tol){
	this->n = n;
	this->dt = dt;
	this->endTime = endTime;
	this->tol = tol;
	
	this->t = 0.0f;
	
	cublasCreate(this->&handle);
	
	this->U = U_0;
	
	
	
	
	char *argv[] = {"myPixmap"};
    int argc = 1;
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); 	
	glutInitWindowPosition(200, 200);		
	glutInitWindowSize(800, 800);	
	
	
	glutCreateWindow("Test of Grid");
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, n, n, 0, -1, 1);	
	
	
	this->waveletGalerkinInit();	
	glutDisplayFunc(draw);
	//glutReshapeFunc(reshape);
	
	//glutTimerFunc(100, timer, 0);
	glutMainLoop();	
	
	
	
	

	this->waveletGalerkinEnd();
};


// This function is called every time GLUT refreshes the display.
void MrMAgs::draw(void)
{
	this->t += this.dt;

 	this->waveletGalerkinIter();
 	this->drawMatrix();

	assert(this->t<this->endTime);
}

void MrMags::drawMatrix(){
	glClear(GL_COLOR_BUFFER_BIT);			//clear buffers to preset values (Indicates the buffers currently enabled for color writing.)
	//int number = 100;
	int r;
	int g;
	int b;	
	int colorarray[3];
	//float m [number];
	
	int n = this->n;
	
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
				
				val_f = this->U[x*n+y];
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




void MrWags::waveletGalerkinInit(){
	
	assert(this.endTime>0.0f);
	assert(this.n>0);
	
	const int n = this.n;
	
	//Setup the Element (should be a wavelet, but isn't yet. FIX!)
	Element phi = Element(0.0f, 0.0f, 1.0f/((float)(n-1)));

	
	//Setup mass and stiffness matrix;
	float *M,*A;
	M = (float*) calloc (n*n*n*n,sizeof(*M));
	A = (float*) calloc (n*n*n*n,sizeof(*A));
	
	std::cout<<"funkar? 1"<<std::endl;
	setMassMat(M, &phi,n);
	std::cout<<"funkar? 2"<<std::endl;
	setStiffMat(A, &phi,n);
	
	std::cout<<"funkar? 3"<<std::endl;

	
	cudaMalloc ((void**)&this->devM, n*n*n*n*sizeof(*this->devM));
	cudaMalloc ((void**)&this->devA, n*n*n*n*sizeof(*this->devA));
	cudaMalloc ((void**)&this->devU, n*n*sizeof(*this->devU));
	cudaMalloc ((void**)&this->devB, n*n*sizeof(*this->devB));
	
	cublasSetMatrix(n*n,n*n, sizeof(*A), A, n*n, this->devA, n*n);
	cublasSetMatrix(n*n,n*n, sizeof(*M), M, n*n, this->devM, n*n);
	cublasSetVector(n*n, sizeof(*U_0), U_0, 1, this->devU, 1);
	
	free(A); free(M);
	
	float k = - dt;
	const float one = 1.0f;
	
	
	//Reassign A as A := M+k*A;
	cublasSgeam(this->handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, &k, this->devA, n, &one, this->devM, n, this->devA, n);
};



void MrWags::waveletGalerkinItter(){
	const int n= this->n;
	const float zero = 0.0f;
	
	//b = A*U_old;
	cublasSgemv(this->handle, CUBLAS_OP_N, n*n, n*n, &one, this->devA, n*n, this->devU, 1, &zero, devB, 1); //M är förmodligen bandmarix. Använd annan funk!! //FIX
	
	//Solve M*u = b
	matSolve(this->devM, this->devU, this->devB, n*n, this->tol ,this->handle);
		
	cublasGetVector(n*n, sizeof(float), this->devU, 1, this->U,1); //Det här är segt. FIX!

};


void MrWags::waveletGalerkinEnd(){
	cudaFree(this->devM);
	cudaFree(this->devA);
	cudaFree(this->devU);
	cudaFree(this->devB);
	
	cublasDestroy(this->handle);
}




