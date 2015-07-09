 /*
 Mr. Wags is a Multi-Resolution WAvelet Galerkin Simulation.
 
 */
 
#ifndef PRECISIT_MRWAGS_H
#define PRECISIT_MRWAGS_H

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <assert.h>
#include <iostream>
#include <stdio.h>


#include <string>

#include "elementClass.h"
//#include "matrixFuncs.h"


#include <time.h>


class MrWags{
	private: 
		float* U, *devU, *devA, *devM, *devB, *devStorage;
		int n;
		float dt, endTime, tol, t;
		cublasHandle_t handle;
		
		void waveletGalerkinInit();
		void waveletGalerkinEnd();
		void prmatgpu(float *devA, int n);
		void prvecgpu(float *devA, int n);
		
		/*static void drawCallback(){
	currentInstance->draw();
}

void setupDrawCallback(){
	currentInstance = this;
	::glutDisplayFunc(MrWags::drawCallback);
}*/
		
		
		
		
		
		
		template<typename T>
		float numIntergStiff(T func1, T func2 , float x0, float x1, float y0, float y1, int N);

		template<typename T>
		float numIntergMass(T func1, T func2 , float x0, float x1, float y0, float y1, int N);

		//Returns the distance (measured in the chess king's norm) between 2 points.
		__host__ __device__
		int dist_l1(const int x1, const int y1, const int x2, const int y2);

		//Sets up the Mass matrix M according to stuff.
		template<typename Element>
		void setMassMat(float *M, const Element *phi, const int n);

		//Sets up the Stiffness matrix A according to stuff.
		template<typename Element>
		void setStiffMat(float *A, const Element *phi, const int n);

		//Solves Ax=b for x. (A is n-by-n).
		//Uses a hella lot of CUBLAS operations to do a Conjugate Gradient iteration.
		//(Should probablably implement a Multigrid solver some time.)
		__host__
		void matSolve(const float *A, float *x, float *b, const int n, const float tol, cublasHandle_t handle);



		void wait(int);
		void prvec(float* a, int n);
		void prmat(float* a, int n);

		
	public:
		void waveletGalerkinIter();
		int getN();
		float getU(const int x, const int y);
		MrWags(float* U_0, const int n, const float dt, const float endTime, const float tol);
};

#endif
