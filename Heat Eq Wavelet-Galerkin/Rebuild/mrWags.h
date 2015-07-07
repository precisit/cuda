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

#include "elementClass.h"
#include "matrixFuncs.h"

class MrWags{
	private: 
		float* U, devU, devA, devM, devB;
		int n;
		float dt, endTime, tol, t;
		cublasHandle_t handle;
		
		void waveletGalerkinInit();
		void waveletGalerkinIter();
		void waveletGalerkinEnd();
		
		
		void draw(void);
		void drawMatrix();
		
		
		static MrWags* currentInstance;

		static void drawCallback(){
			currentInstance->draw();
		}

		void setupDrawCallback(){
			currentInstance = this;
			::glutDisplayFunc(MrWags::drawCallback);
		}
		
		
		
	public:
		MrWags(float* U_0, const int n, const float dt, const float endTime, const float tol);
};

#endif
