#ifndef PRECISIT_WAVEGAL_H
#define PRECISIT_WAVEGAL_H

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <assert.h>
#include <iostream>

#include "elementClass.h"
#include "visualizeheat.h"

//Returns the distance (measured in the chess king's norm) between 2 points.
__host__ __device__
int dist(const int x1, const int y1, const int x2, const int y2);

//Sets up the Mass matrix M according to stuff.
__host__
void setMassMat(float *M, const class Element *phi, const int n);

//Sets up the Stiffness matrix A according to stuff.
__host__
void setStiffMat(float *A, const class Element *phi, const int n);

//Solves Ax=b for x. (A is n-by-n).
//Uses a hella lot of CUBLAS operations to do a Gauss-Seidel iteration.
//(Should probablably implement a Multigrid solver some time.)
__host__
void matSolve(float *A, float *x, float *b, const int n);

//Solves the Heat Eq. using a Wavelet-Galerkin method.
__host__
void waveGal(float *U_0, const int n, const float dt, const float endTime, const float tol);

//Does the same thing but sets the tolerance tol to a default value. (0.001?)
__host__
void waveGal(float *U_0,const int n, const float dt, const float endTime);

#endif
