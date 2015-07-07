#ifndef PRECISIT_MATFUNC_H
#define PRECISIT_MATFUNC_H

template<typename T>
float numIntergStiff(T func1, T func2 , float x0, float x1, float y0, float y1, int N);

template<typename T>
float numIntergMsss(T func1, T func2 , float x0, float x1, float y0, float y1, int N);

//Returns the distance (measured in the chess king's norm) between 2 points.
__host__ __device__
int dist(const int x1, const int y1, const int x2, const int y2);

//Sets up the Mass matrix M according to stuff.
template<typename T>
void setMassMat(float *M, const T *phi, const int n);

//Sets up the Stiffness matrix A according to stuff.
template<typename T>
void setStiffMat(float *A, const T *phi, const int n);

//Solves Ax=b for x. (A is n-by-n).
//Uses a hella lot of CUBLAS operations to do a Gauss-Seidel iteration.
//(Should probablably implement a Multigrid solver some time.)
__host__
void matSolve(float *A, float *x, float *b, const int n);

#endif 
