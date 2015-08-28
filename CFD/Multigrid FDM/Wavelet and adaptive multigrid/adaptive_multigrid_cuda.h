#ifndef ADAPTIVE_MULTIGRID_CUDA_H
#define ADAPTIVE_MULTIGRID_CUDA_H

#include "adaptive_grid.h"
//#include <cuda.h>
#include "node.h"
#include "parameters.h"

typedef datatype (*func_def)(int, int);

__device__ int abs_dev(const int x);

__device__ bool is1StepFromBoundary(const Node * u, const int ind, const int maxGlobIndex);

__device__ bool isInCorner(const Node* u, const int ind, const int maxGlobIndex);

__device__ bool is2StepsFromBoundary(const Node * u, const int ind, const int maxGlobIndex);

__device__ datatype dev_invPow(datatype x, const int n);

__device__ datatype getLaplacianStream(const Node * u, const int index1, const int index2, const datatype h,
  const int maxGlobIndex, const int layerNr, const int maxLayerNr);

__device__ datatype __BC_Function(const int x, const int y, const int maxLayer, const int maxGlobIndex);

__device__ int __findNodeGlobIdx(const Node* u, const int x, const int y, const int len);

__device__ int __pow_2(const datatype x);

__device__ datatype __findValOfClosestPoint(const Node* u, const int x, const int y, const int len);

__device__ datatype __interpolateGhostPoint(const Node* u, const Node* u_coarse, const int len_coarse, const int x, const int y, const int layerNr, const int maxLayer, const int maxGlobIndex, const datatype h);

__global__ void updateBFromInterpolation(Node* b, const Node* u, const int len, const Node* u_coarse, const int len_coarse, const int layerNr, const int maxLayer, const int maxGlobIndex, const datatype h);

__global__ void dev_updateBFromBoundary(Node* b, const Node* u, const int len, const int layerNr, const int maxLayer, const int maxGlobIndex, const datatype h);

void setupBoundaryOnCoarsestGrid(Node* bc, const int len , func_def BC_Function);

__global__ void calculateErrorLaplacian(const Node* b, const Node* u, Node* d, const int len, const datatype h, const int maxGlobIndex, const int layerNr, const int maxLayerNr);

__global__ void restrictMat(Node* u_coarse, const int len_coarse, Node* u_fine, const int len_fine );

__global__ void calculateRHS(const Node* u, const Node* d, Node* b, const int len, const datatype h, const int maxGlobIndex, const int layerNr, const int maxLayerNr);

__global__ void copy(Node* to, const Node* from, const int len);

__global__ void subtract_gpu(Node* a, const Node* b, const Node* c, const int len);

__global__ void add_gpu(Node* to, const Node* from, const int len);

__device__ datatype findNodeIdx(const Node* arr, const int x, const int y, const int n);

__device__ datatype findNodeVal(const Node* arr, const int x, const int y, const int n);

__global__ void interpolate(Node* u_fine, const Node* u_coarse, const int len_fine, const int len_coarse, const int layerNr, const int maxLayerNr, const int maxGlobIndex, const datatype h, Node* b);

__global__ void findNeighbours(Node* array, const int len );

__global__ void jacobiSmootherLaplacianStream(Node *from, Node * to, Node * b, const datatype h, const int len, const int maxGlobIndex, const int layerNr, const int maxLayerNr);

void ERRORCHECK();

__global__ void setVectorsToZero(Node * arr, const int len);

void multigrid_gpu(int k, AdaptiveGrid* grid, int pre, int sol, int post, const int maxGlobIndex, const int maxLayerNr);

__global__ void printAll(Node * arr, const int len);

__global__ void gpuPrint(Node* arr, const int len);

void move2gpu(AdaptiveGrid * grid);

void move2host(AdaptiveGrid * grid);

void adaptive_multigrid(Node* array, int* origoArray, int countTrue, int LAYERS);

__device__ datatype getLaplacianStreamNew(const Node* u, const int index1, const int index2, const datatype h, const int maxGlobIndex, const int layerNr, const int maxLayerNr);










__global__ void u_test(Node* x, const int len);
__global__ void findNeighbours2(Node* array, const int len );


#endif
