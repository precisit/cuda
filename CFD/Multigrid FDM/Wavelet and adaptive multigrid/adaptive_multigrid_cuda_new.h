#ifndef ADAPTIVE_MULTIGRID_CUDA_NEW_H
#define ADAPTIVE_MULTIGRID_CUDA_NEW_H

#include "adaptive_grid.h"
#include "node.h"
#include "parameters.h"

void adaptive_multigrid_new(Node* array, int* origoArray, int countTrue);
__global__ void updateBNew(const Node* u, Node* b, const Node* u_coarse, const int len, const int len_coarse, const int maxGlobIndex);
__global__ void updateCoarsestB(const Node* u, Node* b, const Node* u_coarse, const int len, const int maxGlobIndex);
void multigrid_gpu(int k, AdaptiveGrid* grid, int pre, int sol, int post, const int maxGlobIndex, const int maxLayerNr);
__device__ datatype __BC_Function(const int x, const int y, const int maxGlobIndex);
__device__ bool isInBoundaryNew(const Node* u, const int ind);

bool isLeaf(Node* u, const int first, const int last, const int x_min, const int x_max, const int y_min, const int y_max);
void recursiveGridFill(Node* u, int layer, AdaptiveGrid** gridList, int* firstPtr, const int x_min, const int x_max, const int y_min, const int y_max, const float h);
void move2gpu(AdaptiveGrid * grid);
void move2host(AdaptiveGrid * grid);
void ERRORCHECK();

void calcVortFromStream(AdaptiveGrid* grid);
__global__ void vortInterior(Node* to, Node* from, const int len, const datatype h);
__global__ void vortExterior(Node* to, Node* from, const int len, const datatype h, const int maxGlobIndex);
__global__ void copyVort(Node* to, const Node* from, const int len);
__device__ datatype findNodeVal(const Node* arr, const int x, const int y, const int n);
__device__ datatype findNodeValVort(const Node* arr, const int x, const int y, const int n);
__device__ bool isInterior(const Node* arr, const int i);
__device__ bool is1StepFromBoundary(const Node * u, const int ind, const int maxGlobIndex);
__global__ void jacobiSmootherLaplacianStream(Node *from, Node * to, Node * b, const datatype h, const int len, const int maxGlobIndex, const int layerNr, const int maxLayerNr);
__global__ void calculateErrorLaplacian(const Node* b, const Node* u, Node* d, const int len, const datatype h, const int maxGlobIndex, const int layerNr, const int maxLayerNr);
__global__ void restrictMat(Node* u_coarse, const int len_coarse, Node* u_fine, const int len_fine );
__global__ void copy(Node* to, const Node* from, const int len);
__global__ void subtract_gpu(Node* a, const Node* b, const Node* c, const int len);
__global__ void findNeighbours(Node* array, const int len );
__device__ int abs_dev(const int x);
__device__ bool isVonNeumannBC(const int x, const int y);
__global__ void add_gpu(Node* to, const Node* from, const int len);
__device__ bool isOneStepIn(const Node *u, const int index1, const int index2);
__global__ void interpolate(Node* u_fine, const Node* u_coarse, const int len_fine, const int len_coarse, const int layerNr, const int maxLayerNr, const int maxGlobIndex, const datatype h, Node* b);
__device__ int findNodeIdx(const Node* arr, const int x, const int y, const int n);
__global__ void updateBNew(const Node* u, Node* b, const Node* u_coarse, const int len, const int len_coarse, const int maxGlobIndex, const int layerNr);

__device__ datatype getLaplacianStreamNew(const Node* u, const int index1, const int index2, const datatype h, const int maxGlobIndex, const int layerNr, const int maxLayerNr);

__device__ int __findNodeGlobIdx(const Node* u, const int x, const int y, const int len);
__device__ datatype __findNodeGlobStream(const Node* u, const int x, const int y, const int len);
__global__ void updateFromDirichletBC(Node* u, Node* b, Node* w, const int len, const int maxGlobIndex);

__global__ void gaussSeidelSmootherStream(Node* u, const Node* b, const int len, const int maxGlobIndex, const int layerNr, const datatype h);
datatype findClosestNodeValVortGlob(const Node* arr, const int x, const int y, const int n);
__global__ void copyVortExterior(Node* to, const Node* from, const int len);



#endif