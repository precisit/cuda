#ifndef WAVELET_H
#define WAVELET_H
#include "parameters.h"
#include "node.h"

__device__ bool cornerNode(Node* matrix, int idx, int idy, int row, int colum);

__device__ bool interpolDot(float p1, float p2, float p3, float p4, float p5, float tol);

__global__ void wavelet_kernal(Node* matrix, int row, int colum, float tol, int step, int layers);

void wavelet_compression( Node* __restrict__ matrix, Node* __restrict__  ordedNodelist, int*  __restrict__ origoArray, int* __restrict__ countTrue);

#endif

