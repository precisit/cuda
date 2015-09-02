#ifndef WAVELET_H
#define WAVELET_H

#include "node.h"

__device__ bool corner_node(Node* matrix, int idx, int idy, int LEN_OF_MATRIX);

__device__ bool interpol_dot_c(datatype p1, datatype p2, datatype p3);

__device__ int pow_c(int layer, int steplen);

__global__ void wavelet_kernel(Node* matrix, int LEN_OF_MATRIX, int step);

__global__ void layer_kernel(Node* matrix, int LEN_OF_MATRIX, int step);

void wavelet_compression(Node* matrix, int* countTrue);

#endif

