#ifndef WAVELET_D_H
#define WAVELET_D_H
#include "parameters.h"
#include "node.h"

__device__ float checkInterPol(float p5, float dot);

__device__ float interpolDotEdge(float p1, float p2);

__device__ float interpolDotMiddle(float p1, float p2, float p3, float p4);

__device__ int pow(int layer, int steplen);

__global__ void wavelet_d_kernal(Node* __restrict__ matrix, Node* __restrict__ array, int row, int colum, int countTrue, int layers);

void wavelet_decompression(Node* array, Node* matrix, int *countTrue);


#endif

