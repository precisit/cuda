/*
 * Wavelet decompression is a function decompressing the Wavelet compression
 * by doing the computations in reverse. All computations are done on the GPU.
 */
#ifndef WAVELET_D_H
#define WAVELET_D_H

#include "node.h"


__device__ datatype interpol_dot_d(datatype p1, datatype p2);

__device__ int pow_d(int layer, int steplen);

__global__ void wavelet_d_kernal_y(Node* matrix, int LEN_OF_MATRIX, int layer, int step_in);

__global__ void wavelet_d_kernal_x(Node* matrix, int LEN_OF_MATRIX, int layer, int step_in);

void wavelet_decompression(Node* __restrict__ array, Node* __restrict__ matrix, int *countTrue);


#endif

