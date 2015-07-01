/*
	This code provides the eigPowerMethodGpu-function for calculating the 
	eigenvector corresponding to the greatest eigenvalue of a matrix using
	CUDA. The algorithm used is a naive power method.
	
	It relies heavily on both CUBLAS and Thrust.
		
	@Viktor Wase (30-Jun-2015)
*/


#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/fill.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/host_vector.h>

#include <cmath>
#include <assert.h>
#include <iostream>
#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <vector>

namespace precisit{
	template <typename T>
	struct square;

	__host__
	float _normGpu(thrust::device_vector<float>& d_x);
	void _matVecMultGpu(const float *A, const float *x, float *y, const int n);
	void eigPowerMethodGpu(const std::vector<float> A, std::vector<float> &eigVec, float &eigVal, int len);
}
