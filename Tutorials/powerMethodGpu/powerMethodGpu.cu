#include "powerMethodGpu.h"

namespace precisit{
// square<T> computes the square of a number f(x) -> x*x
template <typename T>
struct square
{
    __host__ __device__
    T operator()(const T& x) const {
        return x * x;
    }
};

__host__
float _normGpu(thrust::device_vector<float>& d_x){
	/*
		Calculates the norm of a vector using the GPU (and Thrust obvs).
	*/
	
	// setup arguments
    square<float>        unary_op;
    thrust::plus<float> binary_op;
    float init = 0.0f;
	
	//calcs norm.
	return  std::sqrt(thrust::transform_reduce(d_x.begin(), d_x.end(), unary_op, init, binary_op) );
};

void _matVecMultGpu(const float *A, const float *x, float *y, const int n) {
	/*
		Multiplies the matrix A by the vector x and stores the result in the
		vector y. This is done on the GPU using CUBLAS.
		Assumes A is n-by-n and stored in column-major form.
	*/
	
	// Create a handle for CUBLAS
     cublasHandle_t handle;
     cublasCreate(&handle);
	
	//Seting up parameters in what seems to be the most stupid
	//way possible. Oh, well!
     const float alf = 1.0f;
     const float bet = 0.0f;
     const float *alpha = &alf;
     const float *beta = &bet;
     
	//Does the actual y:=Ax calculation
	cublasSgemv(handle,CUBLAS_OP_N, n, n, alpha,A, n, x,1, beta, y, 1);
	
	// Destroy the handle
     cublasDestroy(handle);
};

void eigPowerMethodGpu(const std::vector<float> A, std::vector<float> &eigVec, float &eigVal, int len){
	/*
	
		Uses the power method to calculate the eigen-vector (stored in eigVec) 
		that corresponds to the greatest eigen-value (stored in eigVal) of the 
		matrix A.
		A is assumed to be in column-major form. It is further assumed to be an
		len-by-len square matrix.
	
	*/
	
	//Copy the matrix to the GPU (aka device)
	thrust::device_vector<float> mat(A.begin(), A.end());
	
	//Takes a first (sucky) guess at the eigen-vector.
	thrust::device_vector<float> guess(len);
	thrust::fill(guess.begin(), guess.end(), 1.0f);
	thrust::device_vector<float> guess2(len);


    float norm;
	for(int iter=0; iter<200; ++iter){
	
		//Multiply the matrix mat with vector guess and store the result in guess2
		_matVecMultGpu(thrust::raw_pointer_cast(&mat[0]), thrust::raw_pointer_cast(&guess[0]),
				thrust::raw_pointer_cast(&guess2[0]), len);
		
		//Normalize guess2
		norm = _normGpu(guess2);
		thrust::transform(guess2.begin(),guess2.end(),thrust::make_constant_iterator(norm),
						guess2.begin(),thrust::divides<float>());
		
		//Multiply the matrix mat with vector guess2 and store the result in guess
		//Note that the pointers changed positions.
		_matVecMultGpu(thrust::raw_pointer_cast(&mat[0]), thrust::raw_pointer_cast(&guess2[0]),	
					thrust::raw_pointer_cast(&guess[0]), len);
								
		//Normalize guess
		norm = _normGpu(guess);
		thrust::transform(guess.begin(),guess.end(),thrust::make_constant_iterator(norm),
						guess.begin(),thrust::divides<float>());
	}

	//Copy the values back from the GPU.
	//(Den här delen är dum o kan göras stört mkt snabbare. Jaja, det funkar.)
	thrust::host_vector<float> eig=guess;
	for(int i=0; i<len; i++){
		eigVec[i] = eig[i];
	}
	
	eigVal = norm;
}
}
