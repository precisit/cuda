all: adaptive_grid.cpp node.cpp adaptive_multigrid_cuda_new.cu wavelet_compression.cu wavelet_decompression.cu RK4.cu main.cu makefile
	nvcc -arch=sm_50 node.cpp adaptive_grid.cpp adaptive_multigrid_cuda_new.cu wavelet_compression.cu wavelet_decompression.cu RK4.cu main.cu -o MrWags 
