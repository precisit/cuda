all: adaptive_grid.cpp node.cpp adaptive_multigrid_cuda.cu wavelet.cu wavelet_d.cu RK4.cu main.cu makefile
	nvcc -arch=sm_50 -rdc=true -g node.cpp adaptive_grid.cpp adaptive_multigrid_cuda.cu wavelet.cu wavelet_d.cu RK4.cu main.cu -o tillsammans -lcudadevrt

test: node.cpp test.cu adaptive_multigrid_cuda_new.cu adaptive_grid.cpp RK4.cu
	nvcc -arch=sm_50 -g node.cpp adaptive_grid.cpp adaptive_multigrid_cuda_new.cu RK4.cu test.cu -o test

wavelet: node.cpp wavelet_viktor.cu
	nvcc -arch=sm_50 -g node.cpp wavelet_viktor.cu -o wavelet_viktor
