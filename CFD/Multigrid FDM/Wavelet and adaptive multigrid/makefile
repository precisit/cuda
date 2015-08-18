all: adaptive_grid.cpp node.cpp adaptive_multigrid_cuda.cu wavelet.cu wavelet_d_cuda.cu rungiz.cu main.cu makefile
	nvcc -arch=sm_50 -g node.cpp adaptive_grid.cpp adaptive_multigrid_cuda.cu wavelet.cu wavelet_d_cuda.cu rungiz.cu main.cu -o tillsammans

test: node.cpp test.cu
	nvcc -arch=sm_50 -g node.cpp adaptive_grid.cpp test.cu -o test
