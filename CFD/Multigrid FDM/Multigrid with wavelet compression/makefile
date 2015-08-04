all:
	nvcc -arch=sm_50 node.cpp wavelet.cu -o wavelet

adaptive_grid:
	nvcc -arch=sm_50 -g node.cpp adaptive_grid.cpp -o adaptive_grid_test