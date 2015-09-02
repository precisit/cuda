all:
	nvcc -lcublas -arch=sm_50 grid.cpp multigrid.cpp -o MrWags
unit_test: node.cpp adaptive_grid.cpp adaptive_multigrid.cpp adaptive_multigrid_cuda.cu unit_tests.cu
	nvcc -arch=sm_50 -g unit_tests.cu -o unitTests
multigridAdaptive: node.cpp adaptive_grid.cpp
	nvcc -arch=sm_50 -g adaptive_multigrid.cpp -o multigridAdaptive

cudaMultigridAdaptive: node.cpp adaptive_grid.cpp adaptive_multigrid.cpp adaptive_multigrid_cuda.cu
	nvcc -arch=sm_50 -g adaptive_multigrid_cuda.cu -o cudaMultigridAdaptive
