all:
	nvcc -lcublas -arch=sm_50 grid.cpp multigrid.cpp -o multigridTest
unit_test: node.cpp adaptive_grid.cpp
	nvcc -arch=sm_50 -g unit_tests.cpp -o unitTests
multigridAdaptive: node.cpp adaptive_grid.cpp
	nvcc -arch=sm_50 -g adaptive_multigrid.cpp -o multigridAdaptive