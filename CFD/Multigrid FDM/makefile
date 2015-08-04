all:
	nvcc -lcublas -arch=sm_50 grid.cpp multigrid.cpp -o multigridTest
unit_test: node.cpp adaptive_grid.cpp
	nvcc -arch=sm_50  unit_tests.cpp -o UnitTests