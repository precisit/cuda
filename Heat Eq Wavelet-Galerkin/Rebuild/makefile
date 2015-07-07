all: makefile main.cu mrWags.cu elementClass.cpp matrixFuncs.cu
	nvcc -lcublas -lGL -lGLU -lglut -arch=sm_50 main.cu mrWags.cu elementClass.cpp matrixFuncs.cu -o heatEqVis
	
