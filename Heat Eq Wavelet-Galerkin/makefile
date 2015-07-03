all:
	nvcc -lcublas -arch=sm_50 cuda_heat_eq.cu elementClass.cpp -o cuda_heat_eq
	
test:
	nvcc -lcublas -arch=sm_50 unit_tests.cu elementClass.cpp -o unit_tests
	
vis: makefile
	nvcc -lGL -lGLU -lglut visualizeheat.cu -o vis
	
visHeat: makefile
	nvcc -lcublas -lGL -lGLU -lglut -arch=sm_50 visualizeheat.cu cuda_heat_eq.cu elementClass.cpp -o cuda_heat_eq
	
