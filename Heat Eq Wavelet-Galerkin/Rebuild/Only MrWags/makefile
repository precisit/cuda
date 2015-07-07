OBJS = mrWags.o elementClass.o
LFLAGS = -lcublas -lGL -lGLU -lglut -arch=sm_50
CFLAGS = -c
CC = nvcc

p1 : $(OBJS) main.cpp
	$(CC) $(LFLAGS) $(OBJS) main.cpp -o p1

mrWags.o : mrWags.h mrWags.cpp matrixFuncs.h matrixFuncs.cpp elementClass.h elementClass.cpp
	$(CC) $(CFLAGS) $(LFLAGS) mrWags.cpp

matrixFuncs.o : matrixFuncs.h matrixFuncs.cpp elementClass.h elementClass.cpp
	$(CC) $(CFLAGS) $(LFLAGS) matrixFuncs.cpp

elementClass.o : elementClass.h elementClass.cpp
	$(CC) $(CFLAGS) $(LFLAGS) elementClass.cpp

clean:
	\rm *.o *~ p1
	
simple:
	$(CC) $(LFLAGS) elementClass.cpp matrixFuncs.cpp mrWags.cpp main.cpp -o simple

old: makefile main.cu mrWags.cu elementClass.cpp matrixFuncs.cu
	nvcc -lcublas -lGL -lGLU -lglut -arch=sm_50 elementClass.cpp matrixFuncs.cu mrWags.cu main.cu -o heatEqVis
	
