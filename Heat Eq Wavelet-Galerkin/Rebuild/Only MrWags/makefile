OBJS = mrWags.o elementClass.o
LFLAGS = -lcublas -lGL -lGLU -lglut -arch=sm_50
CFLAGS = -c
CC = nvcc

all : $(OBJS) visualizeheat.cu
	$(CC) $(LFLAGS) $(OBJS) visualizeheat.cu -o visualizeheat

mrWags.o : mrWags.h mrWags.cu elementClass.h elementClass.cpp
	$(CC) $(CFLAGS) $(LFLAGS) mrWags.cu

elementClass.o : elementClass.h elementClass.cpp
	$(CC) $(CFLAGS) $(LFLAGS) elementClass.cpp

clean:
	\rm *.o *~visualizeheat

tests: tests.cu
	$(CC) $(LFLAGS) tests.cu -o tests
