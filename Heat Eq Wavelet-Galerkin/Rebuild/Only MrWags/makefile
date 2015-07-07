OBJS = mrWags.o elementClass.o
LFLAGS = -lcublas -lGL -lGLU -lglut -arch=sm_50
CFLAGS = -c
CC = nvcc

all : $(OBJS) visualizeheat.cu
	$(CC) $(LFLAGS) $(OBJS) visualizeheat.cu -o visualizeheat

mrWags.o : mrWags.h mrWags.cpp elementClass.h elementClass.cpp
	$(CC) $(CFLAGS) $(LFLAGS) mrWags.cpp

elementClass.o : elementClass.h elementClass.cpp
	$(CC) $(CFLAGS) $(LFLAGS) elementClass.cpp

clean:
	\rm *.o *~visualizeheat
