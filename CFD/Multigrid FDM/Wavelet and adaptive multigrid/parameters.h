#ifndef PARAMETERS_PRECISIT_H
#define PARAMETERS_PRECISIT_H

//These can be altered.
#define datatype float
#define Re 10.0f
#define KAPPA 1000000.0f
#define TOLERANCE 0.2f
#define THREADS_PER_BLOCK 10
#define LAYERS 4



//LEAVE THESE ALONE!
#define VORTSTREAM
const int row = (1<<LAYERS)+1;
const int colum = (1<<LAYERS)+1;
const int step = 2;
const int layers = LAYERS;
const float tol = TOLERANCE;

#endif
