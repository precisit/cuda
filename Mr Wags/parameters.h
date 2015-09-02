/*
	These are the parameters used in the calculations. Some can be 
	altered, others cannot!
*/


#ifndef PARAMETERS_PRECISIT_H
#define PARAMETERS_PRECISIT_H

//These can be altered.
#define datatype float
#define Re 10.0f
#define KAPPA 1000000.0f
#define TOLERANCE 0.001f
#define THREADS_PER_BLOCK 10
#define LAYERS 4
#define ITERATIONS_UNTIL_GRID_UPDATE 2
#define END_TIME 0.0001f
#define DELTA_T 0.000005f



//LEAVE THESE ALONE!
#define VORTSTREAM
const int LEN_OF_MATRIX = (1<<LAYERS)+1;
#define STEP 2

#endif
