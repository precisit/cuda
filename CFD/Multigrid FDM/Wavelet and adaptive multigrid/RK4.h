#ifndef RK4_H
#define RK4_H

#include "adaptive_multigrid_cuda.h"
#include "node.h"
#include "parameters.h"


void RK4(const datatype dt, Node* y_vector, int* origoArray, int countTrue);


#endif