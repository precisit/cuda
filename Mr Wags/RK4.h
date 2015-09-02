/*
	This is the Runge Kutta 4:th order method. It is used 
	to numerically integrate in time.
*/

#ifndef RK4_H
#define RK4_H

#include "adaptive_multigrid_cuda_new.h"
#include "node.h"
#include "parameters.h"


void RK4(const datatype dt, Node* y_vector, int* origoArray, int countTrue);


#endif