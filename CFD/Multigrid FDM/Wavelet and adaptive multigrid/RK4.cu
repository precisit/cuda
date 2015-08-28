#include "RK4.h"

void RK4(datatype dt, Node* y_vector, int* origoArray, int countTrue) {
	/*
		The Runge Kutta 4 method is used here to solve the 
		incompressible navier stokes eq in 2d (using the 
		call to adaptive_multigrid_new function).
	*/

    int i;
    Node *k1, *k2, *k3, *k4;

    k1 = (Node*) malloc(countTrue*sizeof(Node));
	k2 = (Node*) malloc(countTrue*sizeof(Node));
	k3 = (Node*) malloc(countTrue*sizeof(Node));
	k4 = (Node*) malloc(countTrue*sizeof(Node));


	for (i = 0; i < countTrue; ++i)
	{
		k1[i] = y_vector[i];
		k2[i] = y_vector[i];
		k3[i] = y_vector[i];
		k4[i] = y_vector[i];
	}

	adaptive_multigrid_new(k1, origoArray, countTrue);
    for (i = 0; i < countTrue; i++){
    	k2[i].stream += k1[i].stream*dt/2.0f;
    	k2[i].vort += k1[i].vort*dt/2.0f;
    }
	adaptive_multigrid_new(k2, origoArray, countTrue); 
	for (i = 0; i < countTrue; i++){
    	k3[i].stream += k2[i].stream*dt/2.0f;
    	k3[i].vort += k2[i].vort*dt/2.0f;
    }

	adaptive_multigrid_new(k3, origoArray, countTrue);
	for (i = 0; i < countTrue; i++){
    	k4[i].stream += k3[i].stream*dt;
    	k4[i].vort += k3[i].vort*dt;
    }

	adaptive_multigrid_new(k4, origoArray, countTrue);

	assert(row == colum);
	for (i = 0; i < countTrue; i++){
    	if ((y_vector[i].x_index_global == 0 || y_vector[i].y_index_global == 0 || y_vector[i].x_index_global == row-1 ||
    	 y_vector[i].y_index_global == row-1) == false)
    	{
    		y_vector[i].vort = (dt/6.0f)*(k1[i].vort + 2.0f*k2[i].vort + 2.0f*k3[i].vort + k4[i].vort);
    		
    		//y_vector[i].vort += dt*k1[i].vort;
    	}
    	else{
    		//Since we have Dirichlet BC the boundary values are set exactly.
    		y_vector[i].vort = k4[i].vort;
    	}
    	//The stream function is updated every step, but only as a function of vorticity and not of time.
    	y_vector[i].stream = k4[i].stream;
    }

	free(k1);
	free(k2);
	free(k3);
	free(k4);
}
