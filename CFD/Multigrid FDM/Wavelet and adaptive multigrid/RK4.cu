#include "RK4.h"

/*#define datatype double
#define Re 5.0f

#include <iostream>

void calculateVortInterior(const datatype* stream, datatype* vort_out, const datatype* vort, const int n){

	const datatype dx = 1.0f/((float)n-1.0f);

	//Start with the internal points
	for(int x=1; x < n-1; x++){
		for(int y=1; y < n-1; y++){
			//This is correct. Trust me.
			vort_out[x*n + y] = -(stream[x*n +  y+1]-stream[x*n + y-1])*(vort[(x+1)*n + y]-vort[(x-1)*n + y])/(2.0f*2.0f * dx * dx)+(stream[(x+1)*n + y]-stream[(x-1)*n + y])*(vort[x*n + y+1]-vort[x*n + y-1])/(2.0f*2.0f * dx * dx)+1.0f/(Re * dx * dx) * (-4.0f*vort[x*n + y]+vort[x*n + y+1]+vort[x*n + y-1]+vort[(x+1)*n + y]+vort[(x-1)*n + y]);
		}
	}
}



void calculateVortExterior(const datatype* stream, datatype* vort_out, const int n){

	int x,y;
	const datatype h = 1.0f/((float)n-1.0f);

	x=0;
	for(y=0; y < n; y++){
		vort_out[x*n + y] = -(stream[0*n + y] - stream[1*n + y])*2.0f/(h*h);
	}
	x=n-1;
	for(y=0; y < n; y++){
		vort_out[x*n + y] = (stream[(n-2)*n + y] - stream[(n-1)*n + y])*2.0f/(h*h);
	}
	y=0;
	for(x=0; x < n; x++){
		vort_out[x*n + y] = -(stream[x*n + 0] - stream[x*n + 1])*2.0f/(h*h);
	}
	y=n-1;
	for(x=0; x < n; x++){
		vort_out[x*n + y] = (stream[x*n + n-2] - stream[x*n + n-1])*2.0f/(h*h)+1.0*2.0f/h;
	}

}*/

void RK4(datatype dt, Node* y_vector, int* origoArray, int countTrue) {

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
	std::cout<<"multigrid 1 \n";
	adaptive_multigrid(k1, origoArray, countTrue, LAYERS);


    for (i = 0; i < countTrue; i++){
    	k2[i].stream += k1[i].stream*dt/2.0f;
    	k2[i].vort += k1[i].vort*dt/2.0f;
    }
    std::cout<<"multigrid 2 \n";
	adaptive_multigrid(k2, origoArray, countTrue, LAYERS);


	for (i = 0; i < countTrue; i++){
    	k3[i].stream += k2[i].stream*dt/2.0f;
    	k3[i].vort += k2[i].vort*dt/2.0f;
    }
    std::cout<<"multigrid 3 \n";
	adaptive_multigrid(k3, origoArray, countTrue, LAYERS);
	std::cout<<"LET S GO ! \n";
	for (i = 0; i < countTrue; i++){
    	k4[i].stream += k3[i].stream*dt;
    	std::cout<<" "<<k3[i].stream<<" ";
    	k4[i].vort += k3[i].vort*dt;
    	std::cout<<" "<<k3[i].vort<<" ";
    }
    std::cout<<"multigrid 4 \n";
	adaptive_multigrid(k4, origoArray, countTrue, LAYERS);

	std::cout<<"summation in RK4 \n";
    for (i = 0; i < countTrue; i++){
    	y_vector[i].vort += (dt/6.0f)*(k1[i].vort + 2.0f*k2[i].vort + 2.0f*k3[i].vort + k4[i].vort);
    	y_vector[i].stream += (dt/6.0f)*(k1[i].stream + 2.0f*k2[i].stream + 2.0f*k3[i].stream + k4[i].stream);
    }

	free(k1);
	free(k2);
	free(k3);
	free(k4);
	
	std::cout<<"multigrid & RK4 done! \n";
}
