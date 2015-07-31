#define datatype double
#define Re 5.0f

#include <iostream>
#include <assert.h>

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

}

#ifdef NOT_DEBUGGING
void sofias_bajs_funk(const datatype* stream, datatype *f,const datatype* y_vector,const int n_sqrt){
	/*
		Input is poop on the form of a vector.
		Output is void. In the form om poop.
	*/

	calculateVortInterior(stream, f, y_vector, n_sqrt);
	calculateVortExterior(stream, f, n_sqrt);
}
#else
void sofias_bajs_funk(const datatype* stream, datatype *f,const datatype* y_vector,const int n_sqrt){
	assert(n_sqrt==2);
	//Solves y_t = 7y, x_t = x+y

	f[0] = y_vector[0] + y_vector[1];
	f[1] = 7.0f*y_vector[1];
}

#endif



void RK4_step(const datatype dt, datatype* y_vector, const datatype* stream, const int n, const int n_sqrt) {

    int i;
    datatype *y, *f, *k1, *k2, *k3, *k4;

    y = (datatype*) malloc(n*sizeof(datatype));
    f = (datatype*) malloc(n*sizeof(datatype));
    k1 = (datatype*) malloc(n*sizeof(datatype));
	k2 = (datatype*) malloc(n*sizeof(datatype));
	k3 = (datatype*) malloc(n*sizeof(datatype));
	k4 = (datatype*) malloc(n*sizeof(datatype));

	sofias_bajs_funk(stream, f, y_vector, n_sqrt);


	std::cout<<"y: ";
	for(int i= 0; i<n; i++){
		std::cout<<y_vector[i]<<std::endl;
	}


	std::cout<<"f: ";
	for(int i= 0; i<n; i++){
		std::cout<<f[i]<<std::endl;
	}
    std::cout<<std::endl;

    for (i = 0; i < n; i++){
    	k1[i] = dt * f[i];
    }

    for (i=0; i < n; i++){
    	y[i] = y_vector[i] + 0.5 * k1[i];
    }
    sofias_bajs_funk(stream, f, y_vector, n_sqrt);
    
    for (i = 0; i < n; i++){
    	k2[i] = dt * f[i];
	}
    
    for (i=0; i < n; i++){
    	y[i] = y_vector[i] + 0.5 * k2[i];
    }
    sofias_bajs_funk(stream, f, y_vector, n_sqrt);
    
    for (i = 0; i < n; i++)
    k3[i] = dt * f[i];
    
    for (i=0; i < n; i++){
   		y[i] = y_vector[i] + k3[i];
   	}
    sofias_bajs_funk(stream, f, y_vector, n_sqrt);
    for (i = 0; i < n; i++){
    	k4[i] = dt * f[i];
    }
    
    for (i = 0; i < n; i++){
    	y_vector[i] += (k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]) / 6.0f;
    }


    std::cout<<"k: ";
    std::cout<<k1[1]<<", "<<k2[1]<<", "<<k3[1]<<", "<<k4[1]<<std::endl;
	std::cout<<std::endl;


    std::cout<<"y_vector: ";
    for(int i= 0; i<n; i++){
		std::cout<<y_vector[i]<<std::endl;
	}
	std::cout<<std::endl;

	free(y);
	free(f);
	free(k1);
	free(k2);
	free(k3);
	free(k4);
}

/*
int main(){

	const int n = 25;
	const int n_sqrt = 5;
	const datatype dt = 0.1f;

	datatype *indata, *stream_in;
	indata = (datatype*) calloc(n, sizeof(datatype));

	stream_in = (datatype*) calloc(n, sizeof(datatype));

	for(int i= 0; i<n; i++){
		//indata[i] = ((datatype) (i % n_sqrt)) / ((float) n_sqrt);
		//stream_in[i] = ((datatype) (n_sqrt-i % n_sqrt)) / ((float) n_sqrt);

		indata[i] = rand() % 100;
		stream_in[i] = rand() % 100;
	}

	RK4_step(dt, indata, stream_in, n, n_sqrt);

	for(int i= 0; i<n; i++){
		//std::cout<<indata[i]<<std::endl;
	}

	free(indata);
	free(stream_in);
	return 0;
}
*/
