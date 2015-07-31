#include "rungiz.cu"
#include <math.h>

datatype test_vortInternal(){
	const int n_sqrt = 5;
	const int n = n_sqrt * n_sqrt;

	const datatype h = 1.0f/((float)n_sqrt-1.0f);
	datatype *stream, *vort, *f;

	stream = (datatype* ) malloc(n*sizeof(datatype));
	vort = (datatype* ) malloc(n*sizeof(datatype));
	f = (datatype* ) malloc(n*sizeof(datatype));

	
	datatype x=0, y=0;
	for (int xx = 0; xx < n_sqrt; xx++)
	 {
	 	for (int yy = 0; yy < n_sqrt; yy++)
	 	{
	 		stream[xx*n_sqrt+yy] = x*y/6.0f*(x*x+y*y);
	 		vort[xx*n_sqrt + yy] = -2*x*y;
	 		y += h;
	 	}
	 	x += h;
	 	y = 0;
	 }
	
	calculateVortInterior(stream, f, vort, n_sqrt);

	 datatype error = 0.0f;

	
	x=h, y=h;
	for (int xx = 1; xx < n_sqrt-1; xx++)
	 {
	 	for (int yy = 1; yy < n_sqrt-1; yy++)
	 	{
	 		error += fabs((x*y/3.0f)*(-2.0*x*x+2.0*y*y)- f[xx*n_sqrt + yy]);
	 		y += h;
	 	}
	 	x += h;
	 	y = h;
	 }


	 free(stream);
	 free(vort);
	 free(f);

	 return error;
}

void test_RK4_step(){
	
	const datatype dt = 0.01;
	datatype* y_vector;
	y_vector = (datatype*) malloc(2*sizeof(datatype));
	y_vector[0] = 0.0f;
	y_vector[1] = 1.0f;
	const datatype stream [1]= {0.0};
	const int n = 2;
	const int n_sqrt = 2;

	RK4_step(dt, y_vector, stream, n, n_sqrt);
	free(y_vector);
}

int main(){

	//std::cout<<test_vortInternal()<<std::endl;
	test_RK4_step();
	//run_rungiz();

	return 0;
}
