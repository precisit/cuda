#include "rungiz.cu"
#include <math.h>

datatype test_vortInternal(){
	const int n_sqrt = 17;
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


int main(){
	std::cout<<test_vortInternal()<<std::endl;
	//run_rungiz();

	return 0;
}
