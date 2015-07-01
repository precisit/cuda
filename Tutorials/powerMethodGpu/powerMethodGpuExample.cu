#include "powerMethodGpu.h"

int main(void)
{
	int n=2;
	std::vector<float> mat_vec(n*n);
	mat_vec[0]=4;
	mat_vec[1]=6;
	mat_vec[2]=1;
	mat_vec[3]=3;
	
	
	std::vector<float> v_vec(n);
	float val=0.0f;
	
	precisit::eigPowerMethodGpu(mat_vec, v_vec, val, n);
	std::cout<<std::endl;
	std::cout<<"Eigenvector:"<<std::endl;
	for(int i=0; i<n; i++){
		std::cout<<" "<<v_vec[i]<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Eigenvalue:\n "<<val<<std::endl;
    return 0;
}

