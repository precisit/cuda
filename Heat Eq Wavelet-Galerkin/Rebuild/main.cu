#include "mrWags.h"

int main(int argc, char** argv){
	int n = 4;
	float* U_0;
	U_0 = (float*) calloc (n*n,sizeof(*U_0));
	U_0[n+1] = 1.0f;
	
	MrWags a = MrWags(U_0, n,0.001f, 1.0f, 0.01f);
	return 0;	
}
