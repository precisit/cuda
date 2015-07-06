#include "visualizeheat.h"
#include "cuda_heat_eq.h"

int main(int argc, char** argv){
	int n = 100;
	float* U_0;
	U_0 = (float*) calloc (n*n,sizeof(*U_0));
	std::cout<<"main 1"<<std::endl;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); 	
	glutInitWindowPosition(200, 200);		
	glutInitWindowSize(800, 800);		
				








	
	U_0[n+1] = 1.0f;

	std::cout<<"main 2"<<std::endl;

	waveGal(U_0, n, 0.02f, 0.04f, 0.1f);

	for(int i=0; i<n; i++){
		//std::cout<<U_0[i]<<std::endl;
	}

	std::cout<<"it works?"<<std::endl;
	
	
	glutCreateWindow("Test of Grid");
	init();			
	//glutDisplayFunc(display);
	//glutReshapeFunc(reshape);
	glutTimerFunc(100, timer, 0);
	glutMainLoop();	
	
	std::cout<<"asså, klar på riktigt"<<std::endl;

	return 0;
}
