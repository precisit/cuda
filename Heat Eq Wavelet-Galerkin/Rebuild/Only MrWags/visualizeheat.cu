#ifdef __APPLE_CC__
#include <iostream>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <time.h> 
#include "mrWags.h"
#include <unistd.h>

MrWags *globalMrWags; //GLÖÖÖÖM INTE, SIMBA!

void color(float x, float array[]){ 

   x = (x+0.2f)/1.2f; //Ta bort den här!

   if(x<0.0){
       x = 0.0;
   }

   array[0] = 1.0; array[1] = 1.0 - x; array[2] = 0.0; array[3] = 1.0;
   
}


void display() {
    std::cout<<"display()"<<std::endl;
    
    	glClear(GL_COLOR_BUFFER_BIT);            //clear buffers to preset values (Indicates the buffers currently enabled for color writing.)
    	//int number = 100;
	GLfloat r;
	GLfloat g;
	GLfloat b;
	GLfloat a;
	int i = 0;
    	int j;	
	float num;
	float colorarray[4];    	
	//int n = 50; // FIX!
    	const int n = ::globalMrWags->getN();

	//float U_array[n*n];    
	//float *array;
	//array = start(n, U_array);
	
	::globalMrWags->waveletGalerkinIter();

    glBegin(GL_QUADS);

    for (GLfloat x = 1; x < n-1; x ++){
		j=0;
		for (GLfloat y = 1; y < n-1; y ++){	
		
			num = ::globalMrWags->getU(i+1,j);
			color(num, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			a = colorarray[3];
			glColor4f(r, g, b, a);	
			glVertex2f(x+1, y);
				
			num = ::globalMrWags->getU(i,j);
			color(num, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			a = colorarray[3];
			glColor4f(r, g, b, a);
			glVertex2f(x, y); 
				
			num = ::globalMrWags->getU(i,j+1);
			color(num, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			a = colorarray[3];
			glColor4f(r, g, b, a);			
			glVertex2f(x, y+1);
					
			num = ::globalMrWags->getU(i+1,j+1);
			color(num, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			a = colorarray[3];
			glColor4f(r, g, b, a);
			glVertex2f(x+1, y+1);
			
			
			j++;		
		}				
		i++;	
	}
	glEnd();
	//glutSwapBuffers();
	glFlush();
	unsigned int micro = 100000;
    	usleep(micro);
	//glutPostRedisplay();        
    	display();
}

void init() {
    glClearColor (0.0, 0.0, 0.0, 1.0);

    glMatrixMode(GL_PROJECTION);                //specify which matrix is the current matrix (Applies subsequent matrix operations to the projection matrix stack)
    glLoadIdentity();
    
    const int n = ::globalMrWags->getN();
    glOrtho(0, n, n, 0, -1, 1);
	glShadeModel(GL_SMOOTH); 
}

/*
void timer(int v) {
    glutDisplayFunc(display);                //sets the display callback for the current window. When GLUT determines that the normal plane for the window needs to be redisplayed, the display callback for the window is called.
    
    glutPostRedisplay();                    //marks the current window as needing to be redisplayed.
    glutTimerFunc(v, timer, v);                //registers a timer callback to be triggered in a specified number of milliseconds.
}
/*
void reshape() {
}*/

int main(int argc, char** argv) {

    float *U_0;
    int n = 40;
    U_0 = (float*) calloc (n*n,sizeof(*U_0));
    for(int i=0; i<n*n; i++){
    	U_0[i] = 0.0f;
    }

	for(int x=n/2-4; x<n/2+4; x++){
       for(int y=n/2-4; y<n/2+4; y++){
           U_0[x*n+y]=1.0f;
       }
   }
    
    
    globalMrWags = new MrWags(U_0, n, 0.00001f, 1.0f, 0.0001f);

    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            std::cout<<"U: "<<::globalMrWags->getU(i,j)<<std::endl;
        }
    }

    glutInit(&argc, argv);                    //initialize the GLUT library.
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);         //sets the initial display mode (Bit mask to select a single buffered window; Bit mask to select an RGBA mode window)
    glutInitWindowPosition(0, 0);            //set the initial window position.
    glutInitWindowSize(2000, 1500);                //set the initial window size.
    glutCreateWindow("Simulation of Heat equation");    //Create a top-level window, give de window a name.
    init();            
    glutDisplayFunc(display);
    //glutReshapeFunc(reshape);
    //glutTimerFunc(1000, timer, 0);                //registers a timer callback to be triggered in a specified number of milliseconds.        
    glutMainLoop();                        //enters the GLUT event processing loop.
}
