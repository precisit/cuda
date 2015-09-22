/* 
 *An simple example file for visualization with OpenGL.
 *To compile use gcc -lGL -lGLU -lglut example_visualization.cpp
 */

#include <iostream>
#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <time.h> 
//#include "cuda_heat_eq.h"



void color(float x, float array[]){		
	x=x/10;
	if(x<0.0){
		x = 0.0;
	}

	array[0] = 1.0; array[1] = 1.0 - x; array[2] = 0.0; array[3] = 1.0;
}



void display() {
	glClear(GL_COLOR_BUFFER_BIT);			//clear buffers to preset values (Indicates the buffers currently enabled for color writing.)
	
	GLfloat r;
	GLfloat g;
	GLfloat b;
	GLfloat a;
	int i=0;
	int j;	
	float value;
	float colorarray[4];
	float array [40][40];
	//float num;
	
	//Create a random matrix with size 40x40
	for(int m = 0; m < 40; m++){
		//num=0;
		for(int n = 0; n < 40; n++){
			array[m][n]=rand()%10;
		}}
	
	
	glBegin(GL_QUADS);
	
	for (GLfloat x = 0; x < 40; x ++){
		j=0;
		for (GLfloat y = 0; y < 40; y ++){			
			value = array[i+1][j];
			color(value, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			a = colorarray[3];
			glColor4f(r, g, b, a);	
			glVertex2f(x+1, y);
				
			value = array[i][j];
			color(value, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			a = colorarray[3];
			glColor4f(r, g, b, a);
			glVertex2f(x, y); 
				
			value = array[i][j+1];
			color(value, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			a = colorarray[3];
			glColor4f(r, g, b, a);
			glVertex2f(x, y+1);
					
			value = array[i+1][j+1];
			color(value, colorarray);
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
	glutSwapBuffers();
	glFlush();
	glutPostRedisplay();
	
	
}

void init() {
	glClearColor (1.0, 0.0, 0.0, 1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, 40, 40, 0, -1, 1);
	glShadeModel(GL_SMOOTH); 
}

/*void timer(int v) {
	glutDisplayFunc(display);
	glutPostRedisplay();
	glutTimerFunc(v, timer, v);
}
*/

int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA); 	
	glutInitWindowPosition(200, 200);		
	glutInitWindowSize(1000, 1000);		
	glutCreateWindow("Test of Grid");
	init();			
	glutDisplayFunc(display);
	//glutReshapeFunc(reshape);
	//glutTimerFunc(100, timer, 0);			
	glutMainLoop();				
}
