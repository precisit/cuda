#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <time.h> 
void color(int x, int array[]){
	
		

		if (x==0){
			array[0] = 1.0; array[1] = 0.0; array[2] = 0.0;	
						
		}			
		else if (x==1){
			array[0] = 1.0; array[1] = 0.5; array[2] = 0.0;
			
		}
		else if (x==2){
			array[0] = 1.0; array[1] = 1.0; array[2] = 0.0;
			
		}
		else if (x==3){
			array[0] = 0.5; array[1] = 1.0; array[2] = 0.0;
			
		}
		else if (x==4){
			array[0] = 0.5; array[1] = 1.0; array[2] = 0.5;
			
		}
		else if (x==5){
			array[0] = 0.0; array[1] = 0.5; array[2] = 0.0;
			
		}
		else if (x==6){
			array[0] = 0.0; array[1] = 1.0; array[2] = 1.0;
			
		}
		else {
			array[0] = 0.0; array[1] = 0.0; array[2] = 1.0;
			
		}

	 

	

}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);			//clear buffers to preset values (Indicates the buffers currently enabled for color writing.)
	//int number = 100;
	int r;
	int g;
	int b;	
	int colorarray[3];
	//float m [number];
	
	
	glBegin(GL_QUADS);
	for (GLfloat x = 0; x < 10; x ++){
		for (GLfloat y = 0; y < 10; y ++){
			glVertex2f(x, y);  
			glVertex2f(x+1, y);
			glVertex2f(x+1, y+1);
			glVertex2f(x, y+1);

			//for (int k = 0; k < number; k++) {
				int random = rand()%8;
				//m[j] = random; 				
				color(random, colorarray);
				r = colorarray[0];
				g = colorarray[1];
				b = colorarray[2];
				glColor3f(r, g, b);
				//}
			//glMultMatrixf(m);
		}
	
			
	}
	glEnd();
	glFlush();
		
		
}
		
	

	
	
	

/*glBegin(GL_QUADS);
	for (GLfloat i = 0; i<= 100; i ++){
		for (GLfloat j = 0; j<= 100; j ++){
		glVertex2f(i*2, j*2);
		glVertex2f(i, j);
		color(4, colorarray);
			r = colorarray[0];
			g = colorarray[1];
			b = colorarray[2];
			glColor3f(r, g, b);
	}
	}
	glEnd();*/




void init() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, 100, 100, 0, -1, 1);
}

/*void timer(int v) {
	glutDisplayFunc(display);
	glutPostRedisplay();
	glutTimerFunc(v, timer, v);
}


/*
void reshape() {
}*/

int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); 	
	glutInitWindowPosition(200, 200);		
	glutInitWindowSize(800, 800);		
	glutCreateWindow("Test of Grid");
	init();			
	glutDisplayFunc(display);
	//glutReshapeFunc(reshape);
	//glutTimerFunc(100, timer, 0);			
	glutMainLoop();				
}
