#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <time.h> 
//#include "cuda_heat_eq.h"

#include "mrWags.h"

//extern float* start(int n, float U_array[]);

MrWags *globalMrWags; //GLÖÖÖÖM INTE, SIMBA!


void color(float x, float array[]){

        


        if (x>0.84 && x<=1.0){
            array[0] = 1.0; array[1] = 0.0; array[2] = 0.0;    
                        
        }            
        else if (x>0.70 && x<0.84){
            array[0] = 1.0; array[1] = 0.5; array[2] = 0.0;
            
        }
        else if (x>0.56 && x<=0.70){
            array[0] = 1.0; array[1] = 1.0; array[2] = 0.0;
            
        }
        else if (x>0.42 && x<=0.56){
            array[0] = 0.5; array[1] = 1.0; array[2] = 0.0;
            
        }
        else if (x>0.28 && x<=0.42){
            array[0] = 0.5; array[1] = 1.0; array[2] = 0.5;
            
        }
        else if (x>0.14 && x<=0.28){
            array[0] = 0.5; array[1] = 1.0; array[2] = 1.0;
            
        }
        else if (x>0.0 && x<=0.14){
            array[0] = 0.0; array[1] = 1.0; array[2] = 1.0;
            
        }
        else {
            array[0] = 0.0; array[1] = 0.0; array[2] = 1.0;
            
        }
}

void display() {
    std::cout<<"display()"<<std::endl;
    
    glClear(GL_COLOR_BUFFER_BIT);            //clear buffers to preset values (Indicates the buffers currently enabled for color writing.)
    //int number = 100;
    int r;
    int g;
    int b;    
    float colorarray[3];
    int i = 0;
    int j = 0;

    //int n = 50; // FIX!
    const int n = ::globalMrWags->getN();

    //float U_array[n*n];    
    //float *array;
    //array = start(n, U_array);

    ::globalMrWags->waveletGalerkinIter();

    float num;
    
    glBegin(GL_QUADS);
    for (GLfloat x = 0; x < n; x ++){ 
        i=0;
        for (GLfloat y = 0; y < n; y ++){
            glVertex2f(x, y);  
            glVertex2f(x+1, y);
            glVertex2f(x+1, y+1);
            glVertex2f(x, y+1);
            
            //num = array[i*n+j]; 
            num = ::globalMrWags->getU(i,j); //U[i*n+j];
            std::cout<<"num: "<<num<<std::endl;
            //int random = rand()%8;
            color(num, colorarray);
            r = colorarray[0];
            g = colorarray[1];
            b = colorarray[2];
            glColor3f(r, g, b);

            i++;
    
            //glShadeModel(GL_SMOOTH);
        }
        j++;
    }
    glEnd();
    glFlush();
}

void init() {
    glMatrixMode(GL_PROJECTION);                //specify which matrix is the current matrix (Applies subsequent matrix operations to the projection matrix stack)
    glLoadIdentity();
    
    const int n = ::globalMrWags->getN();
    glOrtho(0, n, n, 0, -1, 1);
}

/*void timer(int v) {
    glutDisplayFunc(display);                //sets the display callback for the current window. When GLUT determines that the normal plane for the window needs to be redisplayed, the display callback for the window is called.
    
    glutPostRedisplay();                    //marks the current window as needing to be redisplayed.
    glutTimerFunc(v, timer, v);                //registers a timer callback to be triggered in a specified number of milliseconds.
}

/*
void reshape() {
}*/

int main(int argc, char** argv) {

    float *U_0;
    int n = 4;
    U_0 = (float*) calloc (n*n,sizeof(*U_0));
    U_0[n*n/2] = 1.0f;
    
    globalMrWags = new MrWags(U_0, n, 0.001f, 1.0f, 0.05f);

    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            std::cout<<"U: "<<::globalMrWags->getU(i,j)<<std::endl;
        }
    }

    glutInit(&argc, argv);                    //initialize the GLUT library.
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);         //sets the initial display mode (Bit mask to select a single buffered window; Bit mask to select an RGBA mode window)
    glutInitWindowPosition(200, 200);            //set the initial window position.
    glutInitWindowSize(800, 800);                //set the initial window size.
    glutCreateWindow("Simulation of Heat equation");    //Create a top-level window, give de window a name.
    init();            
    glutDisplayFunc(display);
    //glutReshapeFunc(reshape);
    //glutTimerFunc(1000, timer, 0);                //registers a timer callback to be triggered in a specified number of milliseconds.        
    glutMainLoop();                        //enters the GLUT event processing loop.
}
