#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <time.h>

#ifndef PRECISIT_VISUAL_H
#define PRECISIT_VISUAL_H

void color(int x, int array[]);
void display();
void init() ;
void timer(int v);

#endif
