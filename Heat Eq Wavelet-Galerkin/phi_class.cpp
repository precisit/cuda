#include "ElementClass.h"

//The constructor
Element::Element(float x, float y, float h_in){
	xc = x;
	yc = y;
	
	h = h_in;
	neg_h_inv = -1.0f/h;
};

//The actual function.
//Lets start with a (weird-looking) simple hat func.
//This should of course be a wavelet function. FIX!
float Element::f(float x, float y){
	float tmp = abs(x-xc)+abs(y-xc);
	if(tmp > h){
		return 0.0f;
	}
	else{
		return tmp*neg_h_inv+1.0f; 
	}
};

//The derivative of the function wrt x.
//Based on central difference (O(h^2)) 
float Element::f_x(const float x, const float y, const float dx){
	return (this->f(x+dx,y)-this->f(x-dx,y))/(2*dx);
};


//The derivative of the function wrt y.
//Based on central difference (O(h^2)) 
float Element::f_y(const float x, const float y, const float dy){
	return (this->f(x,y+dy)-this->f(x,y-dy))/(2*dy);
};


