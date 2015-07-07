#include "elementClass.h"

//The constructor
Element::Element(float x, float y, float h_in){
	xc = x;
	yc = y;
	
	h = h_in;
	h_inv = 1.0f/h;
};

//The actual function.
//Lets start with a (weird-looking) simple hat func.
//This should of course be a wavelet function. FIX!
float Element::f(float x, float y){
	//float tmp = (x-xc)+(y-xc);
	float dist = std::max(std::abs(x-xc),std::abs(y-yc)); //l_inf-norm
	if(dist > h){
		return 0.0f;
	}
	else{
		x = x - xc;
		y = y - yc;
		
		if(y>x){
			if(y>-x){
				return 1.0f-y*h_inv;
			}
			else{
				return 1.0f+x*h_inv;
			}
		}
		else{
			if(y>-x){
				return 1.0f-x*h_inv;
			}
			else{
				return 1.0f+y*h_inv;
			}
		}
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


