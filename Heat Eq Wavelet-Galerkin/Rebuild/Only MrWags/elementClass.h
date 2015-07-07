

#ifndef PRECISIT_ELEMENT_H
#define PRECISIT_ELEMENT_H

#include <cmath>
#include <algorithm>
#include <iostream>

class Element{
	private: 
		float xc,yc;
	public:
		float h, h_inv;
		Element(float x, float y, float h_in);
	
		//The actual function.
		//Lets start with a (weird-looking) simple hat func.
		//This should of course be a wavelet function. FIX!
		float f(float x, float y);
	
		//The derivative of the function wrt x.
		//Based on central difference (O(h^2)) 
		float f_x(const float x, const float y, const float dx);
	
		//The derivative of the function wrt y.
		//Based on central difference (O(h^2)) 
		float f_y(const float x, const float y, const float dy);
};

#endif
