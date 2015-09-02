#ifndef NODE_H
#define NODE_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <cmath>

#include "parameters.h"

//#define datatype float
//#define Re 10.0f
//#define VORTSTREAM

class Node
{

	/*
		This is the data node at the end of the pointers.
		Each grid point should contain one of these.
	*/

public:
	//Position
	datatype x, y;
	int isPicked;
	int layer;

	//Index of Position
	int x_index, y_index;

	//Global index of position
	int x_index_global, y_index_global;

	Node *nodeAbove, *nodeBelow, *nodeRight, *nodeLeft;

	#ifdef VORTSTREAM
	//This is for the vorticity/stream formulation. 
	//Then the vector will be on the form (vort, stream, stream_y)^T 
	datatype vort, stream, stream_y;

	#else
	//Preasure and velocity.
	datatype p, v_x, v_y;

	#endif

	//This sets the vorticity and stream variables
	//void updateVorticityAndStream(Node* left, Node* right, Node* up, Node* down);

	//Constructors
	Node();

	Node(const datatype x_in, const datatype y_in, const int x_index_in, const int y_index_in,
	 const datatype p_in, const datatype v_x_in, const datatype v_y_in);

	Node(const datatype x_in, const datatype y_in, const int x_index_in, const int y_index_in,
		const int x_global_index_in, const int y_global_index_in,
	 const datatype p_in, const datatype v_x_in, const datatype v_y_in);
	
	//Addition operator
	Node operator+(const Node& rhs);

	//Division with a number
	Node operator/(const datatype& rhs);

};

#endif
