#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <cmath>

#define datatype float
#define Re 10.0f
#define VORTSTREAM

class Node
{

	/*
		This is the data node at the end of the pointers.
		Each grid point should contain one of these.
	*/

public:
	//Position
	datatype x, y;
	bool isPicked;
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
	Node(){
		this->isPicked = false;
		this->layer = 0;
		this->x = 0.0f;
		this->y = 0.0f;

		this->x_index = -1;
		this->y_index = -1;

		x_index_global = -1;
		y_index_global = -1;

		this->nodeRight = NULL;
		this->nodeBelow = NULL;
		this->nodeAbove = NULL;
		this->nodeLeft = NULL;

		#ifdef VORTSTREAM
			this->vort = 0.0f;
			this->stream = 0.0f;
			this->stream_y = 0.0f;
		#else
			this->p = 0.0f;
			this->v_x = 0.0f;
			this->v_y = 0.0f;
		#endif
	};

	Node(const datatype x_in, const datatype y_in, const int x_index_in, const int y_index_in,
	 const datatype p_in, const datatype v_x_in, const datatype v_y_in){
		this->isPicked = false;
		this->layer = 0;
		this->x = x_in;
		this->y = y_in;

		this->x_index = x_index_in;
		this->y_index = y_index_in;

		//This should probably be changed? FIX!
		x_index_global = -1;
		y_index_global = -1;

		this->nodeRight = NULL;
		this->nodeBelow = NULL;
		this->nodeAbove = NULL;
		this->nodeLeft = NULL;

	#ifdef VORTSTREAM
		this->vort = p_in;
		this->stream = v_x_in;
		this->stream_y = v_y_in;
	#else
		this->p = p_in;
		this->v_x = v_x_in;
		this->v_y = v_y_in;
	#endif
	};

	Node(const datatype x_in, const datatype y_in, const int x_index_in, const int y_index_in,
		const int x_global_index_in, const int y_global_index_in,
	 const datatype p_in, const datatype v_x_in, const datatype v_y_in){
		this->isPicked = false;
		this->layer = 0;
		this->x = x_in;
		this->y = y_in;

		this->x_index = x_index_in;
		this->y_index = y_index_in;

		//This should probably be changed? FIX!
		x_index_global = x_global_index_in;
		y_index_global = y_global_index_in;

		this->nodeRight = NULL;
		this->nodeBelow = NULL;
		this->nodeAbove = NULL;
		this->nodeLeft = NULL;

	#ifdef VORTSTREAM
		this->vort = p_in;
		this->stream = v_x_in;
		this->stream_y = v_y_in;
	#else
		this->p = p_in;
		this->v_x = v_x_in;
		this->v_y = v_y_in;
	#endif
	};

	/*
	void operator=(const Node& rhs){
			this->x = rhs.x;
			this->y = rhs.y;
			this->x_index = rhs.x_index;
			this->y_index = rhs.y_index;
		#ifdef VORTSTREAM
			this->vort = rhs.vort;
			this->stream = rhs.stream;
			this->stream_y = rhs.stream_y;
		#else
			this->p = rhs.p;
			this->v_x = rhs.v_x;
			this->v_y = rhs.v_y;
		#endif
	};
	*/
	
	//Addition operator
	Node operator+(const Node& rhs){
		#ifdef VORTSTREAM
			return Node(this->x, this->y, this->x_index, this->y_index, this->x_index_global, this->y_index_global, this->vort+rhs.vort, this->stream+rhs.stream, this->stream_y+rhs.stream_y);
		#else
			return Node(this->x, this->y, this->x_index, this->y_index, this->p+rhs.p, this->v_x+rhs.v_x, this->v_y+rhs.v_y);
		#endif
	};

	//Division with a number
	Node operator/(const datatype& rhs){
		#ifdef VORTSTREAM
			return Node(this->x, this->y, this->x_index, this->y_index,this->x_index_global, this->y_index_global, this->vort/rhs, this->stream/rhs, this->stream_y/rhs);
		#else
			return Node(this->x, this->y, this->x_index, this->y_index, this->p/rhs, this->v_x/rhs, this->v_y/rhs);
		#endif
	};

};