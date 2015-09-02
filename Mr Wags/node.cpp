#include "node.h"

//#define datatype float
//#define Re 10.0f
//#define VORTSTREAM


	Node::Node(){
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

	Node::Node(const datatype x_in, const datatype y_in, const int x_index_in, const int y_index_in,
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

	Node::Node(const datatype x_in, const datatype y_in, const int x_index_in, const int y_index_in,
		const int x_global_index_in, const int y_global_index_in,
	 const datatype p_in, const datatype v_x_in, const datatype v_y_in){
		this->isPicked = false;
		this->layer = 0;
		this->x = x_in;
		this->y = y_in;

		this->x_index = x_index_in;
		this->y_index = y_index_in;

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
	
	//Addition operator
	Node Node::operator+(const Node& rhs){
		#ifdef VORTSTREAM
			return Node(this->x, this->y, this->x_index, this->y_index, this->x_index_global, this->y_index_global, this->vort+rhs.vort, this->stream+rhs.stream, this->stream_y+rhs.stream_y);
		#else
			return Node(this->x, this->y, this->x_index, this->y_index, this->p+rhs.p, this->v_x+rhs.v_x, this->v_y+rhs.v_y);
		#endif
	};

	//Division with a number
	Node Node::operator/(const datatype& rhs){
		#ifdef VORTSTREAM
			return Node(this->x, this->y, this->x_index, this->y_index,this->x_index_global, this->y_index_global, this->vort/rhs, this->stream/rhs, this->stream_y/rhs);
		#else
			return Node(this->x, this->y, this->x_index, this->y_index, this->p/rhs, this->v_x/rhs, this->v_y/rhs);
		#endif
	};

