#define datatype float

class Node
{

	/*
		This is the data node at the end of the pointers.
		Each grid point should contain one of these.
	*/

public:
	//Position
	datatype x, y;

	//Index of Position
	int x_index, y_index;

	//Preasure and velocity.
	datatype p, v_x, v_y;


	datatype vort, stream;

	//This sets the vorticity and stream variables
	updateVorticityAndStream(Node* left, Node* right, Node* up, Node* down){
		
	}

	//Constructors
	Node(){
		this->x = 0.0f;
		this->y = 0.0f;

		this->x_index = 0;
		this->y_index = 0;

		this->p = 0.0f;
		this->v_x = 0.0f;
		this->v_y = 0.0f;;
	};

	Node(const datatype x_in, const datatype y_in, const int x_index_in, const int y_index_in,
	 const datatype p_in, const datatype v_x_in, const datatype v_y_in){
		this->x = x_in;
		this->y = y_in;

		this->x_index = x_index_in;
		this->y_index = y_index_in;

		this->p = p_in;
		this->v_x = v_x_in;
		this->v_y = v_y_in;
	};

	//Addition operator
	Node operator+(const Node& rhs){
		return Node(this->x+rhs.x, this->y+rhs.y, this->x_index+rhs.x_index, this->y_index+rhs.y_index, this->p+rhs.p, this->v_x+rhs.v_x, this->v_y+rhs.v_y);
	};

	//Division with a number
	Node operator/(const datatype& rhs){
		return Node(this->x/rhs, this->y/rhs, this->x_index/rhs, this->y_index/rhs, this->p/rhs, this->v_x/rhs, this->v_y/rhs);
	};

};