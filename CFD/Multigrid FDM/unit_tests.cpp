#define UNITTESTING
#include "adaptive_multigrid.cpp"



void test_node(){
	std::cout<<"Testing Node"<<std::endl;

	Node n1;
	n1.x_index = 5;
	Node n2(0.4, 0.4, 1, 1,
	 0.3, 0.2, 0.1);
	n1.nodeLeft = &n2;
	n1 = n1+n2;
	assert(n1.nodeLeft = &n2);

	assert(n1.stream - 0.19999 <0.01);
	assert(n1.x_index == 5);
	assert(n1.y_index == -1);

	assert(n1.x_index_global == n2.x_index_global);
	assert(n1.x_index_global == -1);

	assert(n2.nodeLeft == NULL);

	std::cout<<"	Node seems to work."<<std::endl;
}

void test_adaptive_grid(){

	std::cout<<"Testing AdaptiveGrid"<<std::endl;

	Node *array;
	array = (Node*) malloc(7*sizeof(Node));

	array[0] = Node(0.0f, 0.0f, 0,0, 1.0f, 1.0f, 1.0f);
	array[1] = Node(0.0f, 1.0f, 0,2, 1.0f, 1.0f, 1.0f);
	array[2] = Node(1.0f, 0.0f, 2,0, 1.0f, 1.0f, 1.0f);
	array[3] = Node(1.0f, 1.0f, 2,2, 1.0f, 1.0f, 1.0f);
	array[4] = Node(0.5f, 0.5f, 1,1, 1.0f, 1.0f, 1.0f);

	array[5] = Node(0.0f, 1.0f, 1,3, 1.0f, 2.0f, 1.0f);
	array[6] = Node(0.0f, 1.0f, 1,1, 1.0f, 2.0f, 1.0f);


	AdaptiveGrid grid2(1,2,0,0, &array[0], 5, 0.5);
	AdaptiveGrid grid1(2,2,0,0, &array[5], 2, 0.25);
	assert(grid2.len == 9);

	assert(grid1.len = 15);

	grid1.coarserGrid = &grid2;
	grid1.updateBoundary();

	assert(grid1.findNode(1,1)->nodeAbove->y_index == 2);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			assert(grid1.isInBoundary(i,j) == false);
			assert(grid1.findNode(i,j)->y_index == j);
		}
	}

	assert(grid1.getFromBoundary(3,1) -0.99999f < 0.001f);
	free(array);
	std::cout<<"	AdaptiveGrid seems to work."<<std::endl;

}

void test_smoother(){

	std::cout<<"Testing Smoother"<<std::endl;

	Node *array;

	array = (Node*) malloc(17*17*sizeof(Node));

	for (int x = 0; x < 17; ++x)
	{
		for (int y = 0; y < 17; ++y)
		{
			array[x*17+y] = Node(x*1.0/16.0, y*1.0/16.0, x,y, 1.0f, 1.0f, 1.0f);
		}
	}

	Node *array2;
	array2 = (Node*) malloc(5*sizeof(Node));

	array2[0] = Node(0.0f, 0.0f, 0,0, 1.0f, 1.0f, 1.0f);
	array2[1] = Node(0.0f, 1.0f, 0,2, 1.0f, 1.0f, 1.0f);
	array2[2] = Node(1.0f, 0.0f, 2,0, 1.0f, 1.0f, 1.0f);
	array2[3] = Node(1.0f, 1.0f, 2,2, 1.0f, 1.0f, 1.0f);
	array2[4] = Node(0.5f, 0.5f, 1,1, 1.0f, 1.0f, 1.0f);

	AdaptiveGrid grid2(1,2,0,0, &array2[0], 5, 0.5);

	AdaptiveGrid grid(2,2,0,0, &array[0], 17*17, 1.0/16.0);
	grid.coarserGrid = &grid2;

	grid.updateBoundary();

	for (int i = 0; i < 500; ++i)
	{
		grid.jacobiSmootherLaplacianStream();
	}
	
	double error = 0.0;
	for (int i = 0; i < 17*17; ++i)
	{
		error += grid.u[i].stream*grid.u[i].stream;
	}

	std::cout<<"	error: "<<error<<std::endl;
	if(error<0.001){
		std::cout<<"	Smoother seems to work."<<std::endl;
	}
	else{
		assert(error<0.001);
	}

	free(array);
	free(array2);

}






datatype my_func(datatype x, datatype y){
	if(x>=0.25 && x<=0.75 && y>=0.25 && y<=0.75){
		return 10.0f;
	}
	else{
		return 0.0f;
	}
}

void test_adaptiveMultigrid(){

	std::cout<<"Testing Adaptive Multigrid."<<std::endl;

	datatype (*my_func_ptr)(datatype, datatype) = my_func;



	Node *array;

	array = (Node*) malloc((5+4+16)*sizeof(Node));

	array[0] = Node(0.0f, 1.0f, 0,2, 1.0f, 1.0f, 1.0f);
	array[1] = Node(1.0f, 0.0f, 2,0, 1.0f, 1.0f, 1.0f);
	array[2] = Node(0.0f, 0.0f, 0,0, 1.0f, 1.0f, 1.0f);
	array[3] = Node(1.0f, 1.0f, 2,2, 1.0f, 1.0f, 1.0f);
	array[4] = Node(0.5f, 0.5f, 1,1, 1.0f, 1.0f, 1.0f);

	float h = 0.25f;
	for (int x = 0; x < 2; ++x)
	{
		for (int y = 0; y < 2; ++y)
		{
			array[5 + x*2+y] = Node( (2*x+1)*h, (2*y+1)*h, 2*x+1, 2*y+1, 1.0f, 1.0f, 1.0f);
		}
	}
	h = h/2.0f;

	for (int x = 0; x < 2; ++x)
	{
		for (int y = 0; y < 2; ++y)
		{
			array[5+4 + x*2+y] = Node( (2*x+1)*h+0.25, (2*y+1)*h+0.25, 2*x+1, 2*y+1, 1.0f, 1.0f, 1.0f);
		}
	}

	AdaptiveGrid grid1(1,3,0,0, &array[0], 5, 0.5);
	AdaptiveGrid grid2(2,3,0,0, &array[5], 4, 0.25);
	AdaptiveGrid grid3(3,3,2,2, &array[5+4], 4, 0.125);

	grid2.coarserGrid = &grid1;
	grid3.coarserGrid = &grid2;

	grid1.finerGrid = &grid2;
	grid2.finerGrid = &grid3;

	grid1.updateBoundaryCoarsestGrid();
	grid1.updateBFromFunction(my_func_ptr);

	grid2.updateBoundary();
	//grid2.updateBFromFunction(my_func_ptr);

	grid3.updateBoundary();
	//grid3.updateBFromFunction(my_func_ptr);


	for (int i = 0; i < 100; ++i)
	{
		fullMultigrid(2, &grid1, 50, 100, 50, my_func_ptr);

	}
	
	for (int i = 0; i < grid1.len; ++i)
	{
		std::cout<<"u[i]: "<<grid1.u[i].x_index<<" ; "<<grid1.u[i].y_index<<" , "<<grid1.u[i].stream<<std::endl;
	}
	std::cout<<std::endl;

	for (int i = 0; i < grid2.len; ++i)
	{
		std::cout<<"u[i]: "<<grid2.u[i].x_index<<" ; "<<grid2.u[i].y_index<<" , "<<grid2.u[i].stream<<std::endl;
	}
	std::cout<<std::endl;

	for (int i = 0; i < grid3.len; ++i)
	{
		std::cout<<"u[i]: "<<grid3.u[i].x_index<<" ; "<<grid3.u[i].y_index<<" , "<<grid3.u[i].stream<<std::endl;
	}
	std::cout<<std::endl;

	free(array);
	std::cout<<"	Adaptive Multigrid seems to work."<<std::endl;

}


int main(int argc, char const *argv[])
{
	test_node();
	test_adaptive_grid();
	test_smoother();
	test_adaptiveMultigrid();
	return 0;
} 
