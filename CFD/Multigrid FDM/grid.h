#include "node.cpp"

class Grid
{
    /*
        Contains a structured uniform 2D grid of size len-by-len.

        This should (and to some extent already do) contain all
        data necessary for the multigrid function. That includes 
        stuff like the differential operator and B.C and the ptrs
        to the correct nodes.

        And also: STAGGERED GRIDS SUCK AND WE DO IT OLD SCHOOL!
    */
 
    public:
        Node* ptr;
        int len;
        Node* B;

        //This is the error on the grid. It dumb to store it here. 
        //It's dumb to store it at all. OPT!
        Node* D;

        datatype h;

        Grid* finerGrid;
        Grid* coarserGrid;
 
        void print();
        void printB();
        void printD();
        void printVelx();
        void printVely();
 
        void restrict(Grid* coarseGrid);
        void restrictDtoB(Grid* coarseGrid);
 
        void interpolate(Grid* fineGrid, const int fineN);
        void interpolateB(Grid* fineGrid, const int fineN);
 
        void setValues(Node* ptrIn, int n);
 
        datatype getU(const int index);
        void setU(const datatype u, const int index);
 
        datatype getLaplacian(const int i, const int j);
        datatype getVortTranspDisc(const int i, const int j);

        datatype getB(const int i);
        void setB(const datatype u, const int index);

        datatype getD(const int i);
        void setD(const datatype u, const int index);
 
        int lengthOfFinerGrid();
 
        int lengthOfCoarserGrid();
 
        void jacobiSmootherLaplacian();
        void jacobiSmootherVortTranspDisc();

        void calculateErrorLaplacian();
 
        //Constructors
        Grid(Node* nodeIn, const int nIn);
        Grid(const int n);
        Grid();
 
        //Destructor
        ~Grid();

        bool isInternal(const int x, const int y);

        void addSomeStuffTogether();
};