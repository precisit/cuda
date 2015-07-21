#include "node.cpp"

class Grid
{
    /*
        Contains a structured uniform 2D grid of size len-by-len.

        This should (and to some extent already do) contain all
        data necessary for the multigrid function. That includes 
        stuff like the differential operator and B.C and the ptrs
        to the correct nodes.

        And also: STAGGERD GRIDS SUCK AND WE DO IT OLD SCHOOL!
    */
 
    public:
        Node* ptr;
        int len;
        datatype h;
 
        void print();
 
        void restrict(Grid* coarseGrid);
 
        void interpolate(Grid* fineGrid, const int fineN);
 
        void setValues(Node* ptrIn, int n);
 
        datatype getU(const int index);
        void setU(const datatype u, const int index);
 
        datatype getLaplacian(const int i, const int j);
        datatype getVortTranspDisc(const int i, const int j);
        datatype getB(const int i);
 
        int lengthOfFinerGrid();
 
        int lengthOfCoarserGrid();
 
        void jacobiSmootherLaplacian();
        void jacobiSmootherVortTranspDisc();

        void calculateError();
 
        //Constructors
        Grid(Node* nodeIn, const int nIn);
        Grid(const int n);
        Grid();
 
        //Destructor
        ~Grid();
};