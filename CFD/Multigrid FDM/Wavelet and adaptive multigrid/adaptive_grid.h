#ifndef ADAPTIVE_GRID_H
#define ADAPTIVE_GRID_H

#include <assert.h>
#include "node.h"
//#include "parameters.h"

typedef datatype (*func_t)(datatype, datatype);
typedef datatype (*func_)(int, int);

class AdaptiveGrid{
public:

	Node* u;
    int len;
    Node* b;
    Node* d;
    Node* w;

    Node* boundaryVals;
    int* boundaryIndex;
    int lenOfBoundary;

    Node** pointsChosenByWavelet; //An array of Node pointers
    int numOfPointsChosen;

    int layerNr;
    int numOfLayers;
    int origo_x;
    int origo_y;

    AdaptiveGrid* finerGrid;// = NULL;
    AdaptiveGrid* coarserGrid;// = NULL;
    datatype h;

	bool isInBoundary(const int x, const int y);

	void setB(const datatype val, const int i);

	datatype getB(const int i);

	datatype getFromBoundary(const int x, const int y);

	void calculateErrorLaplacian();

	datatype getLaplacianStream(const int index1, const int index2);

	void resetBoundaryLength();

	bool isInGrid(const int x, const int y);

	void setBoundaryLength();

	void setU(const datatype val, const int i);

	datatype getU(const int i);

	void jacobiSmootherLaplacianStream();

	void restrictUtoU(AdaptiveGrid* coarse);

	void restrictDtoB(AdaptiveGrid* coarse);

	void restrictDtoD(AdaptiveGrid* coarse);

	Node * findNode(const int ind_x, const int ind_y);

	Node * findNodeD(const int ind_x, const int ind_y);

	void interpolateU(AdaptiveGrid *fine);


	void interpolateD(AdaptiveGrid *fine);

	void findNeighboursD();

	void findNeighboursU();


	void setNeighbours();

	void local2global(const int* x_loc, const int* y_loc, int* x_glo, int* y_glo );


	void setupGrid(Node* savedNodes, const int numberOfPoints);

	void setupCoarsestGrid(Node* savedNodes, const int numberOfPoints);

	Node* findNodeGeneral(Node* arr, const int ind_x, const int ind_y);


	void findNeighbours(Node * arr);

	AdaptiveGrid();

	AdaptiveGrid(const int layerIn, const int numLayersIn, const int origo_x_in, const int origo_y_in, 
		Node* savedNodesIn , const int numberOfPointsIn, const datatype h_in);

	Node interpolateGhostPointFromGlobal(int x_glo, int y_glo);

	Node findNodeFromGlobalIndex(const int x, const int y);

	Node interpolateGhostPoint(int x_loc, int y_loc);

	void setBoundary();

	void updateBFromBoundary();

	void updateBFromFunction( func_t externalForceFunction);

	void calculateRHS();

};

#endif