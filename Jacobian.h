#pragma once
#include "Node.h"
//#include "Element2D.h"

class Jacobian2D {
public:
	//friend class Element;

	double detJ;
	double matrixJ[2][2];
	double inv_matrix[2][2];

	Jacobian2D();
	
};

class Jacobian3D {
public:
	Jacobian3D();
};

