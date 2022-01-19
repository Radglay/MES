#pragma once
#include "N.h"
#include "Element.h"


class Element1D : public Element {
int size;

//double* derivativeNksi;


public:
	Element1D(N);

	double integration1D( N, double, double);
	double interpolation1D(double[], double, double);
	double function1D(double);
};

