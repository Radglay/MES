#include "Element1D.h"
#include <iostream>

Element1D::Element1D(N n) {
	this->size = n.size;
}

double Element1D::integration1D(N n, double a, double b) {
	double integral = 0.0;
	double ksi = 0.0;
	for (int i = 0; i < n.size; i++) { //for each integration point
		ksi = n.nodes[i];
		
		//shape functions
		double N[2] = { 0.5 * (1 - ksi) , 0.5 * (1 + ksi) };


		//interpolation
		double x = interpolation1D(N, a, b);

		integral += n.weight[i] * function1D(x);
	}

	double detJ = (b - a) / 2;

	return integral * detJ;
}

double Element1D::interpolation1D(double N[],  double a, double b) {
	return N[0] * a + N[1] * b;
}

double Element1D::function1D(double x) {
	return 0.5 * x * x + 2 * x + 3;
}