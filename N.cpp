#include "N.h"




double function1D(double x) {
	return 5 * x * x + 3 * x + 6;
}


double function2D(double x, double y) {
	return 5 * x * x * y * y + 3 * x * y + 6;
}

double integration1D(N n) {
	double result = 0.0;

	for (int i = 0; i < n.size; i++) {
		result += n.weight[i] * function1D(n.nodes[i]);
	}

	return result;
}

double integration2D(N n) {
	double result = 0.0;

	for (int i = 0; i < n.size; i++) {
		for (int j = 0; j < n.size; j++) {
			result += n.weight[i] * n.weight[j] * function2D(n.nodes[j], n.nodes[i]);
		}
	}

	return result;
}
