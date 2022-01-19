#pragma once
#include <cmath>


class N {
public:
	double* nodes;
	double* weight;
	int size;

	friend double function1D(double);
	friend double function2D(double, double);
	friend double integration1D(N n);
	friend double integration2D(N n);
};

class N2 : public N { //2 punkty ca³kowania
	double nodesValues[2] = { -1.0 / sqrt(3), 1.0 / sqrt(3) };
	double weightValues[2] = { 1.0, 1.0 };

public:
	N2() {
		this->nodes = nodesValues;
		this->weight = weightValues;
		this->size = 2;
	}

};

class N3 : public N { //3 punkty ca³kowania
	double nodesValues[3] = { -sqrt(3.0 / 5.0),0.0 , sqrt(3.0 / 5.0) };
	double weightValues[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

public:
	N3() {
		this->nodes = nodesValues;
		this->weight = weightValues;
		this->size = 3;
	}
};


class N4 : public N { //4 punkty ca³kowania
	double nodesValues[4] = { -0.861136, -0.339981, 0.339981, 0.861136 };
	double weightValues[4] = { 0.347855, 0.652145, 0.652145 ,0.347855 };

public:
	N4() {
		this->nodes = nodesValues;
		this->weight = weightValues;
		this->size = 4;
	}
};

class N5 : public N { //5 punktów ca³kowania
	double nodesValues[5] = { -0.906180, -0.538469, 0.0, 0.538469, 0.906180 };
	double weightValues[5] = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };

public:
	N5() {
		this->nodes = nodesValues;
		this->weight = weightValues;
		this->size = 5;
	}
};

