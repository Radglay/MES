#pragma once
#include "N.h"
#include "Node.h"
#include "Border.h"
#include "Data.h"
#include <string>

class Wall {
public:	
	
	double*** Hbcp;
	double** Hbc;
	double** N_wall;
	double detJ;
	double** pc;
	
	double*** Ppc;
	double** P;
	Border border;
	N n;

	Wall();
	Wall(N, Node*, Border);

	void calculateN();
	void calculateIntegrationPoints();
	void calculateHbc();
	void calculateDetJ(Node*);
	void calculateP();

	std::string enumToString();
	void printN();
	void printHbc();
};