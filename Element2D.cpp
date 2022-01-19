#include "Element2D.h"
#include "Jacobian.h"
#include <iostream>

Element2D::Element2D(){}
Element2D::Element2D(N n, Node* nodes, int* ID) {

	//	if (n.size < 2) {
		//	std::cout << "wrong argument passed to the Element_2D constructor..." << std::endl;
	//		exit(2);
	//	}
	this->nodes = nodes; //element real nodes
	this->size = n.size;	
	this->ID = ID;


	N_matrix = new double* [size * size];
	for (int i = 0; i < size * size; i++) {
		N_matrix[i] = new double[4];

		memset(N_matrix[i], 0, 4 * sizeof(double));
	}


	//array allocation
	derivativeNksi = new double* [size * size];
	for (int i = 0; i < size * size; i++) {
		derivativeNksi[i] = new double[4];

		memset(derivativeNksi[i], 0, 4 * sizeof(double));
	}


	derivativeNeta = new double* [size * size];
	for (int i = 0; i < size * size; i++) {
		derivativeNeta[i] = new double[4];

		memset(derivativeNksi[i], 0, 4 * sizeof(double));
	}




	derivativeNx = new double* [size * size];
	for (int i = 0; i < size * size; i++) {
		derivativeNx[i] = new double[4];

		memset(derivativeNx[i], 0, 4 * sizeof(double));

	}



	derivativeNy = new double* [size * size];
	for (int i = 0; i < size * size; i++) {
		derivativeNy[i] = new double[4];

		memset(derivativeNy[i], 0, 4 * sizeof(double));
	}




	//array initialization
	double eta = 0.0;
	double ksi = 0.0;

	for (int i = 0; i < size * size; i++) {
		int row = i / size; //which row
		int column = i % size; //which column


		ksi = n.nodes[column];
		eta = n.nodes[row];

		//std::cout << "column: " << column << ", row: " << row << std::endl;
		//std::cout << ksi << " " << eta << std::endl << std::endl;

		N_matrix[i][0] = 0.25 * (1 - ksi) * (1 - eta);
		N_matrix[i][1] = 0.25 * (1 + ksi) * (1 - eta);
		N_matrix[i][2] = 0.25 * (1 + ksi) * (1 + eta);
		N_matrix[i][3] = 0.25 * (1 - ksi) * (1 + eta);




		derivativeNksi[i][0] = -0.25 * (1 - eta);
		derivativeNksi[i][1] = 0.25 * (1 - eta);
		derivativeNksi[i][2] = 0.25 * (1 + eta);
		derivativeNksi[i][3] = -0.25 * (1 + eta);



		derivativeNeta[i][0] = -0.25 * (1 - ksi);
		derivativeNeta[i][1] = -0.25 * (1 + ksi);
		derivativeNeta[i][2] = 0.25 * (1 + ksi);
		derivativeNeta[i][3] = 0.25 * (1 - ksi);
	}
	

	//Jacobian initialization
	J = new Jacobian2D[size * size];

	for (int i = 0; i < size * size; i++) {
		calculateJacobian(i);
	}

	calculateDxDy();



	//H_p initialization
	H_p = new double** [size * size];
	for (int i = 0; i < size * size; i++) {
		H_p[i] = new double* [4];
		for (int j = 0; j < 4; j++) {
			H_p[i][j] = new double[4];

			memset(H_p[i][j], 0, 4 * sizeof(double));
		}
	}

	//H initializatison
	H = new double* [4];
	for (int i = 0; i < 4; i++) {
		H[i] = new double[4];
		
		memset(H[i], 0, 4 * sizeof(double));
	}


	Hbc = new double* [4];
	for (int i = 0; i < 4; i++) {
		Hbc[i] = new double[4];

		memset(Hbc[i], 0, 4 * sizeof(double));
	}


	P = new double* [4];
	for (int i = 0; i < 4; i++) {
		P[i] = new double[1]{ 0 };
	}


	C = new double* [4];
	for (int i = 0; i < 4; i++) {
		C[i] = new double[4];

		memset(C[i], 0, 4 * sizeof(double));
	}

	C_p = new double** [size * size];
	for (int i = 0; i < size * size; i++) {
		C_p[i] = new double* [4];
		for (int j = 0; j < 4; j++) {
			C_p[i][j] = new double[4];

			memset(C_p[i][j], 0, 4 * sizeof(double));
		}
	}
	
	//side walls initialization

	N2 n2;
	N3 n3;
	walls = new Wall[4];
	walls[0] = Wall(n2, new Node[2]{ nodes[0], nodes[1] }, Border::DOWN);
	walls[1] = Wall(n2, new Node[2]{ nodes[1], nodes[2] }, Border::RIGHT);
	walls[2] = Wall(n2, new Node[2]{ nodes[3], nodes[2] }, Border::NORTH);
	walls[3] = Wall(n2, new Node[2]{ nodes[0], nodes[3] }, Border::LEFT);


	calculateH(n);
	calculateHbc();
	calculateP();
	calculateC(n);

}