#include "Wall.h"
#include <iostream>

Wall::Wall() {}

Wall::Wall(N n, Node* nodes, Border border) {
	this->border = border;
	this->n = n;

	N_wall = new double* [n.size];
	for (int i = 0; i < n.size; i++) {
		N_wall[i] = new double[4];

		memset(N_wall[i], 0, 4 * sizeof(double));
	}

	pc = new double* [n.size]; 
	for (int i = 0; i < n.size; i++) {
		pc[i] = new double[2];

		memset(pc[i], 0, 2 * sizeof(double));
	}

	Hbcp = new double** [n.size];
	for (int i = 0; i < n.size; i++) {
		Hbcp[i] = new double*[4];
		for (int j = 0; j < 4; j++) {
			Hbcp[i][j] = new double[4];

			memset(Hbcp[i][j], 0, 4 * sizeof(double));
		}
	}

	Hbc = new double* [4];
	for (int i = 0; i < 4; i++) {
		Hbc[i] = new double[4];

		memset(Hbc[i], 0, 4 * sizeof(double));
	}
	

	Ppc = new double** [n.size];
	for (int i = 0; i < n.size; i++) {
		Ppc[i] = new double* [4];
		for (int j = 0; j < 4; j++) {
			Ppc[i][j] = new double[1]{ 0 };
		}
	}

	P = new double* [4];
	for (int i = 0; i < 4; i++) {
		P[i] = new double[1]{ 0 };
	}

	calculateDetJ(nodes);
	//integration points initialization
	calculateIntegrationPoints();
	calculateN();

	calculateHbc();
	calculateP();

}


void Wall::calculateIntegrationPoints() {
	switch (border) {
	case LEFT:
		for (int i = 0; i < n.size; i++) {
			pc[i][0] = -1;
			pc[i][1] = -1 * n.nodes[i];
		}
		break;
	case DOWN:
		for (int i = 0; i < n.size; i++) {
			pc[i][0] = n.nodes[i];
			pc[i][1] = -1;
		}
		break;
	case RIGHT:
		for (int i = 0; i < n.size; i++) {
			pc[i][0] = 1;
			pc[i][1] = -1 * n.nodes[i];
		}
		break;
	case NORTH:
		for (int i = 0; i < n.size; i++) {
			pc[i][0] = n.nodes[i];
			pc[i][1] = 1;
		}
		break;
	default:
		break;
	}
}

void Wall::calculateDetJ(Node* nodes) {
	if (nodes[1].BC == 1 && nodes[0].BC == 1) {
		if (border == LEFT || border == RIGHT) {
			this->detJ = abs(nodes[1].y - nodes[0].y) / (1 - -1);
		}
		else {
			this->detJ = abs(nodes[1].x - nodes[0].x) / (1 - -1);
		}
	}
	else {
		this->detJ = 0;
	}
}

void Wall::calculateN() {

	for (int i = 0; i < n.size; i++) {
		N_wall[i][0] = (border == LEFT || border == DOWN)? 0.25 * (1 - pc[i][0]) * (1 - pc[i][1]) : 0;
		N_wall[i][1] = (border == DOWN || border == RIGHT)? 0.25 * (1 + pc[i][0]) * (1 - pc[i][1]) : 0;
		N_wall[i][2] = (border == RIGHT || border == NORTH)? 0.25 * (1 + pc[i][0]) * (1 + pc[i][1]) : 0;
		N_wall[i][3] = (border == NORTH || border == LEFT)? 0.25 * (1 - pc[i][0]) * (1 + pc[i][1]) : 0;
	}
}

void Wall::calculateHbc() {

	for (int pc = 0; pc < n.size; pc++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Hbcp[pc][i][j] += n.weight[pc]  * N_wall[pc][i] * N_wall[pc][j] * data::alpha * detJ;
				if (Hbcp[pc][i][j] < 0.00000000001) {
					Hbcp[pc][i][j] = 0;
				}

				Hbc[i][j] += Hbcp[pc][i][j];
			}
		}
	}
}


void Wall::calculateP() {

	for (int pc = 0; pc < n.size; pc++) {
		for (int i = 0; i < 4; i++) {
			Ppc[pc][i][0] += data::alpha * n.weight[pc] * N_wall[pc][i] * data::T_ot * detJ;
			
			if (Ppc[pc][i][0] < 0.0000000001) {
				Ppc[pc][i][0] = 0;
			}

			P[i][0] += Ppc[pc][i][0];
		}
	}
}


void Wall::printHbc() {
	std::cout << "wall <" << enumToString() << ">" << std::endl;
	std::cout << "Hbc: " << std::endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << Hbc[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Wall::printN() {
	std::cout << "wall <" << enumToString() << ">" << std::endl;
	std::cout << "detJ: "  << detJ << std::endl;
	std::cout << "N matrix "<< std::endl;
	for (int i = 0; i < n.size; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << N_wall[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

std::string Wall::enumToString() {
	switch (border) {
		case LEFT:
			return "LEFT";
			break;
		case RIGHT:
			return "RIGHT";
			break;
		case NORTH:
			return "NORTH";
			break;
		case DOWN:
			return "DOWN";
			break;
		default:
			break;
	}

	return "";
}