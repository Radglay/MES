#include "Element.h"
#include "Data.h"
#include <iostream>

Element::Element() {}


void Element::printDerivativeNksi() {
    std::cout << "dNdKsi" << std::endl;
	for (int i = 0; i < size * size; i++) {
		std::cout << "(" << i + 1 << ") ";
		for (int j = 0; j < 4; j++) {
			std::cout << derivativeNksi[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Element::printDerivativeNeta() {
    std::cout << "dNdEta" << std::endl;
	for (int i = 0; i < size * size; i++) {
		std::cout << "(" << i + 1 << ") ";
		for (int j = 0; j < 4; j++) {
			std::cout << derivativeNeta[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}



void Element::printDerivativeNx() {
    std::cout << "dNdX" << std::endl;
    for (int i = 0; i < size * size; i++) {
        std::cout << "(" << i + 1 << ")" << " ";
        for (int j = 0; j < 4; j++) {
            std::cout << derivativeNx[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Element::printDerivativeNy() {
    std::cout << "dNdY" << std::endl;
    for (int i = 0; i < size * size; i++) {
        std::cout << "(" << i + 1 << ")" << " ";
        for (int j = 0; j < 4; j++) {
            std::cout << derivativeNy[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Element::printInvJacobian(int point) {
    std::cout << "Inv Jakobian" << std::endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            std::cout << J[point].inv_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
 }



void Element::calculateJacobian(int i) {
    Jacobian2D jacobian;


    //1 row
    double dXdKsi = 0.0;
    double dYdKsi = 0.0;


    //2 row
    double dXdEta = 0.0;
    double dYdEta = 0.0;


    //i-element x values
    double X[4] = {
        nodes[0].x,
        nodes[1].x,
        nodes[2].x,
        nodes[3].x
    };

    //i-element y values
    double Y[4] = {
        nodes[0].y,
        nodes[1].y,
        nodes[2].y,
        nodes[3].y
    };



    //i-integration point
    for (int j = 0; j < 4; j++) {
        dXdKsi += this->derivativeNksi[i][j] * X[j];
        dXdEta += this->derivativeNeta[i][j] * X[j];
        dYdKsi += this->derivativeNksi[i][j] * Y[j];
        dYdEta += this->derivativeNeta[i][j] * Y[j];
    }

    //testing calculation errors
    if (dXdKsi < 0.000000001)
        dXdKsi = 0;

    if (dXdEta < 0.000000001)
        dXdEta = 0;

    if (dYdKsi < 0.000000001)
        dYdKsi = 0;

    if (dYdEta < 0.000000001)
        dYdEta = 0;

    ///

    jacobian.matrixJ[0][0] = dXdKsi;
    jacobian.matrixJ[0][1] = dYdKsi;
    jacobian.matrixJ[1][0] = dXdEta;
    jacobian.matrixJ[1][1] = dYdEta;

    jacobian.detJ = jacobian.matrixJ[0][0] * jacobian.matrixJ[1][1] - jacobian.matrixJ[0][1] * jacobian.matrixJ[1][0];

    //inverted matrix

    jacobian.inv_matrix[0][0] = 1.0 / jacobian.detJ * jacobian.matrixJ[1][1];
    jacobian.inv_matrix[0][1] = -1.0 / jacobian.detJ * jacobian.matrixJ[0][1];
    jacobian.inv_matrix[1][0] = -1.0 / jacobian.detJ * jacobian.matrixJ[1][0];
    jacobian.inv_matrix[1][1] = 1.0 / jacobian.detJ * jacobian.matrixJ[0][0];

    J[i] = jacobian;
}



void Element::calculateDxDy() {
    double integral = 0.0;

  //  int size = this->size * this->size;
    for (int i = 0; i < size * size; i++) { //integration point number


        //for each integration point
        double dNdX;
        double dNdY;
        for (int j = 0; j < 4; j++) { //shape function number

            //test, integration point -> 0
            derivativeNx[i][j] = J[i].inv_matrix[0][0] * derivativeNksi[i][j] + J[i].inv_matrix[0][1] * derivativeNeta[i][j];
            derivativeNy[i][j] = J[i].inv_matrix[1][0] * derivativeNksi[i][j] + J[i].inv_matrix[1][1] * derivativeNeta[i][j];
            ///
        }
    }
}


  

void Element::calculateH(N n) {
    //each integration point
 
    for(int i = 0; i < size *  size; i++) {
        int row = i / size; //which row
        int column = i % size; //which column

        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {

                //dx
                H_p[i][j][k] = derivativeNx[i][j] * derivativeNx[i][k];

                //dy
                H_p[i][j][k] += derivativeNy[i][j] * derivativeNy[i][k];


                H_p[i][j][k] = H_p[i][j][k] * data::k * J[i].detJ * n.weight[row] * n.weight[column];//wagi? 
            }
        }
    }


    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            H[i][j] = 0;
            for (int k = 0; k < size * size; k++) {
                H[i][j] += H_p[k][i][j];
            }
        }
    }

}



void Element::calculateHbc() {
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
            for (int w = 0; w < 4; w++) {
                Hbc[i][j] += walls[w].Hbc[i][j];
            }
        }
    }
}


void Element::calculateP() {
    for (int i = 0; i < 4; i++) {
        for (int w = 0; w < 4; w++) {
            //if (walls[w].P[i][0] > 0) {
                P[i][0] += walls[w].P[i][0];
           // }
        }
    }
}


void Element::calculateC(N n) {
    for (int i = 0; i < size * size; i++) {
        int row = i / size; //which row
        int column = i % size; //which column

        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {

                C_p[i][j][k] += N_matrix[i][j] * N_matrix[i][k];


                C_p[i][j][k] = C_p[i][j][k] * data::c * data::ro * J[i].detJ * n.weight[row] * n.weight[column];
            }
        }
    }


    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            C[i][j] = 0;
            for (int k = 0; k < size * size; k++) {
                C[i][j] += C_p[k][i][j];
            }
        }
    }
}





void Element::printH_p(int point) {
        std::cout << "H for integration point: " << point << std::endl;
        point = point - 1;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                std::cout << H_p[point][j][k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
}


void Element::printH() {

    std::cout << "H matrix: " << std::endl;

    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            std::cout << H[j][k] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


void Element::printC_p(int point) {
    std::cout << "C for integration point: " << point << std::endl;
    point = point - 1;
    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            std::cout << C_p[point][j][k] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Element::printC() {
    std::cout << "C matrix: " << std::endl;

    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            std::cout << C[j][k] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Element::print() {
    printDerivativeNksi();

    printDerivativeNeta();

    printDerivativeNx();

    printDerivativeNy();


    printInvJacobian(0);

    printH_p(1);
    printH_p(2);
    printH_p(3);
    printH_p(4);

    printH();
    
}

    /*
    double Element::interpolation2D(double N[], Grid grid) {
        return 0;
    }


*/