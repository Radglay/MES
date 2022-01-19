#include <iostream>
#include "grid.h"
#include "Element2D.h"


Grid::Grid() {}

Grid::Grid(N n) {
   // if(H <= 0 || B <= 0 || nH <= 0 || nB <= 0) {
   //     std::cout << "Wrong values passed to the constructor!" << std::endl;
   // }

    this->B = data::B;
    this->H = data::H;
    this->nH = data::nH;
    this->nB = data::nB;
    this->nE = (nB - 1) * (nH - 1);
    this->nN = nB * nH;

    //arrays initialization
    elements = new Element[nE];
    nodes = new Node[nN];
    ///

    double dx = B / ((double)nB - 1);
    double dy = H / ((double)nH - 1);

 
    int row = 0;
    int column = 0;

    for (int i = 0; i < nN; i++) {
        //node initialization

        row = i / nH;
        column = (i % nH);
        


        std::cout << row << " " << column << std::endl;

        nodes[i].x = row * dx;
        nodes[i].y = column * dy;

        //temp T0
        nodes[i].T = data::T0;

        if (row == 0 || column == 0) {
            nodes[i].BC = 1;
        }
        else if (row == (nB - 1) || column == (nH - 1)) {
            nodes[i].BC = 1;
        }
        else
        {
            nodes[i].BC = 0;
        }

    }

    int k = 1;
    for (int i = 0; i < nN; i++) {
        if (i < nE) { // element initialization
            int ID[4];
          
          //  elements[i].ID = new int[4];

            if (k % nH == 0) { //new column
                k++;
            }

            //std::cout << k << std::endl;

            ID[0] = k;
            ID[1] = (k) + nH;
            ID[2] = (k + nH) + 1;
            ID[3] = (k) + 1;
            k++;

         
            elements[i] = Element2D(n,
                new Node[4]{nodes[ID[0] - 1], nodes[ID[1] - 1], nodes[ID[2] - 1], nodes[ID[3] - 1]},
                new int[4]{ ID[0], ID[1], ID[2], ID[3] });
        }
    }





   //global H matrix
   H_global = new double* [nN];
   for (int i = 0; i < nN; i++) {
       H_global[i] = new double [nN];

       memset(H_global[i], 0, nN * sizeof(double));
   }

   Hbc_global = new double* [nN];
   for (int i = 0; i < nN; i++) {
       Hbc_global[i] = new double[nN];

       memset(Hbc_global[i], 0, nN * sizeof(double));
   }

   P_global = new double* [nN];
   for (int i = 0; i < nN; i++) {
       P_global[i] = new double[nN];

       memset(P_global[i], 0, nN * sizeof(double));
   }


   C_global = new double* [nN];
   for (int i = 0; i < nN; i++) {
       C_global[i] = new double[nN];

       memset(C_global[i], 0, nN * sizeof(double));
   }

}


///////

void Grid::calculateH_global() {
    for (int e = 0; e < nE; e++) {
        int* ID = elements[e].ID;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                H_global[ID[i] - 1][ID[j] - 1] += elements[e].H[i][j];
                //std::cout << i << " " << j << std::endl;
            }
        }
    }
}

void Grid::calculateHbc_global() {

    for (int e = 0; e < nE; e++) {
        int* ID = elements[e].ID;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Hbc_global[ID[i] - 1][ID[j] - 1] += elements[e].Hbc[i][j];
               
                H_global[ID[i] - 1][ID[j] - 1] += elements[e].Hbc[i][j];
                //std::cout << i << " " << j << std::endl;
            }
        }
    }
}



void Grid::calculateP_global() {
    for (int e = 0; e < nE; e++) {
        int* ID = elements[e].ID;

        for (int i = 0; i < 4; i++) {
            //P_global[ID[i] - 1][0] += elements[e].walls[w].P[i][0];
            P_global[ID[i] - 1][0] += elements[e].P[i][0];
        }
    }
}

void Grid::calculateC_global() {
    for (int e = 0; e < nE; e++) {
        int* ID = elements[e].ID;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                C_global[ID[i] - 1][ID[j] - 1] += elements[e].C[i][j];

            }
        }
    }
}



double* Grid::calculateP_dt() {
    double* P = new double[nN];
    double** C_dt = calculateC_dt();
    for (int i = 0; i < nN; i++) {
        P[i] = 0;
        double Cdt = 0;
      
        for (int j = 0; j < nN; j++) {
           Cdt += C_dt[i][j] * nodes[j].T;
        }
        P[i] = P_global[i][0] + Cdt;
    }

    return P;
}


//??
double** Grid::calculateC_dt() {  //[C]/dt
    double** C = new double* [nN];
    for (int i = 0; i < nN; i++) {
        C[i] = new double[nN];

        for (int j = 0; j < nN; j++) {
            C[i][j] = C_global[i][j] / data::dt;
//            H_global[i][j] += C_global[i][j];
        }
    }

    return C;
}




void Grid::calculateTemp() {
    double* x = new double[nN];


    //matrix initialization
    double** Gauss = new double* [nN];
    for (int i = 0; i < nN; i++) {
        Gauss[i] = new double[nN + 1];
    }


    double** C_dt = calculateC_dt();
    double* P_dt = calculateP_dt();

    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            Gauss[i][j] = H_global[i][j] + C_dt[i][j];
        }
        Gauss[i][nN] = P_dt[i]; //
            //nodes[i].T += (C_global[i][j]) / (H_global[i][j] + C_global[i][j]);
    }




    //calculation
    int i = 1;
    while (i < nN) {
        for (int j = 0; j < i; j++) {
            if (Gauss[i][j] == 0) { //nothing to do
                continue;
            }

            double factor = Gauss[i][j] / Gauss[j][j];
            for (int k = 0; k < (nN + 1); k++) {
                double value_to_substract = 0;

                if (Gauss[j][k] != 0) {
                    value_to_substract = Gauss[j][k] * factor;
                }

                Gauss[i][k] -= value_to_substract;
            }
        }

        i++;
    }



    //calculation of {t1}
    double* temp = new double[nN];

    for (int i = nN - 1; i >= 0; i--) {
        temp[i] = 1;

        double value = Gauss[i][nN];
        for (int j = nN - 1; j >= i; j--) {
            if (i != j) {
                value -= temp[j] * Gauss[i][j];
            }
            else {
                temp[i] = value / Gauss[i][j];
            }
        }
    }

    

    //writing results to the {t1} 
    for (int i = 0; i < nN; i++) {
        nodes[i].T = temp[i];
    }


    data::iteration++;
}


void Grid::printTemp() {

    //printing
    std::cout << "iteration <" << data::iteration << ">" << std::endl;
    for (int i = 0; i < nN; i++) {
        std::cout << nodes[i].T << std::endl;
    }

    double T_min = INT32_MAX;
    double T_max = INT32_MIN;

    for (int i = 0; i < nN; i++) {
        if (nodes[i].T < T_min) {
            T_min = nodes[i].T;
        }
        else {
            T_max = nodes[i].T;
        }
    }

    std::cout << "Time: " << data::iteration * data::dt << "; T_min:  " <<  T_min << "; T_max : " << T_max << std::endl;
    std::cout << std::endl;
}


void Grid::printT_max_min() {
    double T_min = INT32_MAX;
    double T_max = INT32_MIN;

    for (int i = 0; i < nN; i++) {
        if (nodes[i].T < T_min) {
            T_min = nodes[i].T;
        }
        else {
            T_max = nodes[i].T;
        }
    }

    std::cout << "Time: " << data::iteration * data::dt << "; T_min:  " << T_min << "; T_max : " << T_max << std::endl;
    std::cout << std::endl;
}



void Grid::printElements() {

    for(int e = 0; e < nE; e++) {
            std::cout << elements[e].ID[3] << " " << elements[e].ID[2] << "\n"
                << elements[e].ID[0] << " " << elements[e].ID[1] << std::endl;
    }
}

void Grid::printNodes() {
    for (int i = 0; i < nN; i++) {
        std::cout << "(" << i + 1 << ")" << nodes[i].x << " " << nodes[i].y << " BC: " << nodes[i].BC << std::endl;
    }
}


void Grid::printH() {
    std::cout << "H global matrix: " << std::endl;
    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            std::cout << H_global[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Grid::printHbc() {
    std::cout << "Hbc global matrix: " << std::endl;
    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            std::cout << Hbc_global[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Grid::printP() {
    std::cout << "P global matrix: " << std::endl;
    for (int i = 0; i < nN; i++) {
      
        std::cout << P_global[i][0] << " ";
      
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Grid::printC() {
    std::cout << "C global matrix: " << std::endl;
    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            std::cout << C_global[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}