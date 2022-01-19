#pragma once
#include "N.h"
#include "Element.h"
#include "Node.h"
#include "Jacobian.h"
#include "Data.h"
//#include "Element2D.h"


class Grid {
public:
   // friend class Jacobian2D;

    double H; //height
    double B; //length

    //vertical number of nodes
    int nB; //horizontal number of nodes
    int nH;

    int nE; //number of elements
    int nN; //number of nodes

    Element* elements; //array of elements
    Node* nodes; //array of nodes

    double** H_global;
    double** Hbc_global;
    double** P_global;
    double** C_global;


    Grid();
    Grid(N); //height, length, nH, nB
    void printElements();
    void printNodes();
    void printH();
    void printHbc();
    void printP();
    void printC();
 
 //   void calculateDxDy(Element&);
    void calculateHbc_global();
    void calculateH_global();
    void calculateP_global();
    void calculateC_global();

    double* calculateP_dt();
    double** calculateC_dt();


    void calculateTemp();
    void printTemp();
    void printT_max_min();

};

