#pragma once
#include "Jacobian.h"
#include "Wall.h"


class Element {
public:

    double** N_matrix;

    double** derivativeNksi;
    double** derivativeNeta;

    double*** H_p;
    double** H;
    double** Hbc;
    double** P;
    Wall* walls;
    Node* nodes;
    //œciany

    double*** C_p;
    double** C;


    double** derivativeNx;
    double** derivativeNy;
    Jacobian2D* J;

    //
    int* ID;
    int size;

    Element();


    void printDerivativeNksi();
	void printDerivativeNeta();



	void printDerivativeNx();
	void printDerivativeNy();


    void calculateDxDy();
    void calculateJacobian(int);


	void printInvJacobian(int);


    void calculateH(N n);  

    
    void printH_p(int);
    void printH();

    void printC();
    void printC_p(int);
    void calculateHbc();

    void calculateP();

    void calculateC(N n);

    void print();


};