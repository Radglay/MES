#pragma once

class Node {
public:
    //friend class Jacobian2D;

    bool BC;
    double x;
    double y;

    double T;

    Node();
};