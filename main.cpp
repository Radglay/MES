#include "grid.h"
#include "N.h"

#include <iostream>
#include <fstream>
#include <string>


void export_to_file(std::string file_name, Grid grid) {
    std::ofstream myfile;

    std::string directory = "Test_result/";
    file_name = directory + file_name;

    myfile.open(file_name);

    if (myfile.is_open()) {
        myfile << "SimulationTime " << data::t << "\n";
        myfile << "SimulationStepTime " << data::dt << "\n";
        myfile << "Conductivity " << data::k << "\n";
        myfile << "Alfa " << data::alpha << "\n";
        myfile << "Tot " << data::T_ot << "\n";
        myfile << "InitialTemp " << data::T0 << "\n";
        myfile << "Density " << data::ro << "\n";
        myfile << "SpecificHeat " << data::c << "\n";
        myfile << "Nodes number " << (data::nB) * (data::nH) << "\n";
        myfile << "Elements number " << (data::nB - 1) * (data::nH - 1) << "\n";

        myfile << "*Node\n";
        for (int i = 0; i < (data::nB) * (data::nH); i++) {
            myfile << (i + 1) << ", " << grid.nodes[i].x << ", " << grid.nodes[i].y << "\n";
        }

        myfile << "*Element, type=Element2D4\n";
        for (int i = 0; i < (data::nB - 1) * (data::nH - 1); i++) {
            int* ID = grid.elements[i].ID;

            myfile << (i + 1) << ", " << ID[0] << ", " << ID[1] << ", " << ID[2] << ", " << ID[3] << "\n";
        }

        myfile << "*BC\n";
        for (int i = 0; i < (data::nB) * (data::nH); i++) {
            if (grid.nodes[i].BC == 1) {
                myfile << (i + 1);

                if (i < (data::nB) * (data::nH)-1) {
                    myfile << ", ";
                }
            }
        }


        myfile.close();
    }
    else {
        std::cout << "Problem with writing into the file...\n";
    }
}


void import_from_file(std::string file_name) {
    std::ifstream myfile;

    std::string directory = "Test_import/";
    file_name = directory + file_name;


    std::string line;
    myfile.open(file_name);

    if (myfile.is_open()) {
       
        myfile >> data::T0;
        getline(myfile, line);

        myfile >> data::t;
        getline(myfile, line);

        myfile >> data::dt;
        getline(myfile, line);

        myfile >> data::T_ot;
        getline(myfile, line);

        myfile >> data::alpha;
        getline(myfile, line);

        myfile >> data::H;
        getline(myfile, line);

        myfile >> data::B;
        getline(myfile, line);

        myfile >> data::nH;
        getline(myfile, line);

        myfile >> data::nB;
        getline(myfile, line);

        myfile >> data::c;
        getline(myfile, line);

        myfile >> data::k;
        getline(myfile, line);

        myfile >> data::ro;
        getline(myfile, line);
        myfile.close();
    }
    else {
        std::cout << "Problem with reading from the file...\n";
    }
}

int main() {
    import_from_file("test_case1.txt");
   //import_from_file("test_case2.txt");

    N2 n2;
    N3 n3;
    N4 n4;
    N5 n5;
    Grid grid = Grid(n5);
    //Grid grid = Grid(0.2, 0.1, 5, 4, n2);

   // grid.printElements();

    ///*
   // grid.printNodes();

   // grid.printElements();


    //printing test
    {
        /*
        std::cout << "**************************" << std::endl;
        for (int i = 0; i < grid.nE; i++) {
            std::cout << "Element <" << i + 1 << ">" << std::endl;


            grid.elements[i].printDerivativeNksi();


            grid.elements[i].printDerivativeNeta();


            grid.elements[i].printDerivativeNx();

            grid.elements[i].printDerivativeNy();



            grid.elements[i].printInvJacobian(0);
            // grid.elements[i].printInvJacobian(1);
            // grid.elements[i].printInvJacobian(2);
             //grid.elements[i].printInvJacobian(3);


            std::cout << "**************************" << std::endl;
            */
        }

    /*
        std::cout << "H matrices" << std::endl;

        std::cout << "**************************" << std::endl;
        for (int i = 0; i < grid.nE; i++) {
            std::cout << "Element <" << i + 1 << ">" << std::endl;
            // grid.elements[i].calculateH(25.0);

           //  grid.elements[i].printH_p(1);
            // grid.elements[i].printH_p(2);
           //  grid.elements[i].printH_p(3);
           //  grid.elements[i].printH_p(4);


            grid.elements[i].printH();
            std::cout << "**************************" << std::endl;
        }



        std::cout << "**************************" << std::endl;
        std::cout << "**************************" << std::endl;
        std::cout << "**************************" << std::endl;
        std::cout << "**************************" << std::endl;
        
    }
    */
    grid.elements[0].printDerivativeNx();
    grid.elements[0].printDerivativeNy();

  //  int size = grid.elements[0].size;
///for (int i = 0; i < size * size ; i++) {
   //     grid.elements[0].printH_p(i + 1);
  ///  }
    

    {
        //std::cout <<  grid.elements[0].J->detJ;


       //  grid.printElements();
        // grid.printNodes();


         //walls
         /*

         for (int e = 0; e < grid.nE; e++) {
             std::cout << "Element <" << e + 1 << ">" << std::endl;
             grid.elements[e].walls[0].printN();
             grid.elements[e].walls[1].printN();
             grid.elements[e].walls[2].printN();
             grid.elements[e].walls[3].printN();
             std::cout << std::endl;
         }

         */
         // std::cout << grid.elements[0].walls[3].pc[0][0] << " ";
        //  std::cout << grid.elements[0].walls[3].pc[0][1] << std::endl;

        //  std::cout << grid.elements[0].walls[3].pc[1][0] << " ";
        //  std::cout << grid.elements[0].walls[3].pc[1][1] << std::endl;


          //calculate global HBC


          /*
          for (int e = 0; e < grid.nE; e++) {
              std::cout << "Element <" << e + 1 << ">" << std::endl;
              grid.elements[e].walls[0].printHbc();
              grid.elements[e].walls[1].printHbc();
              grid.elements[e].walls[2].printHbc();
              grid.elements[e].walls[3].printHbc();
              std::cout << std::endl;
          }
          */
          //grid.elements[0].walls[4].printHbc();
    }


    /////////////////////////


    grid.elements[0].printH();
    grid.calculateH_global();
  //  grid.printH();



    //add Hbc to the H
    grid.calculateHbc_global();
   // grid.printHbc();

    
   // std::cout << std::endl << "[H] + [Hbc]" << std::endl;
  //  grid.printH();

    grid.calculateP_global();
   // grid.printP();


    grid.calculateC_global();
   // grid.printC();

  //  grid.calculateC_dt();

    {
        /*
        std::cout << "[C]/dt" << std::endl;
        double** C_dt = grid.calculateC_dt();
        for (int i = 0; i < grid.nN; i++) {
            for (int j = 0; j < grid.nN; j++) {
                std::cout << C_dt[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;


        std::cout << "[H] = [H] + [C]/dt" << std::endl;
        for (int i = 0; i < grid.nN; i++) {
            for (int j = 0; j < grid.nN; j++) {
                std::cout << C_dt[i][j] + grid.H_global[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;


        std::cout << "{P} = {P} + {[C]/dt} * {T0}" << std::endl;
        double* P = grid.calculateP_dt();
        for (int i = 0; i < grid.nN; i++) {
            std::cout << P[i] << std::endl;
        }
        std::cout << std::endl;





        //TEMPERTURES
     ///////////////////////////// 
        grid.calculateTemp();
        grid.printTemp();

        std::cout << "{P} = {P} + {[C]/dt} * {T0}" << std::endl;
        P = grid.calculateP_dt();
        for (int i = 0; i < grid.nN; i++) {
            std::cout << P[i] << std::endl;
        }
        std::cout << std::endl;
*/
    }

    for (int i = 0; i < 20; i++) {
        grid.calculateTemp();
        grid.printT_max_min();
    }
    //500S
    /////////////////////////////////////////





    ///
    export_to_file("Test1_4_4.txt", grid);
    return 0;
}
