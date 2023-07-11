#include <iostream>
#include <string>

#include "Case.hpp"

void printProgressBar( double t, double t_end ){

    std::string bar;
    int percent = (t/t_end) * 100;

    for(int i = 0; i < 50; i++){
        if( i < (percent/2)){
            bar.replace(i,1,"=");
        }else if( i == (percent/2)){
            bar.replace(i,1,">");
        }else{
            bar.replace(i,1," ");
        }
    }
    std::cout<< "\r" "[" << bar << "] ";
    std::cout.width( 3 );
    std::cout<< percent << "%     \r" << std::flush;
  }

int main(int argn, char **args) { //argn is the no of command line arguments and **args is the array of strings 
    if (argn > 1) {
        std::string file_name{args[1]};
        int method;
        std::cout<<"Select one of the below methods to run your desired case: "<<std::endl;
        std::cout<<"1. Lattice Boltzman Method"<<std::endl;
        std::cout<<"2. Finite Difference Method - Navier Stokes"<<std::endl;
        std::cin>>method;
        Case problem(file_name, argn, args, method); //problem is an object of Case class which is created to further invoke the member functions and member variables of that class like problem.simulate in the next line
        problem.simulate(method);
    
    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}
