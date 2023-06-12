#pragma once

#include "Datastructures.hpp"
#include "Domain.hpp"
#include <memory>
#include <string>
#include <vector>


class Communication{
       public: 
         	Communication() = default;
         	Communication(int iproc, int jproc, Domain &domain);
		void communicate(Matrix<double> &A);                  //Function to communicate between the field across all parallel neighbouring threads
		double reduce_min(const Matrix<double> &A);
		double reduce_sum(const Matrix<double> &A);
	private:
		
		int x_dim{0};
		int y_dim{0};
		
		int _iproc{0};   //Processes in the x direction
               	int _jproc{0};   //Processes in the y direction
               	int _imax{0};    //Max value of x index in the global domain
               	int _jmax{0};    //Max value of y index in the global domain
               	
		int _il{0}; //i index of the sub-domain left side cells
        	int _ir{0}; //i index of the sub-domain right side cells
        	int _jb{0}; //j index of the sub-domain bottom side cells
        	int _jt{0}; //j index of the sub-domain top side cells
};
