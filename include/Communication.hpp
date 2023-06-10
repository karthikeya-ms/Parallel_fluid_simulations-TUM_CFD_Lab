#pragma once

#include <memory>
#include <string>
#include <vector>


class Communication{
       Public: 
         	Communication() = default;
         	Communication(int iproc, int jproc, int imax, int jmax, int argn, char **args);
		static void communicate(const Matrix<double> &A);                  //Function to communicate between the field across all parallel neighbouring threads
		static double reduce_min(const Matrix<double> &A);
		static double reduce_sum(const Matrix<double> &A);
	Private:
		int _argn;
		char **_args
		
		int _iproc{0};   //Processes in the x direction
               	int _jproc{0};   //Processes in the y direction
               	int _imax{0};    //Max value of x index in the global domain
               	int _jmax{0};    //Max value of y index in the global domain
      };
