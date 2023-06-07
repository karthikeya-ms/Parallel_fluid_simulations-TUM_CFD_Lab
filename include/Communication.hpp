#pragma once

#include <memory>
#include <string>
#include <vector>


class Communication{
       Public: 
         	Communication() = default;
         	Communication(int iproc, int jproc, int imax, int jmax, int argn, char **args);
         	static void init_parallel(int argn, char **args);        //Initialize the communication
		static void finalize();                                  //Finalize the communication
		static void communicate(const Matrix<double> &A);                  //Function to coomunicate between the field across all parallel boundaries
		static double reduce_min(const Matrix<double> &A);
		static double reduce_sum(const Matrix<double> &A);
	Private:
		int _argn;
		char **_args
		
		int _iproc{0};   //Processes in the x direction
               	int _jproc{0};   //Processes in the y direction
               	int _imax{0};    //Max value of x index in the global domain
               	int _jmax{0};    //Max value of y index in the global domain
               	int *_myrank{0}; //Rank of the current process
               	int *_il{0};     //Left cells i index of the sub-domain
               	int *_ir{0};     //Right cells i index of the sub-domain
               	int *_jb{0};     //Bottom cells j index of the sub-domain
               	int *_jt{0};     //Top cells j index of the sub-domain
               	int *_l_rank{0}; //Rank of the left neighbour cells
               	int *_r_rank{0}; //Rank of the right neighbour cells
               	int *_b_rank{0}; //Rank of the bottom neighbour cells
               	int *_t_rank{0}; //Rank of the top neighbour cells
               	int *_om_i{0};   //i index of the sub-domain
               	int *_om_j{0};   //j index of the sub domain
      };
