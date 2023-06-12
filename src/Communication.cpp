#include "Communication.hpp"
#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <numeric>

Communication::Communication(int iproc, int jproc, Domain &domain) : _iproc(iproc), _jproc(jproc){

	x_dim = domain.domain_size_x;
	y_dim = domain.domain_size_y;
	_il = domain.imin; //i index of the sub-domain left side cells
        _ir = domain.imax; //i index of the sub-domain right side cells
        _jb = domain.jmin; //j index of the sub-domain bottom side cells
        _jt = domain.jmax; //j index of the sub-domain top side cells
	
	}
	
void Communication::communicate(Matrix<double> &A) {

	int size = _iproc*_jproc;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int _myrank; //Rank of the current process
	MPI_Comm_rank(MPI_COMM_WORLD, &_myrank);
	
	//Calculating the omega (i index (0m_i)) and omega j (j index (om_j)) of the sub-domain
	int _om_i = ((_myrank) % _iproc) + 1;
        int _om_j = (((_myrank + 1) - _om_i)/_iproc) + 1;
        
  	int l_rank{0}; //Rank of the left neighbour sub domain
  	int r_rank{0}; //Rank of the right neighbour sub domain
  	int b_rank{0}; //Rank of the bottom neighbour sub domain
  	int t_rank{0}; //Rank of the top neighbour sub domain

  	if(_om_i == 1)       {l_rank = MPI_PROC_NULL;}
  	else                 {l_rank = _myrank - 1;}
	if(_om_i == _iproc)  {r_rank = MPI_PROC_NULL;}
  	else                 {r_rank = _myrank + 1;}

	if(_om_j == 1)       {b_rank = MPI_PROC_NULL;}
  	else                 {b_rank = _myrank - _iproc;}
	if(_om_j == _jproc)  {t_rank = MPI_PROC_NULL;}
  	else                 {t_rank = _myrank + _iproc;}
  	
  	
	double SDatax[x_dim] = {0}; //Send Buffer in the x direction for top and bottom cells
	double SDatay[y_dim] = {0}; //Send Buffer in the y direction for left and right cells
	double RDatax[x_dim] = {0}; //Recieve Buffer in the x direction for top and bottom cells
	double RDatay[y_dim] = {0}; //Recieve Buffer in the y direction for left and right cells
	
//############################################# CASES WHERE iproc and jproc ARE NOT EQUAL TO 1 STARTS HERE ###################################################################
	
//########################################################## CORNER SUB-DOMAINS START HERE ###################################################################################

	//LEFT BOTTOM corner domain -> Send and Recieve from the right, Send and recieve from the top
	if (l_rank == MPI_PROC_NULL && b_rank == MPI_PROC_NULL && t_rank != MPI_PROC_NULL && r_rank != MPI_PROC_NULL){
		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(x_dim-1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, r_rank, 0, &RDatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(x_dim,j) = RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, y_dim-1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, t_rank, 1, &RDatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, y_dim) = RDatax[i];
  		}
  	}
  		
  	//RIGHT BOTTOM corner domain -> Send and Recieve from the left, Send and recieve from the top
	else if (l_rank != MPI_PROC_NULL && b_rank == MPI_PROC_NULL && t_rank != MPI_PROC_NULL && r_rank == MPI_PROC_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, l_rank, 0, &RDatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A(0, j) = RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, y_dim-1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, t_rank, 1, &RDatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, y_dim) = RDatax[i];
  		}
  	}
  	
  	//LEFT TOP corner domain -> Send and Recieve from the right, Send and recieve from the bottom
	else if (l_rank == MPI_PROC_NULL && b_rank != MPI_PROC_NULL && t_rank == MPI_PROC_NULL && r_rank != MPI_PROC_NULL){
		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(x_dim-1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, r_rank, 0, &RDatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(x_dim, j) = RDatay[j];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, 1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, b_rank, 1, &RDatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, 0) = RDatax[i];
  		}
  	}
  	
  	//RIGHT TOP corner domain -> Send and Recieve from the left, Send and recieve from the bottom
	else if (l_rank != MPI_PROC_NULL && b_rank != MPI_PROC_NULL && t_rank == MPI_PROC_NULL && r_rank == MPI_PROC_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, l_rank, 0, &RDatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A(0, j) = RDatay[j];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, 1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, b_rank, 1, &RDatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, 0) = RDatax[i];
  		}
  	}
  	
//########################################################## CORNER SUB-DOMAINS END HERE ################################################################################### 		
  		
//########################################################## SIDE SUB-DOMAINS START HERE ###################################################################################  

	//LEFT SIDE domain -> Send and Recieve from the right, Send and recieve from the top, Send and recieve from the bottom
	else if (l_rank == MPI_PROC_NULL && b_rank != MPI_PROC_NULL && t_rank != MPI_PROC_NULL && r_rank != MPI_PROC_NULL){
		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(x_dim-1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, r_rank, 0, &RDatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(x_dim, j) = RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, y_dim-1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, t_rank, 1, &RDatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, y_dim) = RDatax[i];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, 1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, b_rank, 1, &RDatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, 0) = RDatax[i];
  		}
  	}
  	
 	//TOP SIDE domain -> Send and Recieve from the right, Send and recieve from the left, Send and recieve from the bottom
	else if (l_rank != MPI_PROC_NULL && b_rank != MPI_PROC_NULL && t_rank == MPI_PROC_NULL && r_rank != MPI_PROC_NULL){
		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(x_dim-1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, r_rank, 0, &RDatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(x_dim, j) = RDatay[j];
  		}
  		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, l_rank, 0, &RDatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A(0, j) = RDatay[j];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, 1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, b_rank, 1, &RDatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, 0) = RDatax[i];
  		}
  	}
  	
  	//RIGHT SIDE domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the bottom
	else if (l_rank != MPI_PROC_NULL && b_rank != MPI_PROC_NULL && t_rank != MPI_PROC_NULL && r_rank == MPI_PROC_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, l_rank, 0, &RDatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A(0, j) = RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, y_dim-1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, t_rank, 1, &RDatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, y_dim) = RDatax[i];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, 1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, b_rank, 1, &RDatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, 0) = RDatax[i];
  		}
  	}
  	
  	//BOTTOM SIDE domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the right
	else if (l_rank != MPI_PROC_NULL && b_rank == MPI_PROC_NULL && t_rank != MPI_PROC_NULL && r_rank != MPI_PROC_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, l_rank, 0, &RDatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A(0, j) = RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, y_dim-1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, t_rank, 1, &RDatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, y_dim) = RDatax[i];
  		}
  		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(x_dim-1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, r_rank, 0, &RDatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(x_dim, j) = RDatay[j];
  		}
  	}
  	
//########################################################## SIDE SUB-DOMAINS END HERE ###################################################################################  

//########################################################## INNER SUB-DOMAINS START HERE ###################################################################################	

//INNER domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the right, Send and recieve from the bottom
	else if (l_rank != MPI_PROC_NULL && b_rank != MPI_PROC_NULL && t_rank != MPI_PROC_NULL && r_rank != MPI_PROC_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, l_rank, 0, &RDatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A(0, j) = RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, y_dim-1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, t_rank, 1, &RDatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, y_dim) = RDatax[i];
  		}
  		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(x_dim-1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, r_rank, 0, &RDatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(x_dim, j) = RDatay[j];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, 1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, b_rank, 1, &RDatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, 0) = RDatax[i];
  		}
  	}
  	
//########################################################## INNER SUB-DOMAINS END HERE ########################################################################################

//############################################# CASES WHERE iproc and jproc ARE NOT EQUAL TO 1 END HERE #######################################################################	

//############################################# CASES WHERE iproc and jproc ARE EQUAL TO 1 START HERE #########################################################################
//############################################# CASES WHERE iproc IS EQUAL TO 1 START HERE ####################################################################################
	else if (l_rank == MPI_PROC_NULL && b_rank == MPI_PROC_NULL && t_rank != MPI_PROC_NULL && r_rank == MPI_PROC_NULL){
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, y_dim-1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, t_rank, 1, &RDatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, y_dim) = RDatax[i];
  		}
  	}
  	
  	else if (l_rank == MPI_PROC_NULL && b_rank != MPI_PROC_NULL && t_rank != MPI_PROC_NULL && r_rank == MPI_PROC_NULL){
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, y_dim-1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, t_rank, 1, &RDatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, y_dim) = RDatax[i];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, 1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, b_rank, 1, &RDatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, 0) = RDatax[i];
  		}
  	}
  	
  	else if (l_rank == MPI_PROC_NULL && b_rank != MPI_PROC_NULL && t_rank == MPI_PROC_NULL && r_rank == MPI_PROC_NULL){
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			SDatax[i] = A(i, 1);
  		}
  		MPI_Sendrecv(&SDatax, x_dim, MPI_DOUBLE, b_rank, 1, &RDatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A(i, 0) = RDatax[i];
  		}
  	}
//############################################# CASES WHERE iproc IS EQUAL TO 1 END HERE #################################################################################### 
//############################################# CASES WHERE jproc IS EQUAL TO 1 START HERE ##################################################################################	

	else if (l_rank == MPI_PROC_NULL && b_rank == MPI_PROC_NULL && t_rank == MPI_PROC_NULL && r_rank != MPI_PROC_NULL){
  		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(x_dim-1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, r_rank, 0, &RDatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(x_dim, j) = RDatay[j];
  		}
  	}
  	
  	else if (l_rank != MPI_PROC_NULL && b_rank == MPI_PROC_NULL && t_rank == MPI_PROC_NULL && r_rank != MPI_PROC_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, l_rank, 0, &RDatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(0, j) = RDatay[j];
  		}
  		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(x_dim-1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, r_rank, 0, &RDatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A(x_dim, j) = RDatay[j];
  		}
  	}
  	
  	else if (l_rank != MPI_PROC_NULL && b_rank == MPI_PROC_NULL && t_rank == MPI_PROC_NULL && r_rank == MPI_PROC_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			SDatay[j] = A(1, j);
  		}
  		MPI_Sendrecv(&SDatay, y_dim, MPI_DOUBLE, l_rank, 0, &RDatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A(0, j) = RDatay[j];
  		}
  	}
//############################################# CASES WHERE jproc IS EQUAL TO 1 END HERE ######################################################################################
//############################################# CASES WHERE iproc and jproc ARE EQUAL TO 1 END HERE ###########################################################################						
}  						
  	
double Communication::reduce_min(const Matrix<double> &A){
	
	std::vector<double> vec;
	for(int i = 1; i < x_dim; ++i){
		for(int j = 1; j < y_dim; ++j){
			vec.push_back(A(i,j));
		}
	}
	double localMin; 
	localMin = *std::min_element(vec.begin(),vec.end());
	double globalMin;
	MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return globalMin;
}

double Communication::reduce_sum(const Matrix<double> &A){
	
	double localSum = 0;
	std::vector<double> vec{0};
	for(int i = 1; i < x_dim; ++i){
		for(int j = 1; j<y_dim; ++j){
			vec.push_back(A(i,j));
		}
	}
	localSum = std::accumulate(vec.begin(), vec.end(), 0);
	double globalSum;
	MPI_Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalSum;
}


  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
