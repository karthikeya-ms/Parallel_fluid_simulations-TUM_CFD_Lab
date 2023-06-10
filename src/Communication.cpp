
#include <mpi.h>
#include <algorithm>
#include <iostream>

Communication::Communication(int iproc, int jproc, int imax, int jmax, int argn, char **args) : _iproc(iproc), _jproc(jproc), _imax(imax), _jmax(jmax), _argn(argn), _**args(**args) {
	
	}
	
static void Communication::Communicate(const Matrix<double> &A) {

	MPI_Init(&argc, &argv);
	int size = _iproc*_jproc;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int _myrank; //Rank of the current process
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	//Calculating the omega (i index (0m_i)) and omega j (j index (om_j)) of the sub-domain
	int _om_i = (_myrank % _iproc) + 1;
        int _om_j = ((_myrank + 1 - _om_i)/_iproc) + 1;
        
        int _il{0}; //i index of the sub-domain left side cells
        int _ir{0}; //i index of the sub-domain right side cells
        int _jb{0}; //j index of the sub-domain bottom side cells
        int _jt{0}; //j index of the sub-domain top side cells
        
        il = (om_i - 1) * (_imax/_iproc) + 1;
        if(_om_i != _iproc) {_ir = (_om_i) * (_imax/_iproc);}
        else {_ir = _imax;}
	
	_jb = (_om_j - 1) * (_jmax/_jproc) + 1;
  	if(_om_j != _jproc) {_jt = (_om_j) * (_jmax/_jproc);}
  	else {_jt = _jmax;}
  	
  	int l_rank{0}; //Rank of the left neighbour sub domain
  	int r_rank{0}; //Rank of the right neighbour sub domain
  	int b_rank{0}; //Rank of the bottom neighbour sub domain
  	int t_rank{0}; //Rank of the top neighbour sub domain

  	if(_il == 1)      {l_rank = MPI_PROC_NULL;}
  	else              {l_rank = _myrank - 1;}
	if(_ir == _imax)  {r_rank = MPI_PROC_NULL;}
  	else              {r_rank = _myrank + 1;}

	if(_jb == 1)      {b_rank = MPI_PROC_NULL;}
  	else              {b_rank = _myrank - _iproc;}
	if(_jt == _jmax)  {t_rank = MPI_PROC_NULL;}
  	else              {t_rank = _myrank + _iproc;}
  	
  	int x_dim = _ir - _il + 1;
	int y_dim = _jt - _jb + 1;
	double *SDatax[x_dim] = {0}; //Send Buffer in the x direction for top and bottom cells
	double *SDatay[y_dim] = {0}; //Send Buffer in the y direction for left and right cells
	double *RDatax[x_dim] = {0}; //Recieve Buffer in the x direction for top and bottom cells
	double *RDatay[y_dim] = {0}; //Recieve Buffer in the y direction for left and right cells
	
//########################################################## CORNER SUB-DOMAINS START HERE ###################################################################################

	//LEFT BOTTOM corner domain -> Send and Recieve from the right, Send and recieve from the top
	if (l_rank == NPI_PROG_NULL && b_rank == NPI_PROG_NULL && t_rank != NPI_PROG_NULL && r_rank != NPI_PROG_NULL){
		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[x_dim-1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, r_rank, 0, Rdatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A[x_dim][j] = *RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][y_dim-1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, t_rank, 1, Rdatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][y_dm] = *RDatax[i];
  		}
  	}
  		
  	//RIGHT BOTTOM corner domain -> Send and Recieve from the left, Send and recieve from the top
	else if (l_rank != NPI_PROG_NULL && b_rank == NPI_PROG_NULL && t_rank != NPI_PROG_NULL && r_rank == NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, l_rank, 0, Rdatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A[0][j] = *RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][y_dim-1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, t_rank, 1, Rdatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][y_dm] = *RDatax[i];
  		}
  	}
  	
  	//LEFT TOP corner domain -> Send and Recieve from the right, Send and recieve from the bottom
	else if (l_rank == NPI_PROG_NULL && b_rank != NPI_PROG_NULL && t_rank == NPI_PROG_NULL && r_rank != NPI_PROG_NULL){
		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[x_dim-1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, r_rank, 0, Rdatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A[x_dim][j] = *RDatay[j];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, b_rank, 1, Rdatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][0] = *RDatax[i];
  		}
  	}
  	
  	//RIGHT TOP corner domain -> Send and Recieve from the left, Send and recieve from the bottom
	else if (l_rank != NPI_PROG_NULL && b_rank != NPI_PROG_NULL && t_rank == NPI_PROG_NULL && r_rank == NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, l_rank, 0, Rdatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A[0][j] = *RDatay[j];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, b_rank, 1, Rdatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][0] = *RDatax[i];
  		}
  	}
  	
//########################################################## CORNER SUB-DOMAINS END HERE ################################################################################### 		
  		
//########################################################## SIDE SUB-DOMAINS START HERE ###################################################################################  

	//LEFT SIDE domain -> Send and Recieve from the right, Send and recieve from the top, Send and recieve from the bottom
	else if (l_rank == NPI_PROG_NULL && b_rank != NPI_PROG_NULL && t_rank != NPI_PROG_NULL && r_rank != NPI_PROG_NULL){
		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[x_dim-1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, r_rank, 0, Rdatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A[x_dim][j] = *RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][y_dim-1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, t_rank, 1, Rdatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][y_dm] = *RDatax[i];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, b_rank, 1, Rdatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][0] = *RDatax[i];
  		}
  	}
  	
 	//TOP SIDE domain -> Send and Recieve from the right, Send and recieve from the left, Send and recieve from the bottom
	else if (l_rank != NPI_PROG_NULL && b_rank != NPI_PROG_NULL && t_rank == NPI_PROG_NULL && r_rank != NPI_PROG_NULL){
		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[x_dim-1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, r_rank, 0, Rdatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A[x_dim][j] = *RDatay[j];
  		}
  		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, l_rank, 0, Rdatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A[0][j] = *RDatay[j];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, b_rank, 1, Rdatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][0] = *RDatax[i];
  		}
  	}
  	
  	//RIGHT SIDE domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the bottom
	else if (l_rank != NPI_PROG_NULL && b_rank != NPI_PROG_NULL && t_rank != NPI_PROG_NULL && r_rank == NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, l_rank, 0, Rdatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A[0][j] = *RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][y_dim-1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, t_rank, 1, Rdatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][y_dm] = *RDatax[i];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, b_rank, 1, Rdatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][0] = *RDatax[i];
  		}
  	}
  	
  	//BOTTOM SIDE domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the right
	else if (l_rank != NPI_PROG_NULL && b_rank == NPI_PROG_NULL && t_rank != NPI_PROG_NULL && r_rank != NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, l_rank, 0, Rdatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A[0][j] = *RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][y_dim-1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, t_rank, 1, Rdatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][y_dm] = *RDatax[i];
  		}
  		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[x_dim-1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, r_rank, 0, Rdatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A[x_dim][j] = *RDatay[j];
  		}
  	}
  	
//########################################################## SIDE SUB-DOMAINS END HERE ###################################################################################  

//########################################################## INNER SUB-DOMAINS START HERE ###################################################################################	

//INNER domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the right, Send and recieve from the bottom
	else if (l_rank != NPI_PROG_NULL && b_rank != NPI_PROG_NULL && t_rank != NPI_PROG_NULL && r_rank != NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, l_rank, 0, Rdatay, y_dim, MPI_DOUBLE, l_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		
  		for (int j = 1; j < y_dim; ++j){
    			A[0][j] = *RDatay[j];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][y_dim-1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, t_rank, 1, Rdatax, x_dim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][y_dm] = *RDatax[i];
  		}
  		//Send and Recieve from the right
		for (int j = 1; j < y_dim; ++j){
    			*SDatay[j] = A[x_dim-1][j];
  		}
  		MPI_Sendrecv(Sdatay, y_dim, MPI_DOUBLE, r_rank, 0, Rdatay, y_dim, MPI_DOUBLE, r_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int j = 1; j < y_dim; ++j){
    			A[x_dim][j] = *RDatay[j];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i < x_dim; ++i){
    			*SDatax[i] = A[i][1];
  		}
  		MPI_Sendrecv(Sdatax, x_dim, MPI_DOUBLE, b_rank, 1, Rdatax, x_dim, MPI_DOUBLE, b_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  		for (int i = 1; i < x_dim; ++i){
    			A[i][0] = *RDatax[i];
  		}
  	}
  	
//########################################################## INNER SUB-DOMAINS END HERE ########################################################################################		
  MPI_Finalize();				
}  						
  	
static double Communicate::reduce_min(const Matrix<double> &A){
	
	for(int i = 1; i<x_dim; ++i){
		for(int j = 1; j<y_dim; ++j){
			vec.push_back(A[i][j]);
		}
	}
	double localMin = *std::min_element(vec.begin(),vec.end());
	double globalMin;
	MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return globalMin;
}

static double Communicate::reduce_sum(const Matrix<double> &A){
	
	double localSum;
	for(int i = 1; i<x_dim; ++i){
		for(int j = 1; j<y_dim; ++j){
			vec.push_back(A[i][j]);
		}
	}
	*std::accumulate(vec.begin(),vec.end(), localSum);
	double globalSum;
	MPI_Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalSum;
}


  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
