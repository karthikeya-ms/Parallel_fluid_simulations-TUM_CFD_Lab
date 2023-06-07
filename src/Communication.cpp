
#include <mpi.h>
#include <algorithm>
#include <iostream>

Communication::Communication(int iproc, int jproc, int imax, int jmax, int argn, char **args) : _iproc(iproc), _jproc(jproc), _imax(imax), _jmax(jmax), _argn(argn), _**args(**args) {
	
	*_myrank = init_parallel(_argn, _**args);
	*_om_i = (*_myrank % _iproc) + 1;
        *_om_j = ((*_myrank + 1 - *_om_i)/_iproc) + 1;
        
        *_il = (*om_i - 1) * (_imax/_iproc) + 1;
        if(*_om_i != _iproc) {_*ir = (*_om_i) * (_imax/_iproc);}
        else {_*ir = _imax;}
	
	*_jb = (*_om_j - 1) * (_jmax/_jproc) + 1;
  	if(*_om_j != _jproc) {_*jt = (*_om_j) * (_jmax/_jproc);}
  	else {*_jt = _jmax;}

  	if(*_il == 1)     {*l_rank = MPI_PROC_NULL;}
  	else              {*l_rank = *_myrank - 1;}
	if(*_ir == _imax) {*r_rank = MPI_PROC_NULL;}
  	else              {*r_rank = *_myrank + 1;}

	if(*_jb == 1)     {*b_rank = MPI_PROC_NULL;}
  	else              {*b_rank = *_myrank - _iproc;}
	if(*_jt == _jmax) {*t_rank = MPI_PROC_NULL;}
  	else              {*t_rank = *_myrank + _iproc;}
	}
	
static void Communication::Communicate(const Matrix<double> &A) {

	int x_dim = _ir - _il + 1;
	int y_dim = _jt - _jb + 1;
	double SDatax[x_dim] = {0};
	double SDatay[y_dim] = {0};
	double RDatax[x_dim] = {0};
	double RDatax[y_dim] = {0};
	
//########################################################## CORNER SUB-DOMAINS START HERE ###################################################################################

	//LEFT BOTTOM corner domain -> Send and Recieve from the right, Send and recieve from the top
	if (*l_rank == NPI_PROG_NULL && *b_rank == NPI_PROG_NULL && *t_rank != NPI_PROG_NULL && *r_rank != NPI_PROG_NULL){
		//Send and Recieve from the right
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[xdim-1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[x_dim][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][y_dim-1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][y_dm] = RDatax[i-1];
  		}
  	}
  		
  	//RIGHT BOTTOM corner domain -> Send and Recieve from the left, Send and recieve from the top
	else if (*l_rank != NPI_PROG_NULL && *b_rank == NPI_PROG_NULL && *t_rank != NPI_PROG_NULL && *r_rank == NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[0][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][y_dim-1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][y_dm] = RDatax[i-1];
  		}
  	}
  	
  	//LEFT TOP corner domain -> Send and Recieve from the right, Send and recieve from the bottom
	else if (*l_rank == NPI_PROG_NULL && *b_rank != NPI_PROG_NULL && *t_rank == NPI_PROG_NULL && *r_rank != NPI_PROG_NULL){
		//Send and Recieve from the right
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[x_dim-1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[xdim][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][0] = RDatax[i-1];
  		}
  	}
  	
  	//RIGHT TOP corner domain -> Send and Recieve from the left, Send and recieve from the bottom
	else if (*l_rank != NPI_PROG_NULL && *b_rank != NPI_PROG_NULL && *t_rank == NPI_PROG_NULL && *r_rank != NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[0][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][0] = RDatax[i-1];
  		}
  	}
  	
//########################################################## CORNER SUB-DOMAINS END HERE ################################################################################### 		
  		
//########################################################## SIDE SUB-DOMAINS START HERE ###################################################################################  

	//LEFT SIDE domain -> Send and Recieve from the right, Send and recieve from the top, Send and recieve from the bottom
	else if (*l_rank == NPI_PROG_NULL && *b_rank != NPI_PROG_NULL && *t_rank != NPI_PROG_NULL && *r_rank != NPI_PROG_NULL){
		//Send and Recieve from the right
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[xdim-1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[x_dim][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][y_dim-1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][y_dm] = RDatax[i-1];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][0] = RDatax[i-1];
  		}
  	}
  	
 	//TOP SIDE domain -> Send and Recieve from the right, Send and recieve from the left, Send and recieve from the bottom
	else if (*l_rank != NPI_PROG_NULL && *b_rank != NPI_PROG_NULL && *t_rank == NPI_PROG_NULL && *r_rank != NPI_PROG_NULL){
		//Send and Recieve from the right
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[xdim-1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[x_dim][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the left
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[0][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][0] = RDatax[i-1];
  		}
  	}
  	
  	//RIGHT SIDE domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the bottom
	else if (*l_rank != NPI_PROG_NULL && *b_rank != NPI_PROG_NULL && *t_rank != NPI_PROG_NULL && *r_rank == NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[0][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][y_dim-1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][y_dm] = RDatax[i-1];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][0] = RDatax[i-1];
  		}
  	}
  	
  	//BOTTOM SIDE domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the right
	else if (*l_rank != NPI_PROG_NULL && *b_rank == NPI_PROG_NULL && *t_rank != NPI_PROG_NULL && *r_rank != NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[0][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][y_dim-1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][y_dm] = RDatax[i-1];
  		}
  		//Send and Recieve from the right
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[xdim-1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[x_dim][j] = RDatay[j-1];
  		}
  	}
  	
//########################################################## SIDE SUB-DOMAINS END HERE ###################################################################################  

//########################################################## INNER SUB-DOMAINS START HERE ###################################################################################	

//INNER domain -> Send and Recieve from the left, Send and recieve from the top, Send and recieve from the right, Send and recieve from the bottom
	else if (*l_rank != NPI_PROG_NULL && *b_rank != NPI_PROG_NULL && *t_rank != NPI_PROG_NULL && *r_rank != NPI_PROG_NULL){
		//Send and Recieve from the left
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[0][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the top
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][y_dim-1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][y_dm] = RDatax[i-1];
  		}
  		//Send and Recieve from the right
		for (int j = 1; j <=y_dim; ++j){
    			SDatay[j-1] = A[xdim-1][j];
  		}
  		MPI_Send(SDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatay, ydim, MPI_DOUBLE, r_rank, 1, MPI_COMM_WORLD, status);
  		for (int j = 1; j <=y_dim; ++j){
    			A[x_dim][j] = RDatay[j-1];
  		}
  		//Send and Recieve from the bottom
		for (int i = 1; i <=x_dim; ++i){
    			SDatax[i-1] = A[i][1];
  		}
  		MPI_Send(SDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(RDatax, xdim, MPI_DOUBLE, t_rank, 1, MPI_COMM_WORLD, status);
  		for (int i = 1; i <=x_dim; ++i){
    			A[i][0] = RDatax[i-1];
  		}
  	}
  	
//########################################################## INNER SUB-DOMAINS END HERE ########################################################################################		
  				
}  		
  		
static int Communication::init_parallel(int argn, char **args){

	MPI_Init(&argn, &args);
	int size = _iproc*_jproc;
	MPI _Comm_size(MPI_COMM_WORLD, &size);
	int r;
	MPI_Comm_rank(MPI_COMM_WORLD, &r);
	return r;
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
	return globalMin;
}



  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
  		
