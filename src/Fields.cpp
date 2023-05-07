#include "Fields.hpp"

#include <algorithm>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

//Added: Function implemented for second task.
void Fields::calculate_fluxes(Grid &grid, double gamma) {
	
	int i_idx{0};
	int j_idx{0};
	
	//Have to wonder if doing a nested loop is actually faster than going over the fluid cells (just one loop), and asking if their indexes are in the desired limits.
	//Doing the second approach requires for 2*(imax - 1)*(jmax - 1) getters for indexes and same number of comparisons, maybe worse.
	//for (auto const cell : grid.fluid_cells()){
		//i_idx = cell->i();
		//j_idx = cell->j();
	for (int i_idx = 1; i_idx < grid.imax(); i_idx++){
		for (int j_idx = 1; j_idx < grid.jmax(); j_idx++){
			f(i_idx, j_idx) = u(i_idx, j_idx) + _dt*(_nu*(d2udx2(i_idx,j_idx) + d2udy2(i_idx, j_idx)) - du2dx(i_idx, j_idx, gamma) - duvdy(i_idx, j_idx, gamma) + _gx);	
			g(i_idx, j_idx) = v(i_idx, j_idx) + _dt*(_nu*(d2vdx2(i_idx,j_idx) + d2vdy2(i_idx, j_idx)) - dv2dy(i_idx, j_idx, gamma) - duvdx(i_idx, j_idx, gamma) + _gy);
		}
	}

}

//Added: Function implemented for third task.
void Fields::calculate_rs(Grid &grid) {

	int i_idx{0};
	int j_idx{0};
	
	//for (auto const cell : grid.fluid_cells()){
		//i_idx = cell->i();
		//j_idx = cell->j();
	for (int i_idx = 1; i_idx < grid.imax(); i_idx++){
		for (int j_idx = 1; j_idx < grid.jmax(); j_idx++){
			rs(i_idx, j_idx) = ((f(i_idx, j_idx) - f(i_idx - 1, j_idx))/grid.dx() + (g(i_idx, j_idx) - g(i_idx, j_idx - 1))/grid.dy())/_dt;
		}
	}

}

//Added: Function implemented for sixth task.
void Fields::calculate_velocities(Grid &grid) {

	for (int i_idx = 1; i_idx < grid.imax(); i_idx++){
		for (int j_idx = 1; j_idx < grid.jmax(); j_idx++){
			u(i_idx, j_idx) = f(i_idx, j_idx) - _dt*dpdx(i_idx, j_idx);
			v(i_idx, j_idx) = g(i_idx, j_idx) - _dt*dpdy(i_idx, j_idx);
		}
	}
}

//Added: Function implemented for seventh task.
double Fields::calculate_dt(Grid &grid) { 

	//Search for the absolute maximums of both U and V matrices.
	std::vector<double> vector_umax;
	std::vector<double> vector_vmax;
	for (auto i_idx = 0; i_idx < _U.size(); i_idx++){
		std::valarray<double> row_u(_U[i_idx].data(), _U[i_idx].size());
		std::valarray<double> row_v(_V[i_idx].data(), _V[i_idx].size());
		row_u = abs(row_u);
		row_v = abs(row_v);
		
		vector_umax.append(*std::max_element(begin(row_u), end(row_u)));
		vector_vmax.append(*std::max_element(begin(row_v), end(row_v)));
	}
	double umax = *std::max_element(begin(vector_umax), end(vector_umax));
	double vmax = *std::max_element(begin(vector_vmax), end(vector_vmax));

	_dt = _tau*std::min(pow((1/pow(grid.dx(), 2.0) + 1/pow(grid.dy(), 2.0)), -1.0)/(2*_nu),  grid.dx()/umax, grid.dy/vmax);
	return _dt; 
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }


//Added: Implementation of helper functions containing derivative terms to reduce cluttering in flux calculations for second task.
double Fields::d2udx2(int i_idx, int j_idx){

	return (u(i_idx + 1, j_idx) - 2*u(i_idx, j_idx) + u(idx - 1, j_idx))/(pow(grid.dx(), 2.0));
}

double Fields::d2udy2(int i_idx, int j_idx){

	return (u(i_idx, j_idx + 1) - 2*u(i_idx, j_idx) + u(idx, j_idx - 1))/pow(grid.dy(), 2.0);
}

double Fields::du2dx(int i_idx, int j_idx, double gamma){

	return (1/grid.dx())*(pow((u(i_idx, j_idx) + u(idx + 1, j_idx))/2.0, 2.0) - pow((u(i_idx - 1, j_idx) + u(idx, j_idx))/2.0, 2.0)) + (gamma/grid.dx())*(abs(u(i_idx, j_idx) + u(i_idx + 1, j_idx))/2.0 * (u(i_idx, j_idx) - u(i_idx + 1, j_idx))/2.0 - abs(u(i_idx - 1, j_idx) + u(i_idx, j_idx))/2.0 *(u(i_idx - 1, j_idx) - u(i_idx, j_idx))/2.0);
}

double Fields::duvdy(int i_idx, int j_idx, double gamma){

 	return (1/grid.dy())*((v(i_idx, j_idx) + v(i_idx + 1, j_idx))/2.0 * (u(i_idx, j_idx) + u(i_idx, j_idx + 1))/2.0 - (v(i_idx, j_idx - 1) + v(i_idx + 1, j_idx - 1))/2.0 * (u(i_idx, j_idx - 1) + u(i_idx, j_idx))/2.0) +
(gamma/grid.dy())*(abs(v(i_idx,j_idx) + v(i_idx + 1, j_idx))/2.0 * (u(i_idx, j_idx) - u(i_idx, j_idx + 1))/2.0 - abs(v(i_idx, j_idx - 1) + v(i_idx + 1, j_idx - 1))/2 * (u(i_idx, j_idx - 1) - u(i_idx, j_idx))/2.0);
}

double Fields::duvdx(int i_idx, int j_idx, double gamma){

	return (1/grid.dx())*((u(i_idx,j_idx) + u(i_idx, j_idx + 1))/2.0 * (v(i_idx, j_idx) + v(i_idx + 1, j_idx))/2.0 - (u(i_idx - 1,j_idx) + u(i_idx - 1, j_idx + 1))/2.0 * (v(i_idx - 1, j_idx) + v(i_idx, j_idx))/2.0) +
(gamma/grid.dy())*(abs(v(i_idx,j_idx) + v(i_idx + 1, j_idx))/2.0 * (u(i_idx, j_idx)-u(i_idx, j_idx + 1))/2.0 - abs(v(i_idx,j_idx - 1) + v(i_idx + 1, j_idx - 1))/2.0 * (u(i_idx, j_idx - 1) - u(i_idx, j_idx))/2.0);
}

double Fields::dv2dy(int i_idx, int j_idx, double gamma){

	return (1/grid.dy())*(pow((v(i_idx, j_idx) + v(i_idx, j_idx + 1))/2.0, 2.0)) - pow((v(i_idx, j_idx - 1) + v(i_idx, j_idx))/2.0, 2.0)) + 
(gamma*1/grid.dy())*(abs(v(i_idx, j_idx) + v(i_idx , j_idx + 1))/2.0 * (v(i_idx, j_idx) - v(i_idx, j_idx + 1))/2.0 - abs(v(i_idx, j_idx - 1) + v(i_idx , j_idx))/2.0 * (v(i_idx, j_idx - 1) - v(i_idx, j_idx))/2.0);
}

double Fields::d2vdx2(int i_idx, int j_idx){

	return (v(i_idx + 1, j_idx) - 2*v(i_idx, j_idx) + v(i_idx - 1, j_idx))/pow(grid.dx(), 2.0);
}

double Fields::d2vdy2(int i_idx, int j_idx){

	return (v(i_idx, j_idx + 1) - 2*v(i_idx, j_idx) + v(i_idx, j_idx - 1))/pow(grid.dy(), 2.0);
}

double Fields::dpdx(int i_idx, int j_idx){

	return (p(i_idx + 1, j_idx) - p(i_idx, j_idx))/grid.dx();
}

double Fields::dpdy(int i_idx, int j_idx){

	return (p(i_idx, j_idx + 1) - p(i_idx, j_idx))/grid.dy();
}
