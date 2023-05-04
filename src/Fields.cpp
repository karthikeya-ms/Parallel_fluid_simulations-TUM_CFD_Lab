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

//Function implemented for second task.
void Fields::calculate_fluxes(Grid &grid, double gamma) {
	
	int i_idx{0};
	int j_idx{0};
	for (auto const cell : grid.fluid_cells()){
		
		i_idx = cell.i();
		j_idx = cell.j();
		
		f(i_idx, j_idx) = u(i_idx, j_idx) + _dt*(_nu*(d2udx2(i_idx,j_idx) + d2udy2(i_idx, j_idx)) - du2dx(i_idx, j_idx, gamma) - duvdy(i_idx, j_idx, gamma) + gx);
		
		g(i_idx, j_idx) = v(i_idx, j_idx) + _dt*(_nu*(d2vdx2(i_idx,j_idx) + d2vdy2(i_idx, j_idx)) - dv2dy(i_idx, j_idx, gamma) - duvdx(i_idx, j_idx, gamma) + gy);	
	}

}

void Fields::calculate_rs(Grid &grid) {}

void Fields::calculate_velocities(Grid &grid) {}

double Fields::calculate_dt(Grid &grid) { return _dt; }

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }



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
