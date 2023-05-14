#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau), _imax(imax), _jmax(jmax) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

//Added: Function implemented for second task.
void Fields::calculate_fluxes(Grid &grid, Discretization &discretization) {
	int i_idx{0};
	int j_idx{0};

	for (int i_idx = 1; i_idx < grid.imax() + 1; i_idx++){
		for (int j_idx = 1; j_idx < grid.jmax() + 1; j_idx++){
			_F(i_idx, j_idx) = _U(i_idx, j_idx) + _dt*(_nu*discretization.diffusion(_U, i_idx, j_idx) + discretization.convection_U(_U, _V, i_idx, j_idx) + _gx);	
			_G(i_idx, j_idx) = _V(i_idx, j_idx) + _dt*(_nu*discretization.diffusion(_V, i_idx, j_idx) + discretization.convection_V(_U, _V, i_idx, j_idx) + _gy);
		}
	}
}

//Added: Function implemented for third task.
void Fields::calculate_rs(Grid &grid) {
	int i_idx{0};
	int j_idx{0};

	for (int i_idx = 1; i_idx < grid.imax() + 1; i_idx++){
		for (int j_idx = 1; j_idx < grid.jmax() + 1; j_idx++){
			_RS(i_idx, j_idx) = ((f(i_idx, j_idx) - f(i_idx - 1, j_idx))/grid.dx() + (g(i_idx, j_idx) - g(i_idx, j_idx - 1))/grid.dy())/_dt;
		}
	}
}

//Added: Function implemented for sixth task.
void Fields::calculate_velocities(Grid &grid) {

	for (int i_idx = 1; i_idx < grid.imax() + 1; i_idx++){
		for (int j_idx = 1; j_idx < grid.jmax() + 1; j_idx++){
			_U(i_idx, j_idx) = f(i_idx, j_idx) - _dt*dPdx(i_idx, j_idx, grid);
			_V(i_idx, j_idx) = g(i_idx, j_idx) - _dt*dPdy(i_idx, j_idx, grid);
		}
	}
}

//Added: Function implemented for seventh task.
double Fields::calculate_dt(Grid &grid) { 

	//Search for the absolute maximums of both U and V matrices.
	std::vector<double> vector_umax;
	std::vector<double> vector_vmax;
	for (auto i_idx = 0; i_idx < _imax + 1; i_idx++){
		std::vector<double> row_u = _U.get_row(i_idx);
		std::vector<double> row_v = _V.get_row(i_idx);
		
		for(auto j_idx = 0; j_idx < _jmax + 1; j_idx++){
		row_u[j_idx] = abs(row_u[j_idx]);
		row_v[j_idx] = abs(row_v[j_idx]);
		}
		
		vector_umax.push_back(*std::max_element(std::begin(row_u), std::end(row_u)));
		vector_vmax.push_back(*std::max_element(std::begin(row_v), std::end(row_v)));
	}
	double umax = *std::max_element(begin(vector_umax), end(vector_umax));
	double vmax = *std::max_element(begin(vector_vmax), end(vector_vmax));
         double _dta = std::min(grid.dx()/umax, grid.dy()/vmax);
         double _dtb = pow((1/pow(grid.dx(), 2.0) + 1/pow(grid.dy(), 2.0)), -1.0)/(2*_nu);
	_dt = _tau*std::min(_dta, _dtb);
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

double Fields::dPdx(int i_idx, int j_idx, Grid &grid){
	return (p(i_idx + 1, j_idx) - p(i_idx, j_idx))/grid.dx();
}

double Fields::dPdy(int i_idx, int j_idx, Grid &grid){
	return (p(i_idx, j_idx + 1) - p(i_idx, j_idx))/grid.dy();
}
