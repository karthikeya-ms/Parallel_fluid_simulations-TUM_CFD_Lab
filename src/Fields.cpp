#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <math.h>

Fields::Fields(double gx, double gy, double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double UIN, double VIN, double PI, double TI, double alpha, double beta)
    : _gx(gx), _gy(gy), _nu(nu), _dt(dt), _tau(tau), _imax(imax), _jmax(jmax), _alpha(alpha), _beta(beta) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
    
    _T = Matrix<double>(imax + 2, jmax + 2, TI);
}

//Added: Function implemented for second task.
void Fields::calculate_fluxes(Grid &grid, Discretization &discretization, bool energy_eq, Communication &communication) {
	int i_idx{0};
	int j_idx{0};
	
	double hydro_term_x{0};
	double hydro_term_y{0};

	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		
		//Add the necessary corrections coming from the Boussinesq approximation to the momentum equations (already using the new temperatures).
		if (energy_eq == true){
			hydro_term_x = -1*_beta*(_dt/2.0)*(_T(i_idx, j_idx) + _T(i_idx + 1, j_idx))*_gx;
			hydro_term_y = -1*_beta*(_dt/2.0)*(_T(i_idx, j_idx) + _T(i_idx, j_idx + 1))*_gy;
		}	
		else {
			hydro_term_x = _dt*_gx;
			hydro_term_y = _dt*_gy;
		}
		_F(i_idx, j_idx) = _U(i_idx, j_idx) + _dt*(_nu*discretization.diffusion(_U, i_idx, j_idx) + discretization.convection_U(_U, _V, i_idx, j_idx)) + hydro_term_x;
		_G(i_idx, j_idx) = _V(i_idx, j_idx) + _dt*(_nu*discretization.diffusion(_V, i_idx, j_idx) + discretization.convection_V(_U, _V, i_idx, j_idx)) + hydro_term_y;
	}
	communication.communicate(_F);
	communication.communicate(_G);
	
}

//Added: Function implemented for third task.
void Fields::calculate_rs(Grid &grid) {
	int i_idx{0};
	int j_idx{0};

	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		_RS(i_idx, j_idx) = (1/_dt) * ((_F(i_idx, j_idx) - _F(i_idx - 1, j_idx))/grid.dx() + (_G(i_idx, j_idx) - _G(i_idx, j_idx - 1))/grid.dy());	
	}
}

//Added: Function implemented for sixth task.
void Fields::calculate_velocities(Grid &grid, Communication &communication) {
	int i_idx{0};
	int j_idx{0};
	
	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		_U(i_idx, j_idx) = _F(i_idx, j_idx) - _dt*dPdx(i_idx, j_idx, grid);
		_V(i_idx, j_idx) = _G(i_idx, j_idx) - _dt*dPdy(i_idx, j_idx, grid);
	}
	communication.communicate(_U);
	communication.communicate(_V);
}

void Fields::calculate_temperatures(Grid &grid, Discretization &discretization, Communication &communication){
	Matrix<double> new_T = _T;
	int i_idx{0};
	int j_idx{0};
	
	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		new_T(i_idx, j_idx) = t(i_idx, j_idx) + _dt*(_alpha*discretization.diffusion(_T, i_idx, j_idx) - discretization.convection_T(_T, _U, _V, i_idx, j_idx));
	}
	_T = new_T;
	communication.communicate(_T);
}

double Fields::calculate_maxU(Grid &grid){
	
	int i_idx{0};
	int j_idx{0};
	double max_u{0};
	//double max_v{0};
	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		
		if (fabs(_U(i_idx, j_idx)) > max_u){
			max_u = fabs(_U(i_idx, j_idx));
		}
	}
	return max_u;
}
  
double Fields::calculate_maxV(Grid &grid){	

	int i_idx{0};
	int j_idx{0};
	double max_v{0};
	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		
		if (fabs(_V(i_idx, j_idx)) > max_v){
			max_v = fabs(_V(i_idx, j_idx));
		}
	}
	return max_v;
}

//Added: Function implemented for seventh task.
double Fields::calculate_dt(Grid &grid, bool energy_eq, double &max_u, double &max_v) { 

	double dta = std::min(grid.dx()/max_u, grid.dy()/max_v);
	double dtb = pow((1/pow(grid.dx(), 2.0) + 1/pow(grid.dy(), 2.0)), -1.0)/(2*_nu);
	double dt = _tau*std::min(dta, dtb);
	
	if (energy_eq == true){
		double dtc = pow((1/pow(grid.dx(), 2.0) + 1/pow(grid.dy(), 2.0)), -1.0)/(2*_alpha);
		dt = std::min(dt, dtc);
	}
	_dt = dt;
	return _dt; 
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }
double &Fields::t(int i, int j) { return _T(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }

double Fields::dt() const { return _dt; }

double Fields::dPdx(int i_idx, int j_idx, Grid &grid){
	return (1/grid.dx()) * (_P(i_idx + 1, j_idx) - _P(i_idx, j_idx));
}

double Fields::dPdy(int i_idx, int j_idx, Grid &grid){
	return (1/grid.dy()) * (_P(i_idx, j_idx + 1) - _P(i_idx, j_idx));
}
