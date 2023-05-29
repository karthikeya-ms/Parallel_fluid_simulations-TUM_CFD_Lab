#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double UIN, double VIN, double PI, double TI, double alpha, double beta, double gx, double gy)
    : _nu(nu), _dt(dt), _tau(tau), _imax(imax), _jmax(jmax), _alpha(alpha), _beta(beta), _gx(gx), _gy(gy) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
    
    _T = Matrix<double>(imax + 2, jmax + 2, TI);
}

//Added: Function implemented for second task.
void Fields::calculate_fluxes(Grid &grid, Discretization &discretization, bool energy_eq) {
	int i_idx{0};
	int j_idx{0};
	auto imax = grid.imax();
	auto jmax = grid.jmax();
	
	double hydro_term_x{0};
	double hydro_term_y{0};

	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		
		//Add the necessary corrections coming from the Boussinesq approximation to the momentum equations (already using the new temperatures).
			if (energy_eq == true){
				hydro_term_x = -1*0.5*_beta*(_dt)*(_T(i_idx, j_idx) + _T(i_idx + 1, j_idx))*_gx;
				hydro_term_y = -1*0.5*_beta*(_dt)*(_T(i_idx, j_idx) + _T(i_idx, j_idx + 1))*_gy;
			}	
			else {
			hydro_term_x = _dt*_gx;
			hydro_term_y = _dt*_gy;
			}
			_F(i_idx, j_idx) = _U(i_idx, j_idx) + _dt*(_nu*discretization.diffusion(_U, i_idx, j_idx) - discretization.convection_U(_U, _V, i_idx, j_idx)) + hydro_term_x;
			//std::cout<<_F(i_idx, j_idx)<<std::endl;
			_G(i_idx, j_idx) = _V(i_idx, j_idx) + _dt*(_nu*discretization.diffusion(_V, i_idx, j_idx) - discretization.convection_V(_U, _V, i_idx, j_idx)) + hydro_term_y;
			//std::cout<<_G(i_idx, j_idx)<<std::endl;
	}
}
	

//Added: Function implemented for third task.
void Fields::calculate_rs(Grid &grid) {
	int i_idx{0};
	int j_idx{0};

	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		_RS(i_idx, j_idx) = (1.0/_dt) * ((_F(i_idx, j_idx) - _F(i_idx - 1, j_idx))/grid.dx() + (_G(i_idx, j_idx) - _G(i_idx, j_idx - 1))/grid.dy());	
	}
}

//Added: Function implemented for sixth task.
void Fields::calculate_velocities(Grid &grid) {
	int i_idx{0};
	int j_idx{0};
	auto imax = grid.imax();
	auto jmax = grid.jmax();
	
	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		
		if(i_idx != imax){
		//_U(i_idx, j_idx) = _F(i_idx, j_idx) - _dt*dPdx(i_idx, j_idx, grid);
		_U(i_idx, j_idx) = _F(i_idx, j_idx) - _dt / grid.dx() * (_P(i_idx + 1, j_idx) - _P(i_idx, j_idx));
		}
		if(j_idx != jmax){
		//_V(i_idx, j_idx) = _G(i_idx, j_idx) - _dt*dPdy(i_idx, j_idx, grid);
		_V(i_idx, j_idx) = _G(i_idx, j_idx) - _dt / grid.dx() * (_P(i_idx, j_idx + 1) - _P(i_idx, j_idx));
		}
	}
}

void Fields::calculate_temperatures(Grid &grid, Discretization &discretization){
	Matrix<double> new_T(_imax + 2, _jmax + 2);
	new_T = _T;
	int i_idx{0};
	int j_idx{0};
	
	for (const auto currentCell : grid.fluid_cells()){
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		new_T(i_idx, j_idx) = _T(i_idx, j_idx) + _dt*(_alpha*discretization.diffusion(_T, i_idx, j_idx) - discretization.convection_T(_T, _U, _V, i_idx, j_idx));
		//std::cout<<new_T(i_idx, j_idx)<<std::endl;
	}
	_T = new_T;
}

//Added: Function implemented for seventh task.
double Fields::calculate_dt(Grid &grid, bool energy_eq) { 

	int i_idx{0};
	int j_idx{0};
	double max_u{0};
	double max_v{0};
	//for (const auto currentCell : grid.fluid_cells()){
	//	i_idx = currentCell->i();
	//	j_idx = currentCell->j();
		
	//	if (abs(_U(i_idx, j_idx)) > max_u){
	//		max_u = abs(_U(i_idx, j_idx));
	//	}
	//	else { max_u = max_u;}
		
	//	if (abs(_V(i_idx, j_idx)) > max_v){
	//		max_v = abs(_V(i_idx, j_idx));
	//	}
	//	else { max_v = max_v;}
	//}
	max_u = std::abs(*std::max_element(_U.data() + 1, _U.data() + grid.fluid_cells().size() + 1));
        max_v = std::abs(*std::max_element(_V.data() + 1, _V.data() + grid.fluid_cells().size() + 1));
        std::cout<<"umax = "<<max_u<<std::endl;
        std::cout<<"vmax = "<<max_v<<std::endl;
        std::cout<<"dx = "<<grid.dx()<<std::endl;
        std::cout<<"dy = "<<grid.dy()<<std::endl;
	
	double dta = std::min(grid.dx()/max_u, grid.dy()/max_v);
	double dtb = pow((1/pow(grid.dx(), 2.0) + 1/pow(grid.dy(), 2.0)), -1.0)/(2*_nu);
	double dt = _tau*std::min(dta, dtb);
	std::cout<<"dta = "<<dta<<std::endl;
	std::cout<<"dtb = "<<dtb<<std::endl;
	
	if (energy_eq == true){
		double dtc1 = 1/pow(grid.dx(), 2.0) + 1/pow(grid.dy(), 2.0);
		double dtc = 1/(dtc1 * 2 * _alpha);
		dt = std::min(dt, dtc);
		std::cout<<"dtc = "<<dtc<<std::endl;
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

double Fields::dt() const { return _dt; }

double Fields::dPdx(int i_idx, int j_idx, Grid &grid){
	return (1/grid.dx()) * (_P(i_idx + 1, j_idx) - _P(i_idx, j_idx));
}

double Fields::dPdy(int i_idx, int j_idx, Grid &grid){
	return (1/grid.dy()) * (_P(i_idx, j_idx + 1) - _P(i_idx, j_idx));
}
