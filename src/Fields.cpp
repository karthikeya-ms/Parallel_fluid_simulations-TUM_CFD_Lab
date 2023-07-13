#include "Fields.hpp"

#include<cmath>
#include <algorithm>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, double alpha, double beta, int imax, int jmax, double UI, double VI, double PI, double TI, double GX, double GY, Grid &grid, std::string energy_eq)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _energy_eq(energy_eq) {
    // intializing u, v and p
    _U = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _V = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _P = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _T = Matrix<double>(imax + 2, jmax + 2, 0.0);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);

    _gx = GX;
    _gy = GY;

    for(auto &currentCell: grid.fluid_cells()){
        int i = currentCell->i();
        int j = currentCell->j();
        setu(i,j,UI);
        setv(i,j,VI);
        setp(i,j,PI);
        if(_energy_eq == "on"){
            setT(i,j,TI);
        }
    }
    if(_energy_eq == "on"){
        for(auto &currentCell: grid.inflow_cells()){
            setT(currentCell->i(),currentCell->j(),TI);
        }
    }
}//end of Fields constructor

void Fields::calculate_fluxes(Grid &grid) {

    for(auto &currentCell: grid.fluid_cells()){
        int i = currentCell->i();
        int j = currentCell->j();
        if(i != grid.imax()){   //excluding imax as f_imax is part of the fixed boundary and set as 0.0
            setf(i,j,u(i,j)+dt()*(_nu*(Discretization::diffusion(_U, i, j))-Discretization::convection_u(_U, _V, i, j) - 0.5*_beta*_gx*(T(i,j)+T(i+1,j)) ));
        }
        if(j != grid.jmax()){   // excluding jmax as g_jmax is part of the moving boundary and set as 0.0
            setg(i,j,v(i,j)+dt()*(_nu*(Discretization::diffusion(_V, i, j))-Discretization::convection_v(_U, _V, i, j) - 0.5*_beta*_gy*(T(i,j)+T(i,j+1)) ));
        }
    }
}

void Fields::calculate_rs(Grid &grid) {
    for(auto &currentCell: grid.fluid_cells()){
        int i = currentCell->i();
        int j = currentCell->j();
        // calculating the rhs of pressure equation
        double val = ((f(i,j)-f(i-1,j))/grid.dx() + (g(i,j)-g(i,j-1))/grid.dy())/dt();
        setrs(i, j, val);        
    }
}

void Fields::calculate_velocities(Grid &grid) {
    for(auto &currentCell: grid.fluid_cells()){
        int i = currentCell->i();
        int j = currentCell->j();
        if(i != grid.imax()){   //excluding imax as u_imax is part of the fixed boundary and set as 0.0
            setu(i, j, f(i,j) - _dt/grid.dx()*(p(i+1,j)-p(i,j)));
        }
        if(j != grid.jmax()){   // excluding jmax as v_jmax is part of the moving boundary (in x direction) and set as 0.0
            setv(i, j, g(i,j) - _dt/grid.dy()*(p(i,j+1)-p(i,j)));       
        }
    }
}

//calculating temperature at new timestep
void Fields::calculate_Temperature(Grid &grid){
    Matrix<double> tempT = _T;
    for(auto &currentCell: grid.fluid_cells()){
        int i = currentCell->i();
        int j = currentCell->j();
        tempT(i,j) = _T(i,j) + _dt*(_alpha * (Discretization::diffusion(_T,i,j)) - Discretization::convection_T(_T,_U,_V,i,j));
    }
    _T = tempT;
}

//calculating the timestep keeping in mind the stability crtiteria
double Fields::calculate_dt(Grid &grid) { 
    double dt;
    if((_nu < _alpha) && (_energy_eq != "on")){
        dt = 0.5/_alpha/(1/grid.dx()/grid.dx() + 1/grid.dy()/grid.dy())*_tau;
    }else{
        dt = 0.5/_nu/(1/grid.dx()/grid.dx() + 1/grid.dy()/grid.dy())*_tau;
    }
    
    for(auto &currentCell: grid.fluid_cells()){
        int i = currentCell->i();
        int j = currentCell->j();
        double dt1 = std::abs(grid.dx()/u(i,j))*_tau;
        double dt2 = std::abs(grid.dy()/v(i,j))*_tau;
        if(dt > dt1){  
            dt = dt1;
        }
        if(dt > dt2){  
            dt = dt2;
        }
    }
    _dt = dt;
    return dt; 
}

// get functions
double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }
double &Fields::T(int i, int j) { return _T(i, j); }
std::string &Fields::Energy() { return _energy_eq; }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }

// set functions
void Fields::setp(int i, int j, double val) { _P(i, j) = val; }
void Fields::setT(int i, int j, double val) { _T(i, j) = val; }
void Fields::setu(int i, int j, double val) { _U(i, j) = val; }
void Fields::setv(int i, int j, double val) { _V(i, j) = val; }
void Fields::setf(int i, int j, double val) { _F(i, j) = val; }
void Fields::setg(int i, int j, double val) { _G(i, j) = val; }
void Fields::setrs(int i, int j, double val) { _RS(i, j) = val; }