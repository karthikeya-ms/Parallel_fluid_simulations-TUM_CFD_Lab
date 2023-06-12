#include "Discretization.hpp"

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
	double d2Adx2 = (A(i+1,j) - 2*A(i,j) + A(i-1,j))/pow(_dx,2.0);
	double d2Ady2 = (A(i,j+1) - 2*A(i,j) + A(i,j-1))/pow(_dy,2.0);
	
	return d2Adx2 + d2Ady2;
}

double Discretization::convection_U(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
	double dU2dx = (1/_dx)*(pow((U(i,j) + U(i+1,j))/2.0, 2.0) - pow((U(i-1,j) + U(i,j))/2.0, 2.0)) + (_gamma/_dx)*(abs(U(i,j) + U(i+1,j))/2.0 * (U(i,j) - U(i+1,j))/2.0 - abs(U(i-1,j) + U(i,j))/2.0 * (U(i-1,j) - U(i,j))/2.0);
	double dUVdy = (1/_dy)*((V(i,j) + V(i+1,j))/2.0 * (U(i,j) + U(i,j+1))/2.0 - (V(i,j-1) + V(i+1,j-1))/2.0 * (U(i,j-1) + U(i,j))/2.0) + (_gamma/_dy)*(abs(V(i,j) + V(i+1,j))/2.0 * (U(i,j) - U(i,j+1))/2.0 - abs(V(i,j-1) + V(i+1,j-1))/2 * (U(i,j-1) - U(i,j))/2.0);
	
	return -dU2dx - dUVdy;
}

double Discretization::convection_V(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
	double dV2dy = (1/_dy)*(pow((V(i,j) + V(i,j+1))/2.0, 2.0) - pow((V(i,j-1) + V(i,j))/2.0, 2.0)) + (_gamma/_dy)*((abs(V(i,j) + V(i,j+1))/2.0 * (V(i,j) - V(i,j+1))/2.0) - (abs(V(i,j-1) + V(i,j))/2.0 * (V(i,j-1) - V(i,j))/2.0));
	double dUVdx = (1/_dx)*((U(i,j) + U(i,j+1))/2.0 * (V(i,j) + V(i+1,j))/2.0 - (U(i-1,j) + U(i-1,j+1))/2.0 * (V(i-1,j) + V(i,j))/2.0) + (_gamma/_dx)*(abs(U(i,j) + U(i,j+1))/2.0 * (V(i,j) - V(i+1,j))/2.0 - abs(U(i-1,j) + U(i-1,j+1))/2.0 * (V(i-1,j) - V(i,j))/2.0);
	
	return -dV2dy - dUVdx;
}

double Discretization::convection_T(const Matrix<double> &T, const Matrix<double> &U, const Matrix<double> &V, int i, int j){
	double dUTdx = (1/_dx)*(U(i,j)*(T(i,j) + T(i+1,j))/2.0 - U(i-1,j)*(T(i-1,j) + T(i,j))/2.0) + (_gamma/_dx)*(abs(U(i,j))*(T(i,j) - T(i+1,j))/2.0 - abs(U(i-1,j))*(T(i-1,j) - T(i,j))/2.0);
	double dVTdy = (1/_dy)*(V(i,j)*(T(i,j) + T(i,j+1))/2.0 - V(i,j-1)*(T(i,j-1) + T(i,j))/2.0) + (_gamma/_dy)*(abs(V(i,j))*(T(i,j) - T(i,j+1))/2.0 - abs(V(i,j-1))*(T(i,j-1) - T(i,j))/2.0);
	return dUTdx + dVTdy;
}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {}

double Discretization::gamma(){
    return _gamma;
}
