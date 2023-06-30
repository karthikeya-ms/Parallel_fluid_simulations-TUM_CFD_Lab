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

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double u2x = ( Discretization::interpolate(U,i,j,1,0)  * Discretization::interpolate(U,i,j,1,0)   - 
                   Discretization::interpolate(U,i,j,-1,0) * Discretization::interpolate(U,i,j,-1,0)   )/_dx 
                    + _gamma * (std::abs(Discretization::interpolate(U,i,j,1,0)  ) * (U(i,j)-U(i+1,j))/2 - 
                                std::abs(Discretization::interpolate(U,i,j,-1,0) ) * (U(i-1,j)-U(i,j))/2)/_dx;
    double uvy = (Discretization::interpolate(U,i,j,0,1)  * Discretization::interpolate(V,i,j,1,0)  - 
                 Discretization::interpolate(V,i,j-1,1,0) * Discretization::interpolate(U,i,j,0,-1))/_dy 
                 +_gamma*(std::abs(Discretization::interpolate(V,i,j,1,0))  * (U(i,j)-U(i,j+1))/2 - 
                          std::abs(Discretization::interpolate(V,i,j-1,1,0))* (U(i,j-1)-U(i,j))/2)/_dy;
    return u2x+uvy;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double v2y = (Discretization::interpolate(V,i,j,0,1) * Discretization::interpolate(V,i,j,0,1) - 
                  Discretization::interpolate(V,i,j,0,-1) * Discretization::interpolate(V,i,j,0,-1))/_dy 
                  + _gamma * (std::abs(Discretization::interpolate(V,i,j,0,1)) * (V(i,j)-V(i,j+1))/2 - 
                              std::abs(Discretization::interpolate(V,i,j,0,-1)) * (V(i,j-1)-V(i,j))/2)/_dy;
    double uvx = (Discretization::interpolate(U,i,j,0,1) * Discretization::interpolate(V,i,j,1,0) - 
                  Discretization::interpolate(U,i-1,j,0,1) * Discretization::interpolate(V,i,j,-1,0))/_dx
                  + _gamma * (std::abs(Discretization::interpolate(U,i,j,0,1)) * (V(i,j)-V(i+1,j))/2 - 
                              std::abs(Discretization::interpolate(U,i-1,j,0,1)) * (V(i-1,j)-V(i,j))/2)/_dx;
    return uvx+v2y;
}

double Discretization::convection_T(const Matrix<double> &T, const Matrix<double> &u, const Matrix<double> &v, int i, int j) {
        double uTx = ( u(i,j) * Discretization::interpolate(T,i,j,1,0) - u(i-1,j) * Discretization::interpolate(T,i,j,-1,0) )/_dx 
                    + _gamma * ( std::abs(u(i,j)) * (T(i,j)-T(i+1,j))/2 - std::abs(u(i-1,j)) * (T(i-1,j)-T(i,j))/2 )/_dx; 
        double vTy = ( v(i,j) * Discretization::interpolate(T,i,j,0,1)- v(i,j-1) * Discretization::interpolate(T,i,j,0,-1) )/_dy 
                    + _gamma * ( std::abs(v(i,j)) * (T(i,j)-T(i,j+1))/2 - std::abs(v(i,j-1)) * (T(i,j-1)-T(i,j))/2 )/_dy; 
        return uTx + vTy;
} 
 
double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
    return laplacian(A, i ,j);
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

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    return 0.5 * ( A(i,j) + A(i+i_offset, j+j_offset) );
}