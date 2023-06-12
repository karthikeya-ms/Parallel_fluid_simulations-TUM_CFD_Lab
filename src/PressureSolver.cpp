#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries, Communication &communication) {

    double dx = grid.dx();
    double dy = grid.dy();
    _p = Matrix<double>(imax + 2, jmax + 2, 0);

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
        _p(i, j) = field.p(i, j); 
    }
    communication.Communicate(&_p); //Communicate and exchange the pressure values at the ghost layers

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    {
        res = rloc / (grid.fluid_cells().size());
        res = std::sqrt(res);
    }

    return res;
}


