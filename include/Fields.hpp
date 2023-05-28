#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"

/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields
     *
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     *
     */
    Fields(double _nu, double _dt, double _tau, int imax, int jmax, double UI, double VI, double UIN, double VIN, double PI, double TI = 0.0, double alpha = 0.0, double beta = 0.0);

    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     * @param[in] discretization used to help calculate difussion and convection terms
     * @param[in] energy parameter to determine if the temperature is being calculated
     *
     */
    void calculate_fluxes(Grid &grid, Discretization &discretization, bool energy_eq);

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_velocities(Grid &grid);
    
    /**
     * @brief Temperature calculation using previous time step temperature values
     *
     * @param[in] grid in which the calculations are done
     * @param[in] discretization used to help calculate difussion and convection terms
     *
     */
    void calculate_temperatures(Grid &grid, Discretization &discretization);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition
     *
     * @param[in] grid in which the calculations are done
     *
     */
    double calculate_dt(Grid &grid, bool energy_eq);

    /// x-velocity index based access and modify
    double &u(int i, int j);

    /// y-velocity index based access and modify
    double &v(int i, int j);

    /// pressure index based access and modify
    double &p(int i, int j);

    /// RHS index based access and modify
    double &rs(int i, int j);

    /// x-momentum flux index based access and modify
    double &f(int i, int j);

    /// y-momentum flux index based access and modify
    double &g(int i, int j);
    
    /// temperature index based access and modify
    double &t(int i, int j);

    /// get timestep size
    double dt() const;

    /// pressure matrix access and modify
    Matrix<double> &p_matrix();

    //Added: Headers of helper functions containing derivative terms to compute pressure derivatives.
    double dPdx(int i_idx, int j_idx, Grid &grid);
    double dPdy(int i_idx, int j_idx, Grid &grid);


  private:
    /// x-velocity matrix
    Matrix<double> _U;
    /// y-velocity matrix
    Matrix<double> _V;
    /// pressure matrix
    Matrix<double> _P;
    /// x-momentum flux matrix
    Matrix<double> _F;
    /// y-momentum flux matrix
    Matrix<double> _G;
    /// right hand side matrix
    Matrix<double> _RS;
    //  temperature matrix
    Matrix<double> _T;

    /// kinematic viscosity
    double _nu;
    /// gravitional acceleration in x direction
    double _gx{0.0};
    /// gravitional acceleration in y direction
    double _gy{0.0};
    /// timestep size
    double _dt;
    /// adaptive timestep coefficient
    double _tau;
    int _imax;
    int _jmax;
    double _beta;
    double _alpha;
};
