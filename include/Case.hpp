#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Boundary.hpp"
#include "Discretization.hpp"
#include "Domain.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"

/**
 * @brief Class to hold and orchestrate the simulation flow.
 *
 */
class Case {
  public:
    /**
     * @brief Parallel constructor for the Case.
     *
     * Reads input file, creates Fields, Grid, Boundary, Solver and sets
     * Discretization parameters Creates output directory
     *
     * @param file_name for the input parameters
     * @param argn
     * @param args
     * @param my_rank the rank of the parallel process
     * @param size total number of processes
     */
    Case(std::string file_name, int argn, char **args, int my_rank = 0, int size = 1);

    /**
     * @brief Main function to simulate the flow until the end time.
     *
     * Calculates the fluxes
     * Calculates the right hand side
     * Solves pressure
     * Calculates velocities
     * Outputs the solution files
     */
    void simulate();

  private:
    /// Plain case name without paths
    std::string _case_name;
    /// Output directory name
    std::string _dict_name;
    /// Geometry file name
    std::string _geom_name{"NONE"};
    /// Relative input file path
    std::string _prefix;
    
    bool _energy_eq{false};
    
    int _iproc;
    int _jproc;
    int _my_rank;
    int _size;

    /// Simulation time
    double _t_end;
    /// Solution file outputting frequency
    double _output_freq;

    Fields _field;
    Grid _grid;
    Discretization _discretization;
    std::unique_ptr<PressureSolver> _pressure_solver;
    std::vector<std::unique_ptr<Boundary>> _boundaries;

    /// Solver convergence tolerance
    double _tolerance;

    /// Maximum number of iterations for the solver
    int _max_iter;

    /**
     * @brief Creating file names from given input data file
     *
     * Extracts path of the case file and creates code-readable file names
     * for outputting directory and geometry file.
     *
     * @param[in] input data file
     */
    void set_file_names(std::string file_name);

    /**
     * @brief Solution file outputter
     *
     * Outputs the solution files in .vtk format. Ghost cells are excluded.
     * Pressure is cell variable while velocity is point variable while being
     * interpolated to the cell faces
     *
     * @param[in] Timestep of the solution
     */
    void output_vtk(int t, int my_rank = 0);

    void build_domain(Domain &domain, int imax_domain, int jmax_domain);

    void output_log(std::string dat_file_name, double nu, double UI, double VI, double PI, double GX, double GY,
                    double xlength, double ylength, double dt, double imax, double jmax, double gamma, double omg,
                    double tau, double itermax, double eps, double TI, double alpha, double beta, double num_walls
                    );
    std::ofstream simulation_log_file();
};
