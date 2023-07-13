#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"

/**
 * @brief Abstact of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
  public:
    /**
     * @brief Main method to patch the boundary conditons to given field and
     * grid
     *
     * @param[in] Field to be applied
     */
    virtual void apply(Fields &field) = 0;
    virtual ~Boundary() = default;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure.
 * When energy equation on, used as Adiabatic Wall.
 */
class FixedWallBoundary : public Boundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature);
    virtual ~FixedWallBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    double _wall_temperature = -1.0;
};

/**
 * @brief Inflow wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is specified through input file, Neumann for pressure.
 */
class InflowBoundary : public Boundary {
  public:
    InflowBoundary(std::vector<Cell *> cells, double _inlet_velocity_x, double _inlet_velocity_y);
    virtual ~InflowBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    double _inlet_velocity_x;
    double _inlet_velocity_y;
};

/**
 * @brief Outflow wall boundary condition for the outer boundaries of the domain.
 * Neumann for velocities, Dirichlet for pressure, which is zero if not specified in the input file.
 */
class OutflowBoundary : public Boundary {
  public:
    OutflowBoundary(std::vector<Cell *> cells, double outflow_pressure);
    virtual ~OutflowBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    double _outflow_pressure;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                       std::map<int, double> wall_temperature);
    virtual ~MovingWallBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_velocity;
    std::map<int, double> _wall_temperature;
};
