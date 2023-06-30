#pragma once
#include <cmath>

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention, which is:
// 0: fluid, 3: fixed wall, 4: moving wall
namespace LidDrivenCavity {
const int moving_wall_id = 8;
const int fixed_wall_id = 4;
const double wall_velocity = 1.0;
} // namespace LidDrivenCavity

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

namespace border {
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;
} // namespace border

enum class cell_type {
    FLUID,
    INFLOW,
    OUTFLOW,
    ADIABATIC_WALL,
    HOT_WALL,
    COLD_WALL,
    FIXED_WALL,
    MOVING_WALL,
    DEFAULT
};

static const int LATTICEVELOCITIES[9][2] = {
    {0, 0},   // (0, 0)
    {-1, 0},  // (-1, 0)
    {1, 0},   // (1, 0)
    {0, -1},  // (0, -1)
    {0, 1},   // (0, 1)
    {-1, -1}, // (-1, -1)
    {-1, 1},  // (-1, 1)
    {1, -1},  // (1, -1)
    {1, 1}    // (1, 1)
};

static const double LATTICEWEIGHTS[9] = {
    4.0 / 9.0,   // (0, 0)
    1.0 / 9.0,  // (-1, 0)
    1.0 / 9.0,   //(1, 0)
    1.0 / 9.0,// (0, -1)
    1.0 / 9.0,    //(0, 1)
    1.0 / 36.0,// (-1, -1)
    1.0 / 36.0,  // (-1, 1)
    1.0 / 36.0, //(1, -1)
    1.0 / 36.0,  // (1, 1) 
};

static const double C_S = 1.0 / sqrt(3.0);

// static const int FLUID = 0;
// static const int NO_SLIP = 1;
// static const int MOVING_WALL = 2;

static const int D = 2;
static const int Q = 9;
