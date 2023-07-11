
#include <math.h>

static const int FLUID = 0;
static const int NO_SLIP = 1;
static const int MOVING_WALL = 2;
static const int FREE_SLIP = 3;
static const int INFLOW = 4;
static const int OUTFLOW = 5;
static const int PRESSURE_IN = 6;
static const int Q_NUMBER = 19;



static const int LATTICEVELOCITIES[Q_NUMBER][3] = {{0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 0}, {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0}, {0, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
static const double LATTICEWEIGHTS[Q_NUMBER] = {1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36, 1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36, 2.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36};

static const double C_S = 0.577350269189626;

	// Precalculations for performance
static const double C_S_sq = 0.577350269189626 * 0.577350269189626;

