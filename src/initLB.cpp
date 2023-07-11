#include "initLB.hpp"
#include <iostream>


int initialiseFields(double* collideField, double* streamField, int* flagField, int* xlength, std::string problem, int* initxyzXYZ) {
    int x, y, z, i;
    int xlen2 = xlength[0] + 2;
    int ylen2 = xlength[1] + 2;
    int zlen2 = xlength[2] + 2;
    int xyz_field_domain = xlen2 * ylen2 * zlen2 * Q_NUMBER;
    int xylen2 = xlen2 * ylen2;
    int xsize, ysize;

    if (problem == "lbm_plane_shear_flow") {
        // Fluid init (inner part of flagField).
        for (z = 1; z <= xlength[2]; ++z) {
            for (y = 1; y <= xlength[1]; ++y) {
                for (x = 1; x <= xlength[0]; ++x) {
                    flagField[x + y * xlen2 + z * xylen2] = FLUID;
                }
            }
        }
    }
    else if (problem == "lbm_flow_over_step") {
        // Fluid init (inner part of flagField).
        for (z = 1; z <= xlength[2]; ++z) {
            for (y = 1; y <= xlength[1]; ++y) {
                for (x = 1; x <= xlength[0]; ++x) {
                    flagField[x + y * xlen2 + z * xylen2] = FLUID;
                }
            }
        }

        // add "step" initialization - "cube" has z = x = xlen/2, y = ylen ... x-> dimension starts from top left!
        if (2 * xlength[2] < xlength[0]) {
            std::cout << "Error: for the flow over a step scenario, z-dimension has to be greater than xlength/2.\n";
            return 1;
        } else {
            for (z = 1; z <= xlength[0] / 2; ++z) {
                for (y = 1; y <= xlength[1]; ++y) {
                    for (x = 1; x <= xlength[0] / 2; ++x) { // cube is in the 2nd half of the x dimension
                        flagField[x + y * xlen2 + z * xylen2] = NO_SLIP;
                    }
                }
            }
        }
        std::cout << "init done!\n";
    }
    else {
        std::cout << "\n\nUnknown / Misspelled scenario name!\n\n";
        std::cout << "Known scenarios are:\n\n";
        std::cout << "Plane shear flow\n";
        std::cout << "Flow over a step\n";
        return 1;
    }

    // identical for all scenarios => called once outside of if-statement.

    // stream & collide Fields initialization.
    for (x = 0; x < xyz_field_domain; x += Q_NUMBER) {
        for (i = 0; i < Q_NUMBER; ++i) {
            streamField[i + x] = LATTICEWEIGHTS[i];
            collideField[i + x] = LATTICEWEIGHTS[i];
        }
    }

    /** flagField init:
    * 0 - FLUID
    * 1 - NO_SLIP
    * 2 - MOVING_WALL
    * 3 - FREE-SLIP
    * 4 - INFLOW
    * 5 - OUTFLOW
    * 6 - PRESSURE_IN
    **/

    y = 0;
    for (int x = 0; x < xlen2; x++) {
        for (int z = 0; z < zlen2; z++) {
            flagField[x + y * xlen2 + z * xylen2] = FREE_SLIP;
        }
    }

    y = 2;
    for (int x = 0; x < xlen2; x++) {
        for (int z = 0; z < zlen2; z++) {
            flagField[x + y * xlen2 + z * xylen2] = FREE_SLIP;
        }
    }

    /* Boundary initialization: using input parameters */
    for (int y = 0; y < ylen2; ++y) {
        for (x = 0; x < xlen2; ++x) {
            flagField[x + y * xlen2] = initxyzXYZ[2]; // z- dimension
            flagField[x + y * xlen2 + (zlen2 - 1) * xylen2] = initxyzXYZ[5]; // z+ dimension
        }
    }
    for (int z = 0; z < zlen2; ++z) {
        for (x = 0; x < xlen2; ++x) {
            flagField[x + z * xylen2] = initxyzXYZ[1]; // y- dimension
            flagField[x + (ylen2 - 1) * xlen2 + z * xylen2] = initxyzXYZ[4]; // y+ dimension
        }
    }
    for (int z = 0; z < zlen2; ++z) {
        for (y = 0; y < ylen2; ++y) {
            flagField[y * xlen2 + z * xylen2] = initxyzXYZ[0]; // x- dimension
            flagField[xlen2 - 1 + y * xlen2 + z * xylen2] = initxyzXYZ[3]; // x+ dimension
        }
    }

    // Check for forbidden boundary cells (no_slip or free_slip)
    // fluidNeighbors: counter of neighbor cells that are fluid
    // Although we have 3D problems, we assume that boundaries have no fluid neighbors
    // in the y-direction (foreground-background).
    // This check has a big optimization potential but it is executed only once at the start.
    int fluidNeighbors;
    for (z = 1; z <= xlength[2]; ++z) {
        for (y = 1; y <= xlength[1]; ++y) {
            for (x = 1; x <= xlength[0]; ++x) {
                // How many fluid neighbors does the current cell have?
                fluidNeighbors = 0;
                if (flagField[x + y * xlen2 + z * xylen2] != FLUID) {
                    if (flagField[(x - 1) + y * xlen2 + z * xylen2] == FLUID) ++fluidNeighbors;
                    if (flagField[(x + 1) + y * xlen2 + z * xylen2] == FLUID) ++fluidNeighbors;
                    if (flagField[x + y * xlen2 + (z - 1) * xylen2] == FLUID) ++fluidNeighbors;
                    if (flagField[x + y * xlen2 + (z + 1) * xylen2] == FLUID) ++fluidNeighbors;
                }

                // Depending on the number of neighbors, is it an inner cell, an edge, a valid corner or something forbidden?
                if (fluidNeighbors > 2) {
                    std::cout << "The domain contains at least one boundary cell with more than 2 fluid neighbors.\n";
                    std::cout << "Coordinates: x = " << x << ", y = " << y << ", z = " << z << std::endl;
                    return 1;
                } else if (fluidNeighbors == 0) {
                    continue; // it is an inner obstacle cell
                } else if (fluidNeighbors == 1) {
                    continue; // it is an edge
                } else if (fluidNeighbors == 2) {
                    // Check if the two fluid neighbors share a corner (check all the possibilities, remember we don't have general 3D cases)
                    if (flagField[(x - 1) + y * xlen2 + z * xylen2] == FLUID && flagField[x + y * xlen2 + (z + 1) * xylen2] == FLUID) {
                        continue;
                    }
                    if (flagField[(x + 1) + y * xlen2 + z * xylen2] == FLUID && flagField[x + y * xlen2 + (z + 1) * xylen2] == FLUID) {
                        continue;
                    }
                    if (flagField[(x + 1) + y * xlen2 + z * xylen2] == FLUID && flagField[x + y * xlen2 + (z - 1) * xylen2] == FLUID) {
                        continue;
                    }
                    if (flagField[(x - 1) + y * xlen2 + z * xylen2] == FLUID && flagField[x + y * xlen2 + (z - 1) * xylen2] == FLUID) {
                        continue;
                    }
                    // The execution didn't move to the next iteration, so this is not a valid corner!
                    std::cout << "The domain contains at least one boundary cell with maximum 2 fluid neighbors that is forbidden.\n";
                    std::cout << "Coordinates: x = " << x << ", y = " << y << ", z = " << z << std::endl;
                    return 1;
                }
            }
        }
    }
    return 0;
}
