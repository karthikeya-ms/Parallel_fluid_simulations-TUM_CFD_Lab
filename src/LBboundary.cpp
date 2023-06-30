#include "LBboundary.hpp"
#include "Enums.hpp"
#include "computeCellValues.hpp"
#include "LBstreaming.hpp"

void treatBoundary(double *collideField, int* flagField,
                   const double velocityWallX, const double velocityWallY, int xlength){
	int x, y, dx, dy, i;
	int len = xlength + 2;
	double cu = 0;
	double density;
	double *currentCell;
	int index;

	for (x = 0; x < len; x++){
        for (y = 0; y < len; y++) {
            index = y*len + x;
            currentCell = collideField + Q*index;
            for (i = 0; i < Q; ++i) {
                // save Lattice velocities for easier manipulation
                dx = LATTICEVELOCITIES[i][0];
                dy = LATTICEVELOCITIES[i][1];

                // check if index is valid (if particle is streamed to fluid)
                if ((x+dx > 0 && x+dx<len-1) && (y+dy > 0 && y+dy< len-1)){
                    // deal with moving wall  boundary condition
                    if (flagField[index] == MOVING_WALL){
                        cu = velocityWallX*LATTICEVELOCITIES[i][0] + velocityWallY*LATTICEVELOCITIES[i][1];
                        computeDensity(currentCell, &density);

                        *(currentCell + i) = *(currentCell + Q*(dy*len + dx) + Q-1-i) +
                                        2*LATTICEWEIGHTS[i]*density*cu/(C_S*C_S);
                    }
                    // deal with no-slip boundary condition
                    else if (flagField[index] == NO_SLIP){
                        *(currentCell + i) = *(currentCell + Q*(dy*len + dx) + Q-1-i );
                    }
                }
            }
		}
	}
}