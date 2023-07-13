#include "LBcollision.hpp"
#include "LBDefinitions.hpp"
#include <iostream>
#include "LBcomputeCellValues.hpp"

void computePostCollisionDistributions(double *currentCell, const double * const tau_lbm, const double *const feq){
        for (int i=0; i<Q_NUMBER; i++)
                *(currentCell+i) = *(currentCell+i) - ( *(currentCell+i)-(*(feq+i)) ) / (*tau_lbm);
}


void doCollision(double *collideField, int *flagField,const double * const tau_lbm,int xlength_lbm){

	double density;
	double velocity[D];
	double feq[Q_NUMBER];
	double *currentCell = NULL; // currentCell points to the first distribution function within the respective cell
        for (int iz=1; iz<=xlength_lbm; iz++){
		for (int iy=1; iy<=xlength_lbm; iy++){
			for (int ix=1; ix<=xlength_lbm; ix++){
				// set pointer to current cell
				currentCell = collideField + Q_NUMBER*(iz*(xlength_lbm+2)*(xlength_lbm+2) + iy*(xlength_lbm+2) + ix);

				// compute density, velocity and equilibrium prob. distrib. for this cell
				computeDensity (currentCell, &density);
				computeVelocity (currentCell, &density, velocity);
				computeFeq (&density, velocity, feq);

				computePostCollisionDistributions (currentCell, tau_lbm, feq);
			}
		}
	}
}

