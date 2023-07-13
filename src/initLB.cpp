#include "initLB.hpp"
#include <iostream>


void initializeFields(double *collideField, double *streamField, int *flagField, int xlength_lbm){
	// initialization of particle distribution func fields
	int xlen = xlength_lbm + 2;
	int xlen2 = xlen*xlen;

	for (int a=0; a<=xlength_lbm+1; a++){
		for (int b=0; b<=xlength_lbm+1; b++){
			for (int c=0; c<=xlength_lbm+1; c++){
				for (int i=0; i<Q_NUMBER; i++){
					// initialize streamField and collideField arrays
					streamField [ Q_NUMBER*(c*xlen2 + b*xlen + a) + i ] = LATTICEWEIGHTS [i];
					collideField [ Q_NUMBER*(c*xlen2 + b*xlen + a) + i ] = LATTICEWEIGHTS [i];
				}
				// set all as inner points
				flagField [ c*xlen2 + b*xlen + a ] = FLUID;
			}
			// if c==xlength_lbm+1, overwrite as moving wall:
			flagField [ xlen2*(xlength_lbm+1) + b*xlen + a ] = MOVING_WALL;
		}
	}
	// overwrite fluid to no-slip at other boundaries:
	for (int k=0; k<xlen; k++){
		for (int j=0; j<xlength_lbm+1; j++){
			// if a||b||c==0 :
			flagField [ j*xlen2 + k*xlen ] = NO_SLIP;
			flagField [ j*xlen2 + k] = NO_SLIP;
			flagField [ k*xlen +j ] = NO_SLIP;

			// if a||b==xlength_lbm+1 :
			flagField [ j*xlen2 + k*xlen + xlen - 1 ] = NO_SLIP;
			flagField [ j*xlen2 + xlen2 - xlen + k ] = NO_SLIP;
		}
		// what remains for (if a||b==0) - case j==xlength_lbm+1:
		flagField [ k*xlen + xlen -1 ] = NO_SLIP;
	}
}
