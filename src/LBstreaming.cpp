#include "LBstreaming.hpp"
#include "LBDefinitions.hpp"

void doStreaming(double *collideField, double *streamField,int *flagField,int *xlength_lbm){

	unsigned int xylen2 = (xlength_lbm[0] + 2) * (xlength_lbm[1] + 2);

	for(unsigned int z = 1; z <= xlength_lbm[2]; ++z){					// iterate through all FLUID cells
		unsigned int izz = xylen2 * z;											// improve code readability and optimize computations

		for(unsigned int y = 1; y <= xlength_lbm[1]; ++y){
			unsigned int jy = (xlength_lbm[0] + 2) * y;							// improve code readability and optimize computations

			for(unsigned int x = 1; x <= xlength_lbm[0]; ++x){
				unsigned int current_idx = (x + jy + izz) * Q_NUMBER;
				unsigned int current_idx_flag = x + jy + izz;
				for(unsigned int Q_iter = 0; Q_iter < Q_NUMBER; ++Q_iter){	// iterate through all directions
					// this is very cache ineficient, because ~9 different cache lines are accessed for each cell update(3 per z-plane for each row)

					// - find neighboring cells and take their respective directions ... tricky...
					// - we need the opposite cell to the Q - direction(hence the "-" sign in front of the offset with Lattice velocities)
					// with the velocity in the same direction(hence the 2nd Q_iter without minus)

					if (flagField[current_idx_flag] == FLUID)
					streamField[current_idx + Q_iter] = collideField[current_idx -(LATTICEVELOCITIES[Q_iter][0]
                                                                                +  LATTICEVELOCITIES[Q_iter][1] * (xlength_lbm[0] + 2)
                                                                                +  LATTICEVELOCITIES[Q_iter][2] * xylen2) * Q_NUMBER + Q_iter] ;
				}
			}
		}
	}
}