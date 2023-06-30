#include "LBstreaming.hpp"
#include <iostream>
#include "Enums.hpp"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
    int dx, dy;
    double fi;
    
    for (int y = 1; y < xlength + 1; ++y) {
        for (int x = 1; x < xlength + 1; ++x) {
            for (int i = 0; i < Q; ++i) {

                // dx = c_i_x*dt, where dt = 1, etc.
                dx = LATTICEVELOCITIES[i][0];
                dy = LATTICEVELOCITIES[i][1];
                

                fi = collideField [ Q * ((y-dy)*(xlength+2) + x - dx) + i ];
                streamField [Q * (y*(xlength+2) + x) + i ] = fi;
                
            }
        }
    }
    // std::cout<<"streamField "<<*streamField<<std::endl;
}