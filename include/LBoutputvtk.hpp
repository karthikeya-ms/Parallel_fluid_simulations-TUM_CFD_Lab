#include <string>
#include <stdio.h>
/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by '_case_name' and timestep 't'. 
 */
void writeVtkOutput(const double* const collideField, const int* const flagField, const std::string _case_name, unsigned int t, int xlength_lbm, std::string _dict_name);
