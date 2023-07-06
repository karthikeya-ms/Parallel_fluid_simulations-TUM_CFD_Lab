

// static const int FLUID = 0;
// static const int NO_SLIP = 1;
// static const int MOVING_WALL = 2;
/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(double *collideField, double *streamField,int *flagField,int xlength);