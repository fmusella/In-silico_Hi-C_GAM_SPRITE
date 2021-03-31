#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../settings.h"
#include "utilities.h"
#include "io.h"
#include "clustering.h"
#include "insilico_gam_utilities.h"

// This header file contains the function that implements in-silico GAM on SBS model structures

#ifndef GAM_HEADER_H
#define GAM_HEADER_H

double **GAM_EXPERIMENT(int N_tubes, double p_detection);

#endif
