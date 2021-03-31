#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../settings.h"
#include "utilities.h"
#include "io.h"
#include "clustering.h"

// This header file contains the function that implements in-silico Hi-C on SBS model structures

#ifndef HIC_HEADER_H
#define HIC_HEADER_H

double **HIC_EXPERIMENT(int N_cells, double p_crosslinking, double p_biotinylation, double p_ligation, double p_detection);

#endif
