#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../settings.h"
#include "utilities.h"
#include "io.h"
#include "clustering.h"

// This header file contains the function that implements in-silico SPRITE on SBS model structures

#ifndef SPRITE_HEADER_H
#define SPRITE_HEADER_H

double **SPRITE_EXPERIMENT(int N_cells, double p_crosslinking, double p_tagging, double p_detection);

#endif
