#include <stdio.h>
#include <stdlib.h>
#include "../settings.h"
#include "utilities.h"
#include "io.h"
#include "clustering.h"

// This header file contains some utility functions necessary for the implementation of in-silico GAM

#ifndef GAM_UTILITIES_HEADER_H
#define GAM_UTILITIES_HEADER_H

// These functions define the data structure for the in-silico Nuclear Profile (NP)

struct NP
{
    double *unit_vector;
    double *point;
};

struct NP generate_NP();

void free_NP(struct NP NP);

double distance_bead_NP(struct bead bead, struct NP NP);

void perform_NP_cut(struct cluster *segregated_beads, struct bead *chain_1, struct bead *chain_2);

// These functions define the Segregation Table and the Co-Segregation Matrix for in-silico GAM

double **initialize_segregation_table(int N_tubes);

void free_segregation_table(double **segregation_table);

double **compute_cosegregation_matrix(double **segregation_table, int N_tubes);

#endif
