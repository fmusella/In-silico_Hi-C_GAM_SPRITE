#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../settings.h"
#include "utilities.h"

// This header file contains the functions for the Input-Output (IO) handling

#ifndef IO_HEADER_H
#define IO_HEADER_H

// These functions define the data structure for a polymer chain of the SBS model structures

struct bead
{
    int id;
    int bin;
    double coordinates[3];
    double color;
};

struct bead *read_structure_data(char *file_name);

void read_cell_data(struct bead *chain_1, struct bead *chain_2, char *file_name);

// These functions define the data structure for the contact matrices of in-silico Hi-C, SPRITE and GAM

double **initialize_contact_matrix();

void update_contact_matrix(double **matrix_to_update, double **matrix_to_sum);

void save_contact_matrix(double **contact_matrix, char *output_file_name);

void free_contact_matrix(double **contact_matrix);

#endif
