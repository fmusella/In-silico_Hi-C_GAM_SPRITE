#include <stdlib.h>
#include <math.h>
#include "../settings.h"

// This header file contains utility functions used throughout the codes

#ifndef UTILITIES_HEADER_H
#define UTILITIES_HEADER_H

void *my_malloc(size_t size);

void *my_calloc(size_t count, size_t size);

double distance(double x[3], double y[3]);

void update_structure_index(int *structure_index);

void update_pair_index(int *pair_index);

#endif
