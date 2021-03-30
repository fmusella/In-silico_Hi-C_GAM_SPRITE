#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../settings.h"
#include "../headers/utilities.h"

void *my_malloc(size_t size)
{
    void *ptr = NULL;
    while (ptr == NULL) ptr = malloc(size);
    return ptr;
}

void *my_calloc(size_t count, size_t size)
{
    void *ptr = NULL;
    while (ptr == NULL) ptr = calloc(count, size);
    return ptr;
}

double distance(double x[3], double y[3])
{
    int i;
    double d = 0.;
    for (i = 0; i < 3; i++) d = d + pow(x[i] - y[i], 2);
    d = sqrt(d);
    return d;
}

void update_structure_index(int *structure_index)
{
    if(*structure_index == N_structures) *structure_index = 1; else *structure_index = *structure_index + 1;
}

void update_cell_index(int *cell_index)
{
    if(*cell_index == (int) N_structures / 2) *cell_index = 1; else *cell_index = *cell_index + 1;
}
