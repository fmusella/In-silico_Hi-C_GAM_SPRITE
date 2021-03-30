#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../settings.h"
#include "../headers/utilities.h"
#include "../headers/io.h"

struct bead *read_structure_data(char *file_name)
{
    struct bead *chain = my_malloc(beads_number * sizeof(struct bead));
    FILE *file_pointer = fopen(file_name, "r");
    for(int a=0; a<beads_number; a++)
    {
        chain[a].id = a;
        chain[a].bin = (int) floor(a / beads_to_bins);
        fscanf(file_pointer, "%lf %lf %lf %lf", &chain[a].coordinates[0], &chain[a].coordinates[1], &chain[a].coordinates[2], &chain[a].color);
    }
    fclose(file_pointer);
    return chain;
}

void read_cell_data(struct bead *chain_1, struct bead *chain_2, char *file_name)
{
    FILE *file_pointer = fopen(file_name, "r");
    for(int a=0; a<beads_number; a++)
    {
        chain_1[a].id = a;
        chain_2[a].id = a;
        chain_1[a].bin = (int) floor(a / beads_to_bins);
        chain_2[a].bin = (int) floor(a / beads_to_bins);
        fscanf(file_pointer, "%lf %lf %lf %lf %lf %lf %lf %lf", &chain_1[a].coordinates[0], &chain_1[a].coordinates[1], &chain_1[a].coordinates[2], &chain_1[a].color, &chain_2[a].coordinates[0], &chain_2[a].coordinates[1], &chain_2[a].coordinates[2], &chain_2[a].color);
    }
    fclose(file_pointer);
}

double **initialize_contact_matrix()
{
    double **contact_matrix = my_malloc(bins_number * sizeof(double *));
    for(int i=0; i<bins_number; i++)
    {
        contact_matrix[i] = my_malloc(bins_number * sizeof(double));
        for(int j=0; j<bins_number; j++)
        {
            contact_matrix[i][j] = 0.;
        }
    }
    return contact_matrix;
}

void update_contact_matrix(double **matrix_to_update, double **matrix_to_sum)
{
    for(int i=0; i<bins_number; i++)
    {
        for(int j=0; j<bins_number; j++)
        {
            matrix_to_update[i][j] = matrix_to_update[i][j] + matrix_to_sum[i][j];
        }
    }
}

void save_contact_matrix(double **contact_matrix, char *output_file_name)
{
    FILE *output_file_pointer;
    output_file_pointer = fopen(output_file_name, "w");
    for(int i=0; i<bins_number; i++)
    {
        for(int j=0; j<bins_number; j++)
        {
            fprintf(output_file_pointer, "%lf ", contact_matrix[i][j]);
        }
        fprintf(output_file_pointer, "\n");
    }
    fclose(output_file_pointer);
}

void free_contact_matrix(double **contact_matrix)
{
    for(int i=0; i<bins_number; i++)
    {
        free(contact_matrix[i]);
    }
    free(contact_matrix);
}
