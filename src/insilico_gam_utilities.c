#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../settings.h"
#include "../headers/utilities.h"
#include "../headers/io.h"
#include "../headers/clustering.h"
#include "../headers/insilico_gam_utilities.h"

struct NP generate_NP()
{
    struct NP NP;

    double v1 = drand48();
    double theta1 = acos(1 - 2 * v1);
    double phi1 = 2 * M_PI * drand48();
    double *x = my_malloc(3 * sizeof(double));
    x[0] = sin(theta1) * cos(phi1);
    x[1] = sin(theta1) * sin(phi1);
    x[2] = cos(theta1);
    NP.unit_vector = x;

    double u = drand48();
    double rho = R * pow(u, 1./3.);
    double v2 = drand48();
    double theta2 = acos(1 - 2 * v2);
    double phi2 = 2 * M_PI * drand48();
    double *y = my_malloc(3 * sizeof(double));
    y[0] = rho * sin(theta2) * cos(phi2);
    y[1] = rho * sin(theta2) * sin(phi2);
    y[2] = rho * cos(theta2);
    NP.point = y;

    return NP;
}

void free_NP(struct NP NP)
{
    free(NP.unit_vector);
    free(NP.point);
}

double distance_bead_NP(struct bead bead, struct NP NP)
{
    double distance = 0.;
    for(int i=0; i<3; i++) distance = distance + (bead.coordinates[i] - NP.point[i]) * NP.unit_vector[i];
    distance = fabs(distance);
    return distance;
}

void perform_NP_cut(struct cluster *segregated_beads, struct bead *chain_1, struct bead *chain_2)
{
    struct NP NP = generate_NP();
    for(int a=0; a<beads_number; a++)
    {
        if(distance_bead_NP(chain_1[a], NP) < h / 2.) add_to_cluster(segregated_beads, chain_1[a]);
        if(distance_bead_NP(chain_2[a], NP) < h / 2.) add_to_cluster(segregated_beads, chain_2[a]);
    }
    free_NP(NP);
}

double **initialize_segregation_table(int N_tubes)
{
    double **segregation_table = my_malloc(bins_number * sizeof(double *));
    for(int i=0; i<bins_number; i++)
    {
        segregation_table[i] = my_malloc(N_tubes * sizeof(double));
        for(int t=0; t<N_tubes; t++)
        {
            segregation_table[i][t] = 0.;
        }
    }
    return segregation_table;
}

void free_segregation_table(double **segregation_table)
{
    for(int i=0; i<bins_number; i++)
    {
        free(segregation_table[i]);
    }
    free(segregation_table);
}

double **compute_cosegregation_matrix(double **segregation_table, int N_tubes)
{
    double **cosegregation_matrix = initialize_contact_matrix();
    for(int i=0; i<bins_number; i++)
    {
        for(int j=i+1; j<bins_number; j++)
        {
            cosegregation_matrix[i][j] = 0.;
            for(int t=0; t<N_tubes; t++) cosegregation_matrix[i][j] = cosegregation_matrix[i][j] + segregation_table[i][t] * segregation_table[j][t];
            cosegregation_matrix[i][j] = cosegregation_matrix[i][j] / (double) N_tubes;
            cosegregation_matrix[j][i] = cosegregation_matrix[i][j];
        }
    }
    return cosegregation_matrix;
}
