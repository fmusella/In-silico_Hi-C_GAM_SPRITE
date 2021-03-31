#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../headers/insilico_hic.h"
#include "../headers/insilico_sprite.h"
#include "../headers/insilico_gam.h"

double **perform_experiment(int which_experiment, int N, double efficiency)
// This function performs the in-silico experiment over a population of SBS polymer structures
{
    double **output;
    switch(which_experiment)
    {
        // The in-silico efficiencies of all processes are interchangeable
        case 1:
            output = HIC_EXPERIMENT(N, efficiency, 1., 1., 1.);
            break;
        case 2:
            output = SPRITE_EXPERIMENT(N, efficiency, 1., 1.);
            break;
        case 3:
            output = GAM_EXPERIMENT(N, efficiency);
            break;
    }
    return output;
}

void write_output_name(char output_name[100], int which_experiment, int N, double efficiency)
{
    char efficiency_string[100];
    snprintf(efficiency_string, 100, "%lf", efficiency);
    char info_to_print[100];
    snprintf(info_to_print, 100, "%d_%c%c%c", N, efficiency_string[0], efficiency_string[2], efficiency_string[3]);
    switch(which_experiment)
    {
        case 1:
            snprintf(output_name, 100, "hic_mat_%s.txt", info_to_print);
            break;
        case 2:
            snprintf(output_name, 100, "sprite_mat_%s.txt", info_to_print);
            break;
        case 3:
            snprintf(output_name, 100, "gam_mat_%s.txt", info_to_print);
            break;
    }
}

void print_info(int which_experiment, int N, double efficiency)
{
    printf("\n\n");
    switch(which_experiment)
    {
        case 1:
            printf("In-silico Hi-C experiment.\n");
            printf("Number of in-silico cells: %d,\n", N);
            break;
        case 2:
            printf("In-silico SPRITE experiment.\n");
            printf("Number of in-silico cells: %d,\n", N);
            break;
        case 3:
            printf("In-silico GAM experiment.\n");
            printf("Number of in-silico cells: %d,\n", N);
            break;
    }
    printf("Efficiency: %lf.\n", efficiency);
    printf("\n\n");
}

void save_output(int which_experiment, double **output, int N, char output_name[100])
{
    switch(which_experiment)
    {
        case 1:
            save_contact_matrix(output, output_name);
            break;
        case 2:
            save_contact_matrix(output, output_name);
            break;
        case 3:
            save_contact_matrix(output, output_name);
            break;
    }
}

void free_output(int which_experiment, double **output)
{
    switch(which_experiment)
    {
        case 1:
            free_contact_matrix(output);
            break;
        case 2:
            free_contact_matrix(output);
            break;
        case 3:
            free_contact_matrix(output);
            break;
    }
}


/**
 * argc: 3.
 * @argv[1]: which_experiment (int). 1 for In-silico Hi-C,
 *                                   2 for In-silico SPRITE,
 *                                   3 for In-silico GAM,
 * @argv[2]: N (int). Number of in-silico cells
 * @argv[3]: efficiency (double).
 */
int main(int argc, char *argv[])
{
    srand48(time(NULL));

    int which_experiment;
    sscanf(argv[1], "%d", &which_experiment);
    int N;
    sscanf(argv[2], "%d", &N);
    double efficiency;
    sscanf(argv[3], "%lf", &efficiency);

    print_info(which_experiment, N, efficiency);

    double **output = perform_experiment(which_experiment, N, efficiency);
    
    char output_name[100];
    write_output_name(output_name, which_experiment, N, efficiency);
    save_output(which_experiment, output, N, output_name);
    
    free_output(which_experiment, output);
    
    return 0;
}
