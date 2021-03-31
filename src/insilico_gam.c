#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../settings.h"
#include "../headers/utilities.h"
#include "../headers/io.h"
#include "../headers/clustering.h"
#include "../headers/insilico_gam_utilities.h"
#include "../headers/insilico_gam.h"

void NP_CUTTING(struct cluster *segregated_beads, char *pair_name)
// This function implements a proxy of the process of cutting a Nuclear Profile on a simulated cell
{
    struct bead *chain_1 = my_malloc(beads_number * sizeof(struct bead));
    struct bead *chain_2 = my_malloc(beads_number * sizeof(struct bead));
    read_pair_data(chain_1, chain_2, pair_name);
    perform_NP_cut(segregated_beads, chain_1, chain_2);
    free(chain_1);
    free(chain_2);
}

struct cluster *GAM_SINGLE_TUBE(double p_detection, char *pair_name)
// This function implements a proxy of a GAM experiment with a single tube
{
    struct cluster *segregated_beads = initialize_cluster();
    NP_CUTTING(segregated_beads, pair_name);
    struct cluster *detected_beads = random_select(p_detection, segregated_beads);
    struct cluster *unique_beads = unique_bins(detected_beads);
    free_cluster(segregated_beads);
    free_cluster(detected_beads);
    return unique_beads;
}

void fill_segregation_table(int tube_index, double **segregation_table, struct cluster *tube_beads)
{
    if(tube_beads->size == 0) return;
    struct cluster_item *p = tube_beads->start;
    while(p != NULL)
    {
        segregation_table[p->bead.bin][tube_index] = 1.;
        p = p->next;
    }
    free(p);
}

double **GAM_EXPERIMENT(int N_tubes, double p_detection)
// This function implements a proxy of a GAM experiment over a population of in-silico cells
{
    double **segregation_table = initialize_segregation_table(N_tubes);
    int pair_index = 1 + (int) floor(drand48() * (double) N_structures / 2.);
    for(int tube_index=0; tube_index<N_tubes; tube_index++)
    {
        char pair_name[100];
        snprintf(pair_name, 100, "../data/pairs/pair_%d.txt", pair_index);
        struct cluster *tube_beads = GAM_SINGLE_TUBE(p_detection, pair_name);
        fill_segregation_table(tube_index, segregation_table, tube_beads);
        free_cluster(tube_beads);
        update_pair_index(&pair_index);
    }
    double **cosegregation_matrix = compute_cosegregation_matrix(segregation_table, N_tubes);
    free_segregation_table(segregation_table);
    return cosegregation_matrix;
}
