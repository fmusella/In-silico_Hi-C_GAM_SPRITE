#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../settings.h"
#include "../headers/utilities.h"
#include "../headers/io.h"
#include "../headers/clustering.h"
#include "../headers/insilico_hic.h"

struct set *BIOTINYLATION(double p_biotinylation, struct set *crosslinked_set)
// This function implements a proxy of the biotinylation process on the polymers
{
    struct set *biotinylated_set = initialize_set();
    struct set_item *sp = crosslinked_set->start;
    while(sp != NULL)
    {
        struct cluster *biotinylated_cluster = random_select(p_biotinylation, sp->cluster);
        if(biotinylated_cluster->size > 1) add_to_set(biotinylated_set, biotinylated_cluster);
        else free_cluster(biotinylated_cluster);
        sp = sp->next;
    }
    return biotinylated_set;
}

void shuffle(struct bead *array, size_t n)
{
    if(n > 1)
    {
        size_t i;
        for(i = n-1; i > 0; i--)
        {
            size_t j = (unsigned int) (drand48() * (i+1));
            struct bead temp = array[j];
            array[j] = array[i];
            array[i] = temp;
        }
    }
}

void LIGATION(double p_ligation, double p_detection, struct cluster *cluster, double **contact_matrix)
// This function implements a proxy of the ligation process on the polymers
{
    struct cluster *can_ligate_cluster = random_select(p_ligation * p_detection, cluster);
    if (can_ligate_cluster->size <= 1)
    {
        free_cluster(can_ligate_cluster);
        return;
    }
    struct bead *can_ligate_array = cluster_to_array(can_ligate_cluster);
    shuffle(can_ligate_array, can_ligate_cluster->size);
    int *ligation_label = my_calloc(can_ligate_cluster->size, sizeof(int));
    for(int a=0; a<can_ligate_cluster->size; a++)
    {
        if(ligation_label[a] == 1) continue;
        for(int b=0; b<can_ligate_cluster->size; b++)
        {
            struct bead bead_a = can_ligate_array[a];
            struct bead bead_b = can_ligate_array[b];
            if(a != b && ligation_label[b] == 0 && distance(bead_a.coordinates, bead_b.coordinates) < max_distance_to_link)
            {
                ligation_label[a] = 1;
                ligation_label[b] = 1;
                if(bead_a.bin != bead_b.bin)
                {
                    contact_matrix[bead_a.bin][bead_b.bin] = contact_matrix[bead_a.bin][bead_b.bin] + 1.;
                    contact_matrix[bead_b.bin][bead_a.bin] = contact_matrix[bead_b.bin][bead_a.bin] + 1.;
                }
                break;
            }
        }
    }
    free(can_ligate_array);
    free(ligation_label);
    free_cluster(can_ligate_cluster);
}

double **HIC_SINGLE_STRUCTURE(double p_crosslinking, double p_biotinylation, double p_ligation, double p_detection, char *structure_name)
// This function implements a proxy of an Hi-C experiment over a single polymer structure
{
    double **hic_matrix = initialize_contact_matrix();
    struct bead *beads_chain = read_structure_data(structure_name);
    struct set *crosslinked_set = CROSSLINKING(p_crosslinking, beads_chain);
    struct set *biotinylated_set = BIOTINYLATION(p_biotinylation, crosslinked_set);
    struct set_item *sp = biotinylated_set->start;
    while(sp != NULL)
    {
        LIGATION(p_ligation, p_detection, sp->cluster, hic_matrix);
        sp = sp->next;
    }
    free(beads_chain);
    free_set(crosslinked_set);
    free_set(biotinylated_set);
    return hic_matrix;
}

double **HIC_EXPERIMENT(int N_cells, double p_crosslinking, double p_biotinylation, double p_ligation, double p_detection)
// This function implements a proxy of an Hi-C experiment over a population of simulated cells
{
    double **hic_matrix = initialize_contact_matrix();
    int structure_index = 1 + (int) floor(drand48() * N_structures);
    for(int s=0; s<2*N_cells; s++)
    {
        char structure_name[100];
        snprintf(structure_name, 100, "../data/structures/structure_%d.txt", structure_index);
        double **single_structure_hic_matrix = HIC_SINGLE_STRUCTURE(p_crosslinking, p_biotinylation, p_ligation, p_detection, structure_name);
        update_contact_matrix(hic_matrix, single_structure_hic_matrix);
        free_contact_matrix(single_structure_hic_matrix);
        update_structure_index(&structure_index);
    }
    return hic_matrix;
}
