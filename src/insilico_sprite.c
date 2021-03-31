#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../settings.h"
#include "../headers/utilities.h"
#include "../headers/io.h"
#include "../headers/clustering.h"
#include "../headers/insilico_sprite.h"

void TAGGING(double p_tagging, double p_detection, struct cluster *cluster, double **contact_matrix)
// This function implements a proxy of the split-pool tagging process over a polymer structure
{
    struct cluster *tagged_cluster = random_select(p_tagging * p_detection, cluster);
    struct cluster *unique_cluster = unique_bins(tagged_cluster);
    if(unique_cluster->size <= 1)
    {
        free_cluster(tagged_cluster);
        free_cluster(unique_cluster);
        return;
    }
    struct bead *unique_array = cluster_to_array(unique_cluster);
    double weight = 2. / tagged_cluster->size;
    for(int a=0; a<unique_cluster->size; a++)
    {
        for(int b=0; b<a; b++)
        {
            struct bead bead_a = unique_array[a];
            struct bead bead_b = unique_array[b];
            contact_matrix[bead_a.bin][bead_b.bin] = contact_matrix[bead_a.bin][bead_b.bin] + weight;
            contact_matrix[bead_b.bin][bead_a.bin] = contact_matrix[bead_b.bin][bead_a.bin] + weight;
        }
    }
    free(unique_array);
    free_cluster(tagged_cluster);
    free_cluster(unique_cluster);
}

double **SPRITE_SINGLE_STRUCTURE(double p_crosslinking, double p_tagging, double p_detection, char *structure_name)
// This function implements a proxy of a SPRITE experiment over a single polymer structure
{
    double **sprite_matrix = initialize_contact_matrix();
    struct bead *beads_chain = read_structure_data(structure_name);
    struct set *crosslinked_set = CROSSLINKING(p_crosslinking, beads_chain);
    struct set_item *sp = crosslinked_set->start;
    while(sp != NULL)
    {
        TAGGING(p_tagging, p_detection, sp->cluster, sprite_matrix);
        sp = sp->next;
    }
    free(beads_chain);
    free_set(crosslinked_set);
    return sprite_matrix;
}

double **SPRITE_EXPERIMENT(int N_cells, double p_crosslinking, double p_tagging, double p_detection)
// This function implements a proxy of a SPRITE experiment over a population of in-silico cells
{
    double **sprite_matrix = initialize_contact_matrix();
    int structure_index = 1 + (int) floor(drand48() * N_structures);
    for(int s=0; s<2*N_cells; s++)
    {
        char structure_name[100];
        snprintf(structure_name, 100, "../data/structures/structure_%d.txt", structure_index);
        double **single_structure_sprite_matrix = SPRITE_SINGLE_STRUCTURE(p_crosslinking, p_tagging, p_detection, structure_name);
        update_contact_matrix(sprite_matrix, single_structure_sprite_matrix);
        free_contact_matrix(single_structure_sprite_matrix);
        update_structure_index(&structure_index);
    }
    return sprite_matrix;
}
