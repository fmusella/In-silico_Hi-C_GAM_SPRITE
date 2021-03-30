#include <stdio.h>
#include <stdlib.h>
#include "../settings.h"
#include "../headers/utilities.h"
#include "../headers/io.h"
#include "../headers/clustering.h"

struct cluster *initialize_cluster()
{
    struct cluster *cluster = my_malloc(sizeof(struct cluster));
    cluster->size = 0;
    cluster->start = NULL;
    return cluster;
}

void add_to_cluster(struct cluster *cluster, struct bead bead)
{
    struct cluster_item *p = my_malloc(sizeof(struct cluster_item));
    p->bead = bead;
    p->next = cluster->start;
    cluster->start = p;
    cluster->size = cluster->size + 1;
}

struct cluster *random_select(double prob, struct cluster *cluster)
{
    struct cluster *selected_cluster = initialize_cluster();
    struct cluster_item *p = cluster->start;
    while(p != NULL)
    {
        double random_variable = drand48();
        if(random_variable <= prob) add_to_cluster(selected_cluster, p->bead);
        p = p->next;
    }
    return selected_cluster;
}

struct cluster *unique_bins(struct cluster *cluster)
{
    struct cluster *unique_cluster = initialize_cluster();
    struct cluster_item *p = cluster->start;
    while(p != NULL)
    {
        int add = 1;
        struct cluster_item *q = unique_cluster->start;
        while(q != NULL)
        {
            if(q->bead.bin == p->bead.bin)
            {
                add = 0;
                break;
            }
            q = q->next;
        }
        if(add == 1) add_to_cluster(unique_cluster, p->bead);
        p = p->next;
    }
    return unique_cluster;
}

struct bead *cluster_to_array(struct cluster *cluster)
{
    struct bead *array = my_malloc(cluster->size * sizeof(struct bead));
    int i = 0;
    struct cluster_item *p = cluster->start;
    while(p != NULL)
    {
        array[i] = p->bead;
        p = p->next;
        i = i + 1;
    }
    return array;
}

struct set *initialize_set()
{
    struct set *set = my_malloc(sizeof(struct set));
    set->size = 0;
    set->start = NULL;
    return set;
}

void add_to_set(struct set *set, struct cluster *cluster)
{
    struct set_item *sp = my_malloc(sizeof(struct set_item));
    sp->cluster = cluster;
    sp->next = set->start;
    set->start = sp;
    set->size = set->size + 1;
}

void free_cluster(struct cluster *cluster)
{
    struct cluster_item *p = cluster->start;
    while (p != NULL)
    {
        struct cluster_item *q;
        q = p;
        p = p->next;
        free(q);
    }
    free(cluster);
}

void free_set(struct set *set)
{
    struct set_item *sp = set->start;
    while (sp != NULL)
    {
        struct set_item *sq;
        sq = sp;
        sp = sp->next;
        free_cluster(sq->cluster);
    }
    free(set);
}

void fill_cluster(double prob, int a, struct bead *beads_chain, struct cluster *cluster, int *bead_visited_label)
{
    bead_visited_label[a] = 1;
    double random_variable = drand48();
    if (random_variable <= prob) add_to_cluster(cluster, beads_chain[a]);
    for (int b=0; b<beads_number; b++)
    {
        if (a != b && bead_visited_label[b] == 0 && beads_chain[a].color == beads_chain[b].color && distance(beads_chain[a].coordinates, beads_chain[b].coordinates) < max_distance_to_link)
        {
            fill_cluster(prob, b, beads_chain, cluster, bead_visited_label);
        }
    }
}

struct set *CROSSLINKING(double prob, struct bead *beads_chain)
// This function implements a proxy of the crosslinking process for the SBS model structures
{
    struct set *set = initialize_set();
    int *bead_visited_label = my_calloc(beads_number, sizeof(int));
    for (int a=0; a<beads_number; a++)
    {
        if (bead_visited_label[a] == 0 && beads_chain[a].color != 1.)
        {
            struct cluster *cluster = initialize_cluster();
            fill_cluster(prob, a, beads_chain, cluster, bead_visited_label);
            if (cluster->size > 1) add_to_set(set, cluster);
            else free_cluster(cluster);
        }
    }
    free(bead_visited_label);
    return set;
}
