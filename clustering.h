#include <stdio.h>
#include <stdlib.h>
#include "../settings.h"
#include "utilities.h"
#include "io.h"

// This header file contains the functions for the handling of clusters (lists of beads) and sets (lists of clusters)

#ifndef CLUSTERING_HEADER_H
#define CLUSTERING_HEADER_H

// These functions define the data structures cluster and set

struct cluster_item
{
    struct bead bead;
    struct cluster_item *next;
};

struct cluster
{
    int size;
    struct cluster_item *start;
};

struct set_item
{
    struct cluster *cluster;
    struct set_item *next;
};

struct set
{
    int size;
    struct set_item *start;
};

struct cluster *initialize_cluster();

void add_to_cluster(struct cluster *cluster, struct bead bead);
    
struct cluster *random_select(double prob, struct cluster *cluster);

struct cluster *unique_bins(struct cluster *cluster);

struct bead *cluster_to_array(struct cluster *cluster);

struct set *initialize_set();

void add_to_set(struct set *set, struct cluster *cluster);

void free_cluster(struct cluster *cluster);

void free_set(struct set *set);

// This function implements a proxy of the crosslinking process for the SBS model structures
struct set *CROSSLINKING(double prob, struct bead *beads_chain);

#endif
