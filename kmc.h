#ifndef KMC_H
#define KMC_H

#include "poscar_reader.h"
#include "structure.h"
#include "xdatcar_writer.h"

typedef struct {
    double rate;
    double barrier;
    double frequency;
    int old_index;
    int new_index;
    double time_step;
    double dist;
    double displacement[3];
    char fromA[MAX_TYPE_LENGTH];
    char toB[MAX_TYPE_LENGTH];
} event;

typedef struct {
    double barrier_ev;
    double distance;
    double frequency;
    char fromA[MAX_TYPE_LENGTH];
    char toB[MAX_TYPE_LENGTH];
} barrier;

event select_event(int *occupation, atomic_site *sites, barrier *barrier_list, int num_sites, int num_barriers, double temperature);
int get_barrier(barrier *barrier_list, double dist, int num_barriers, const char *fromA, const char *toB);
barrier* build_barrier_list(const char *filename, int *num_barriers);
void calculate_D_coef(const double x[], const double y[], int n, double *slope, double *intercept, double *r_value, double *stderr_val);


#endif
