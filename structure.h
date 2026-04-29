#ifndef STRUCTURE_H
#define STRUCTURE_H
#include "poscar_reader.h"

typedef struct {
    int initial_index;
    int current_index;
    double displacements[3];
} Atom;

void initialize_occupation(int* occupation, Topology *crystall);
int count_atoms(int* occupation, int n_sites);
void initialize_atoms(int* occupation, Atom* atoms, atomic_site* sites, int n_sites);
void move_atom(int* occupation, Atom* atoms, int num_atoms, int old_index, int new_index, double *displacement);
double calculate_msd(Atom* atoms, int num_atoms);

#endif // POSCAR_READER_H
