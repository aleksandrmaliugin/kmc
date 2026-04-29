#include "poscar_reader.h"
#include "structure.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void start_writing_xdatcar(const char *filename, Topology *structure, atomic_site *sites, Atom *atoms, int num_atoms)
{
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error! %s.\n", filename);
        return;
    }

    fprintf(file, "KMC simulation\n");
    fprintf(file, "%lf\n", structure->scale_factor);
    for (int i = 0; i < 3; i++) {
        fprintf(file, "%15.10f %15.10f %15.10f\n", structure->lattice[i][0], structure->lattice[i][1], structure->lattice[i][2]);
    }

    for (int i = 0; i < structure->num_types; i++) {
        fprintf(file, "%s ", structure->types[i]);
    }

    fprintf(file, "\n");

    for (int i = 0; i < structure->num_types; i++) {
        fprintf(file, "%d ", structure->count_types[i] - (int)round(structure->vacancy_concentration[i] * structure->count_types[i] ) );
    }



    fprintf(file, "\n");

    fprintf(file, "Direct configuration=     1\n");
    

    for (int i = 0; i < num_atoms; i++)
    {
        int index = atoms[i].initial_index;
        fprintf(file, "%15.10f %15.10f %15.10f\n", sites[index].direct_coordinates[0], sites[index].direct_coordinates[1], sites[index].direct_coordinates[2]);
    }

    fclose(file);
}

void write_config(const char *filename, int* counter_xdatcar, atomic_site *sites, Atom *atoms, int num_atoms )
{
    FILE *file = fopen(filename, "a");
    if (!file) {
        fprintf(stderr, "Error! %s.\n", filename);
        return;
    }

    fprintf(file, "Direct configuration=     %d\n", *counter_xdatcar);

    for (int i = 0; i < num_atoms; i++)
    {
        int index = atoms[i].current_index;
        fprintf(file, "%15.10f %15.10f %15.10f\n", sites[index].direct_coordinates[0], sites[index].direct_coordinates[1], sites[index].direct_coordinates[2]);
    }

    (*counter_xdatcar)++;
    fclose(file);
}