#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "structure.h"
#include "poscar_reader.h"
#include <math.h>

void initialize_occupation_old(int* occupation, Topology *crystall)
{
    //unsigned int seed = (unsigned int)time(NULL);

    for (int i = 0; i < crystall->num_types; i++) {
        printf("Vacancy concentration for %s: ", crystall->types[i]);
        if (scanf("%lf", &crystall->vacancy_concentration[i]) != 1) {
            fprintf(stderr, "Invalid input for vacancy concentration!\n");
            exit(EXIT_FAILURE);
        }
    }

    int count = 0;

    for (int i = 0; i < crystall->num_types; i++) {
        
        for (int j = 0; j < crystall->count_types[i]; j++) {

            //double randomValue = (double)rand_r(&seed) / (RAND_MAX + 1.0);
            double randomValue = (double)rand() / (RAND_MAX + 1.0);

            if (randomValue < crystall->vacancy_concentration[i]) {
                occupation[count] = 0;
            } else {
                occupation[count] = 1;
            }
            count++;

        }
    }

}

void initialize_occupation(int* occupation, Topology *crystall)
{
    srand((unsigned int)time(NULL));

    for (int i = 0; i < crystall->num_types; i++) {
        printf("Vacancy concentration for %s: ", crystall->types[i]);
        if (scanf("%lf", &crystall->vacancy_concentration[i]) != 1) {
            fprintf(stderr, "Invalid input for vacancy concentration!\n");
            exit(EXIT_FAILURE); 
        }
    }

    for (int i = 0; i < crystall->num_sites; i++) {
        occupation[i] = 1;
    }

    int count = 0;

    for (int i = 0; i < crystall->num_types; i++) {
        int total_sites_of_type = crystall->count_types[i];
        int num_vacancies = (int)round(crystall->vacancy_concentration[i] * total_sites_of_type);

        int created_vacancies = 0;
        while (created_vacancies < num_vacancies) {
            int random_index = count + (int)((double)rand() / (RAND_MAX + 1.0) * total_sites_of_type);

            if (occupation[random_index] == 1) {
                occupation[random_index] = 0;
                created_vacancies++;
            }
        }

        count += total_sites_of_type;
    }
}

int count_atoms(int* occupation, int n_sites)
{
    int count = 0;

    for (int i = 0; i < n_sites; i++) {
        if (occupation[i] == 1) {
            count++;
        }
    }

    return count;
}

void initialize_atoms(int* occupation, Atom* atoms, atomic_site* sites, int n_sites)
{
    int count = 0;
    for (int i = 0; i < n_sites; i++) {
        if (occupation[i] == 1) {
            atoms[count].initial_index = sites[i].index;
            atoms[count].current_index = sites[i].index;
            atoms[count].displacements[0] = 0.0;
            atoms[count].displacements[1] = 0.0;
            atoms[count].displacements[2] = 0.0;
            count++;
        }
    }

}

void move_atom(int* occupation, Atom* atoms, int num_atoms, int old_index, int new_index, double *displacement)
{
    int atom_index = 0;

    for (int i = 0; i < num_atoms; i++) {
        if (atoms[i].current_index == old_index) {atom_index = i; break;}
    }

    if (occupation[old_index] == 0) printf("Critical error! An unoccupaied site was used for move_atom\n");
    occupation[old_index] = 0;
    atoms[atom_index].current_index = new_index;
    occupation[new_index] = 1;

    atoms[atom_index].displacements[0] = atoms[atom_index].displacements[0] + displacement[0];
    atoms[atom_index].displacements[1] = atoms[atom_index].displacements[1] + displacement[1];
    atoms[atom_index].displacements[2] = atoms[atom_index].displacements[2] + displacement[2];

   //*msd = *msd + atoms[atom_index].displacements[0] * atoms[atom_index].displacements[0] + atoms[atom_index].displacements[1] * atoms[atom_index].displacements[1] + atoms[atom_index].displacements[2] * atoms[atom_index].displacements[2];

}

double calculate_msd(Atom* atoms, int num_atoms)
{
    double msd = 0;

    for (int i = 0; i < num_atoms; i ++) {
        msd = msd + atoms[i].displacements[0] * atoms[i].displacements[0] + atoms[i].displacements[1] * atoms[i].displacements[1] + atoms[i].displacements[2] * atoms[i].displacements[2];
    }

    msd = msd / num_atoms;

    return msd;
}
