#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#include "structure.h"
#include "poscar_reader.h"
#include "xdatcar_writer.h"
#include "kmc.h"


int main() {

    srand((unsigned int)time(NULL));

    int n_sites;
    double cutoff_radius = 4;
    double temperature = 1000;
    int kmc_steps = 10;
    int dump = 1;

    Topology structure;

    printf("Temperature (K): ");
    scanf("%lf", &temperature);

    printf("Cutoff radius (A): ");
    scanf("%lf", &cutoff_radius);

    printf("Number of KMC steps: ");
    scanf("%d", &kmc_steps); 
    
    printf("Dump: ");
    scanf("%d", &dump); 

    atomic_site *sites = read_poscar("POSCAR", &structure, &n_sites);

    // Neighbours
    int max_neighbors = MAX_NEIGHBORS;
    find_neighbors(sites, max_neighbors, cutoff_radius, &structure);
    write_topology("topology.txt", sites, &structure);  

    // Occupation
    int* occupation = (int*)malloc(n_sites * sizeof(int));
    initialize_occupation(occupation, &structure);
    int num_atoms = count_atoms(occupation, structure.num_sites);
    printf("Vacancy_concentration %.4lf\n", (double)num_atoms/structure.num_sites);

    Atom* atoms = malloc(num_atoms * sizeof(Atom));
    initialize_atoms(occupation, atoms, sites, n_sites);

    start_writing_xdatcar("XDATCAR", &structure, sites, atoms, num_atoms);
    int counter_xdatcar = 2;

    int num_barriers;
    barrier *barrier_list = build_barrier_list("barrier_list", &num_barriers);

    FILE *file = fopen("out.txt", "w");
    if (!file) {
        fprintf(stderr, "Error! %s.\n", "out.txt");
        return -1;
    }

    fprintf(file,"#Temperature: %.1lf K\n", temperature);
    fprintf(file,"#Cutoff radius: %.1lf A\n",cutoff_radius);
    fprintf(file,"#Number of KMC steps: %d\n",kmc_steps);
    fprintf(file,"#Dump: %d\n", dump);
    fprintf(file,"#Vacancy_concentration %.4lf\n", (double)num_atoms/structure.num_sites);

    fprintf(file, "#%-9s%-16s%-16s%-15s%-15s%-14s%-14s%-14s\n", "N", "From", "To", "Distance, A", "Barrier, eV", "Time, ps", "MSD, A**2", "Diffusion coefficient A**2/ps");

    double msd = 0;
    double timer = 0;
    double* msd_array = (double*)calloc(kmc_steps, sizeof(double));
    double* time_array = (double*)calloc(kmc_steps, sizeof(double));
    double* D_coef_array = (double*)calloc(kmc_steps + 1, sizeof(double));

    for (int i = 0; i < kmc_steps; i ++) {

        event current_event = select_event(occupation, sites, barrier_list, n_sites, num_barriers, temperature);
        move_atom(occupation, atoms, num_atoms, current_event.old_index, current_event.new_index, current_event.displacement);
        if ((i+1)%dump == 0) write_config("XDATCAR", &counter_xdatcar, sites, atoms, num_atoms);
        msd = calculate_msd(atoms, num_atoms);
        timer = timer + current_event.time_step;

        msd_array[i] = msd;
        time_array[i] = timer;
        D_coef_array[i+1] = msd_array[i] / time_array[i] / 6;

        fprintf(file, "%-10d%6s% -10d%6s% -10d%-15.4lf%-15.4lf%-14.4e%-14.4lf%-14.4e\n", i+1, current_event.fromA, current_event.old_index, 
        current_event.toB, current_event.new_index, current_event.dist, current_event.barrier, timer, msd, D_coef_array[i+1]);

    }

    double D_coef, intercept, r_value, std_D_coef;
    calculate_D_coef(time_array, msd_array, kmc_steps, &D_coef, &intercept, &r_value, &std_D_coef);

    fprintf(file, "#Diffusion coefficient = %.4e m**2/s\n", D_coef * 1E-8);
    fprintf(file, "#Standard deviation = %.4e m**2/s\n", std_D_coef * 1E-8);
    //fprintf(file, "#r-value = %.4lf\n", r_value);

    // End

    free(atoms);
    free(msd_array);
    free(time_array);
    free(occupation);

    return 0;
}
