#ifndef POSCAR_READER_H
#define POSCAR_READER_H

#define MAX_SITES 1000
#define MAX_LINE_LENGTH 256
#define MAX_NEIGHBORS 20
#define MAX_TYPE_LENGTH 10
#define MAX_TYPES 10

typedef struct {
    double cart_coordinates[3];
    double direct_coordinates[3];
    char type[MAX_TYPE_LENGTH];
    int index;
    int neighbor_list[MAX_NEIGHBORS];
    double neighbor_dist[MAX_NEIGHBORS];
    double neighbor_disp[MAX_NEIGHBORS][3];
} atomic_site;

typedef struct {
    double scale_factor;
    double lattice[3][3];
    int num_types;
    int count_types[MAX_TYPES];
    char types[MAX_TYPES][MAX_TYPE_LENGTH];
    int num_sites;
    int num_atoms;
    double vacancy_concentration[MAX_TYPES];
} Topology;

atomic_site* read_poscar(const char *filename, Topology *crystall, int* n_sites);
void convert_to_cartesian(atomic_site *sites, int num_sites, Topology *crystall);
double calculate_distance(atomic_site a, atomic_site b, Topology *crystall, double* temp_disp);
void find_neighbors(atomic_site *sites, int max_neighbors, double cutoff_radius, Topology *crystall);
void write_topology(const char *filename, atomic_site *atoms, Topology *crystall);

#endif // POSCAR_READER_H
