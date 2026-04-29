#include "poscar_reader.h"
#include "structure.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double calculate_distance(atomic_site a, atomic_site b, Topology *crystall, double* temp_disp) {
    double dx, dy, dz;

    dx = a.direct_coordinates[0] - b.direct_coordinates[0];
    dy = a.direct_coordinates[1] - b.direct_coordinates[1];
    dz = a.direct_coordinates[2] - b.direct_coordinates[2];

    dx -= round(dx);
    dy -= round(dy);
    dz -= round(dz);

    double x = dx * crystall->lattice[0][0] + dy * crystall->lattice[1][0] + dz * crystall->lattice[2][0];
    double y = dx * crystall->lattice[0][1] + dy * crystall->lattice[1][1] + dz * crystall->lattice[2][1];
    double z = dx * crystall->lattice[0][2] + dy * crystall->lattice[1][2] + dz * crystall->lattice[2][2];

    temp_disp[0] = x;
    temp_disp[1] = y;
    temp_disp[2] = z;

    return sqrt(x * x + y * y + z * z);
}

double calculate_distance_old(atomic_site a, atomic_site b, Topology *crystall) {
    double dx, dy, dz;

    // Вычисляем разницу координат между атомами в декартовых координатах
    dx = a.cart_coordinates[0] - b.cart_coordinates[0];
    dy = a.cart_coordinates[1] - b.cart_coordinates[1];
    dz = a.cart_coordinates[2] - b.cart_coordinates[2];

    // Создаем вектор смещений в координатах базиса решетки
    double relative_disp[3];
    
    // Преобразование из декартовых координат в координаты базиса решетки
    // Для этого используем обратную матрицу решетки
    double determinant = crystall->lattice[0][0] * (crystall->lattice[1][1] * crystall->lattice[2][2] - crystall->lattice[1][2] * crystall->lattice[2][1])
                       - crystall->lattice[0][1] * (crystall->lattice[1][0] * crystall->lattice[2][2] - crystall->lattice[1][2] * crystall->lattice[2][0])
                       + crystall->lattice[0][2] * (crystall->lattice[1][0] * crystall->lattice[2][1] - crystall->lattice[1][1] * crystall->lattice[2][0]);

    double inv_lattice[3][3]; // Обратная матрица решетки
    inv_lattice[0][0] = (crystall->lattice[1][1] * crystall->lattice[2][2] - crystall->lattice[1][2] * crystall->lattice[2][1]) / determinant;
    inv_lattice[0][1] = (crystall->lattice[0][2] * crystall->lattice[2][1] - crystall->lattice[0][1] * crystall->lattice[2][2]) / determinant;
    inv_lattice[0][2] = (crystall->lattice[0][1] * crystall->lattice[1][2] - crystall->lattice[0][2] * crystall->lattice[1][1]) / determinant;
    inv_lattice[1][0] = (crystall->lattice[1][2] * crystall->lattice[2][0] - crystall->lattice[1][0] * crystall->lattice[2][2]) / determinant;
    inv_lattice[1][1] = (crystall->lattice[0][0] * crystall->lattice[2][2] - crystall->lattice[0][2] * crystall->lattice[2][0]) / determinant;
    inv_lattice[1][2] = (crystall->lattice[0][2] * crystall->lattice[1][0] - crystall->lattice[0][0] * crystall->lattice[1][2]) / determinant;
    inv_lattice[2][0] = (crystall->lattice[1][0] * crystall->lattice[2][1] - crystall->lattice[1][1] * crystall->lattice[2][0]) / determinant;
    inv_lattice[2][1] = (crystall->lattice[0][1] * crystall->lattice[2][0] - crystall->lattice[0][0] * crystall->lattice[2][1]) / determinant;
    inv_lattice[2][2] = (crystall->lattice[0][0] * crystall->lattice[1][1] - crystall->lattice[0][1] * crystall->lattice[1][0]) / determinant;

    // Преобразование в базис решетки
    relative_disp[0] = inv_lattice[0][0] * dx + inv_lattice[0][1] * dy + inv_lattice[0][2] * dz;
    relative_disp[1] = inv_lattice[1][0] * dx + inv_lattice[1][1] * dy + inv_lattice[1][2] * dz;
    relative_disp[2] = inv_lattice[2][0] * dx + inv_lattice[2][1] * dy + inv_lattice[2][2] * dz;

    // Применяем минимальное изображение с учетом периодических граничных условий
    relative_disp[0] -= round(relative_disp[0]);
    relative_disp[1] -= round(relative_disp[1]);
    relative_disp[2] -= round(relative_disp[2]);

    // Преобразуем относительное смещение обратно в декартовы координаты
    double x = relative_disp[0] * crystall->lattice[0][0] + relative_disp[1] * crystall->lattice[1][0] + relative_disp[2] * crystall->lattice[2][0];
    double y = relative_disp[0] * crystall->lattice[0][1] + relative_disp[1] * crystall->lattice[1][1] + relative_disp[2] * crystall->lattice[2][1];
    double z = relative_disp[0] * crystall->lattice[0][2] + relative_disp[1] * crystall->lattice[1][2] + relative_disp[2] * crystall->lattice[2][2];

    // Вычисляем евклидово расстояние
    return sqrt(x * x + y * y + z * z);
}

void convert_to_cartesian(atomic_site *sites, int num_sites, Topology *crystall) {
    for (int i = 0; i < num_sites; i++) {
        double x = sites[i].direct_coordinates[0] * crystall->lattice[0][0] + sites[i].direct_coordinates[1] * crystall->lattice[1][0] + sites[i].direct_coordinates[2] * crystall->lattice[2][0];
        double y = sites[i].direct_coordinates[0] * crystall->lattice[0][1] + sites[i].direct_coordinates[1] * crystall->lattice[1][1] + sites[i].direct_coordinates[2] * crystall->lattice[2][1];
        double z = sites[i].direct_coordinates[0] * crystall->lattice[0][2] + sites[i].direct_coordinates[1] * crystall->lattice[1][2] + sites[i].direct_coordinates[2] * crystall->lattice[2][2];

        sites[i].cart_coordinates[0] = x;
        sites[i].cart_coordinates[1] = y;
        sites[i].cart_coordinates[2] = z;
    }
}

void convert_to_direct(atomic_site *sites, int num_sites, Topology *crystall) {
    double determinant = crystall->lattice[0][0] * (crystall->lattice[1][1] * crystall->lattice[2][2] - crystall->lattice[1][2] * crystall->lattice[2][1])
                       - crystall->lattice[0][1] * (crystall->lattice[1][0] * crystall->lattice[2][2] - crystall->lattice[1][2] * crystall->lattice[2][0])
                       + crystall->lattice[0][2] * (crystall->lattice[1][0] * crystall->lattice[2][1] - crystall->lattice[1][1] * crystall->lattice[2][0]);

    double inv_lattice[3][3];
    inv_lattice[0][0] = (crystall->lattice[1][1] * crystall->lattice[2][2] - crystall->lattice[1][2] * crystall->lattice[2][1]) / determinant;
    inv_lattice[0][1] = (crystall->lattice[0][2] * crystall->lattice[2][1] - crystall->lattice[0][1] * crystall->lattice[2][2]) / determinant;
    inv_lattice[0][2] = (crystall->lattice[0][1] * crystall->lattice[1][2] - crystall->lattice[0][2] * crystall->lattice[1][1]) / determinant;
    inv_lattice[1][0] = (crystall->lattice[1][2] * crystall->lattice[2][0] - crystall->lattice[1][0] * crystall->lattice[2][2]) / determinant;
    inv_lattice[1][1] = (crystall->lattice[0][0] * crystall->lattice[2][2] - crystall->lattice[0][2] * crystall->lattice[2][0]) / determinant;
    inv_lattice[1][2] = (crystall->lattice[0][2] * crystall->lattice[1][0] - crystall->lattice[0][0] * crystall->lattice[1][2]) / determinant;
    inv_lattice[2][0] = (crystall->lattice[1][0] * crystall->lattice[2][1] - crystall->lattice[1][1] * crystall->lattice[2][0]) / determinant;
    inv_lattice[2][1] = (crystall->lattice[0][1] * crystall->lattice[2][0] - crystall->lattice[0][0] * crystall->lattice[2][1]) / determinant;
    inv_lattice[2][2] = (crystall->lattice[0][0] * crystall->lattice[1][1] - crystall->lattice[0][1] * crystall->lattice[1][0]) / determinant;

    for (int i = 0; i < num_sites; i++) {
        double x = sites[i].cart_coordinates[0];
        double y = sites[i].cart_coordinates[1];
        double z = sites[i].cart_coordinates[2];

        sites[i].direct_coordinates[0] = inv_lattice[0][0] * x + inv_lattice[0][1] * y + inv_lattice[0][2] * z;
        sites[i].direct_coordinates[1] = inv_lattice[1][0] * x + inv_lattice[1][1] * y + inv_lattice[1][2] * z;
        sites[i].direct_coordinates[2] = inv_lattice[2][0] * x + inv_lattice[2][1] * y + inv_lattice[2][2] * z;
    }
}

atomic_site* read_poscar(const char *filename, Topology *crystall, int* n_sites) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "%s doesn't exist.\n", filename);
        return NULL;
    }

    char line[MAX_LINE_LENGTH];

    // Comment
    fgets(line, MAX_LINE_LENGTH, file);

    // Scale factor
    fgets(line, MAX_LINE_LENGTH, file);
    if (sscanf(line, "%lf", &crystall->scale_factor) != 1) {
        fprintf(stderr, "Error! Wrong scale factor format.\n");
        fclose(file);
        return NULL;
    }

    // Lattice
    for (int i = 0; i < 3; i++) {
        fgets(line, MAX_LINE_LENGTH, file);
        sscanf(line, "%lf %lf %lf", &crystall->lattice[i][0], &crystall->lattice[i][1], &crystall->lattice[i][2]);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            crystall->lattice[i][j] *= crystall->scale_factor;
        }
    }

    // Atom types
    fgets(line, MAX_LINE_LENGTH, file);
    char site_types[MAX_TYPES][MAX_TYPE_LENGTH];
    int num_site_types = 0;
    char *token = strtok(line, " \t\n");
    while (token != NULL && num_site_types < MAX_TYPES) {
        strcpy(site_types[num_site_types], token);
        strcpy(crystall->types[num_site_types], token);
        num_site_types++;
        token = strtok(NULL, " \t\n");
    }

    crystall->num_types = num_site_types;

    // Atom counts
    fgets(line, MAX_LINE_LENGTH, file);
    int site_counts[MAX_TYPES];
    int total_sites = 0;
    token = strtok(line, " \t\n");
    for (int i = 0; i < num_site_types; i++) {
        site_counts[i] = atoi(token);
        total_sites += site_counts[i];
        crystall->count_types[i] = site_counts[i];
        token = strtok(NULL, " \t\n");
    }

    // Direct or Cartesian
    fgets(line, MAX_LINE_LENGTH, file);
    int is_direct = 0;
    if (strstr(line, "Direct") != NULL || strstr(line, "direct") != NULL) {
        is_direct = 1;
    }

    *n_sites = total_sites;
    crystall->num_sites = total_sites;


    atomic_site *sites = (atomic_site *)malloc(total_sites * sizeof(atomic_site));
    if (!sites) {
        fprintf(stderr, "Memory allocation failed.\n");
        fclose(file);
        return NULL;
    }

    int site_index = 0;

    // Coordinates
    if (is_direct) {
        for (int i = 0; i < num_site_types; i++) {
            for (int j = 0; j < site_counts[i]; j++) {
                if (fgets(line, MAX_LINE_LENGTH, file) != NULL) {
                    sscanf(line, "%lf %lf %lf", &sites[site_index].direct_coordinates[0], &sites[site_index].direct_coordinates[1], &sites[site_index].direct_coordinates[2]);
                    strcpy(sites[site_index].type, site_types[i]);
                    sites[site_index].index = site_index;
                    site_index++;
                }
            }
        }
    }
    else {
        for (int i = 0; i < num_site_types; i++) {
            for (int j = 0; j < site_counts[i]; j++) {
                if (fgets(line, MAX_LINE_LENGTH, file) != NULL) {
                    sscanf(line, "%lf %lf %lf", &sites[site_index].cart_coordinates[0], &sites[site_index].cart_coordinates[1], &sites[site_index].cart_coordinates[2]);
                    strcpy(sites[site_index].type, site_types[i]);
                    sites[site_index].index = site_index;
                    site_index++;
                }
            }
        }
    }

    fclose(file);

    if (is_direct) {
        convert_to_cartesian(sites, total_sites, crystall);
    }
    else convert_to_direct(sites, total_sites, crystall);

    return sites;
}

void find_neighbors(atomic_site *sites, int max_neighbors, double cutoff_radius, Topology *crystall) {
    for (int i = 0; i < crystall->num_sites; i++) {
        double distances[MAX_SITES];
        double neighbor_disp[MAX_SITES][3];
        int neighbor_indices[MAX_SITES];
        int num_neighbors = 0;

        for (int j = 0; j < crystall->num_sites; j++) {
            if (i != j) {
                double temp_disp[3] = {0};
                double dist = calculate_distance(sites[i], sites[j], crystall, temp_disp);
                if (dist > cutoff_radius) continue;
                distances[num_neighbors] = dist;
                neighbor_indices[num_neighbors] = j;
                neighbor_disp[num_neighbors][0] = temp_disp[0];
                neighbor_disp[num_neighbors][1] = temp_disp[1];
                neighbor_disp[num_neighbors][2] = temp_disp[2];
                num_neighbors++;
            }
        }

        for (int m = 0; m < num_neighbors - 1; m++) {
            for (int n = 0; n < num_neighbors - m - 1; n++) {
                if (distances[n] > distances[n + 1]) {
                    double temp_dist = distances[n];
                    distances[n] = distances[n + 1];
                    distances[n + 1] = temp_dist;

                    int temp_index = neighbor_indices[n];
                    neighbor_indices[n] = neighbor_indices[n + 1];
                    neighbor_indices[n + 1] = temp_index;

                    double temp_disp[3];
                    temp_disp[0] = neighbor_disp[n][0];
                    temp_disp[1] = neighbor_disp[n][1];
                    temp_disp[2] = neighbor_disp[n][2];

                    neighbor_disp[n][0] = neighbor_disp[n + 1][0];
                    neighbor_disp[n][1] = neighbor_disp[n + 1][1];
                    neighbor_disp[n][2] = neighbor_disp[n + 1][2];

                    neighbor_disp[n + 1][0] = temp_disp[0];
                    neighbor_disp[n + 1][1] = temp_disp[1];
                    neighbor_disp[n + 1][2] = temp_disp[2];
                }
            }
        }

        for (int k = 0; k < max_neighbors && k < num_neighbors; k++) {
            sites[i].neighbor_list[k] = sites[neighbor_indices[k]].index;
            sites[i].neighbor_dist[k] = distances[k];
            sites[i].neighbor_disp[k][0] = neighbor_disp[k][0];
            sites[i].neighbor_disp[k][1] = neighbor_disp[k][1];
            sites[i].neighbor_disp[k][2] = neighbor_disp[k][2];
        }

    }
}

void write_topology(const char *filename, atomic_site *atoms, Topology *crystall)
{
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error: %s.\n", filename);
        return;
    }

    fprintf(file, "General information:\n");
    fprintf(file, "Number of sites: %d\n", crystall->num_sites);
    fprintf(file, "Lattice vectors:\n");
    for (int i = 0; i < 3; i++) {
        fprintf(file, "%.6f %.6f %.6f\n", crystall->lattice[i][0], crystall->lattice[i][1], crystall->lattice[i][2]);
    }
    fprintf(file, "\n");

    for (int i = 0; i < crystall->num_sites; i++) {
        fprintf(file, "Atom %d (%s): Coordinates (%.6f, %.6f, %.6f)\n",
                atoms[i].index, atoms[i].type,
                atoms[i].cart_coordinates[0], atoms[i].cart_coordinates[1], atoms[i].cart_coordinates[2]);

        fprintf(file, "Nearest neighbours: ");
        fprintf(file, "\n");


        for (int j = 0; j < MAX_NEIGHBORS; j++) {
            fprintf(file, "%8d ", atoms[i].neighbor_list[j]);
        }

        fprintf(file, "\n");


        fprintf(file, "Distances ");
        fprintf(file, "\n");


        for (int j = 0; j < MAX_NEIGHBORS; j++) {
            fprintf(file, "%8.4lf ", atoms[i].neighbor_dist[j]);
        }

        fprintf(file, "\n\n");
        fprintf(file, "Displacements ");
        fprintf(file, "\n");


        for (int j = 0; j < MAX_NEIGHBORS; j++) {

                fprintf(file, "%8.4lf ", atoms[i].neighbor_disp[j][0]);
        }

        fprintf(file, "\n");
                for (int j = 0; j < MAX_NEIGHBORS; j++) {

                fprintf(file, "%8.4lf ", atoms[i].neighbor_disp[j][1]);
        }

        fprintf(file, "\n");
                for (int j = 0; j < MAX_NEIGHBORS; j++) {

                fprintf(file, "%8.4lf ", atoms[i].neighbor_disp[j][2]);
        }

        fprintf(file, "\n\n");
    }

    fclose(file);
    printf("File %s has been built.\n", filename);
}