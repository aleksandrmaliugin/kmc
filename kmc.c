#include "poscar_reader.h"
#include "structure.h"
#include "xdatcar_writer.h"
#include "kmc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define MAX_NEIGHBORS 20
#define MAX_EVENTS 45000
#define MAX_BARRIERS 30
#define k_b 0.000086173303

event select_event(int *occupation, atomic_site *sites, barrier *barrier_list, int num_sites, int num_barriers, double temperature)
{
    event selected_event;
    event event_list[MAX_EVENTS];
    int count_events = 0;

    for (int i = 0; i < num_sites; i++){

        if (occupation[i] == 0) {
            
            for (int j = 0; j < MAX_NEIGHBORS; j++){
                if (sites[i].neighbor_dist[j] != 0.0 && occupation[sites[i].neighbor_list[j]] != 0) {
                    event_list[count_events].old_index = sites[i].neighbor_list[j];
                    event_list[count_events].new_index = i;
                    event_list[count_events].dist = sites[i].neighbor_dist[j];
                    int barrier_num = get_barrier(barrier_list, sites[i].neighbor_dist[j], 
                                                                   num_barriers, sites[sites[i].neighbor_list[j]].type,
                                                                   sites[i].type);
                    event_list[count_events].barrier = barrier_list[barrier_num].barrier_ev;
                    event_list[count_events].frequency = barrier_list[barrier_num].frequency;
                    event_list[count_events].rate = event_list[count_events].frequency * exp(-event_list[count_events].barrier / (k_b * temperature));
                    event_list[count_events].displacement[0] = sites[i].neighbor_disp[j][0];
                    event_list[count_events].displacement[1] = sites[i].neighbor_disp[j][1];
                    event_list[count_events].displacement[2] = sites[i].neighbor_disp[j][2];
                    strcpy(event_list[count_events].fromA, sites[sites[i].neighbor_list[j]].type);
                    strcpy(event_list[count_events].toB, sites[i].type);
                    count_events++;
                }
            }

        }

    }

    double *cumulative_function = (double *)calloc(count_events+1, sizeof(double));

    for (int i = 1; i < count_events + 1; i++) {   // First element = 0
        cumulative_function[i] = cumulative_function[i - 1] + event_list[i-1].rate;
    }

    selected_event.time_step = 1 / cumulative_function[count_events]; //redefine


    //for (int i = 1; i < count_events + 1; i++) {
    //    cumulative_function[i] = cumulative_function[i] / cumulative_function[count_events];
    //}

    //unsigned int seed = (unsigned int)time(NULL);
    //double randomValue = (double)rand_r(&seed) / (RAND_MAX + 1.0);
    double randomValue = (double)rand() / (RAND_MAX + 1.0);


    for (int i = 1; i < count_events + 1; i++) {
        if (cumulative_function[i] > randomValue*cumulative_function[count_events] && cumulative_function[i-1] < randomValue*cumulative_function[count_events]) {
            
            selected_event.old_index = event_list[i-1].old_index;
            selected_event.new_index = event_list[i-1].new_index;
            selected_event.dist = event_list[i-1].dist;
            selected_event.barrier = event_list[i-1].barrier;
            selected_event.frequency = event_list[i-1].frequency;
            selected_event.rate = event_list[i-1].rate;
            selected_event.displacement[0] = event_list[i-1].displacement[0];
            selected_event.displacement[1] = event_list[i-1].displacement[1];
            selected_event.displacement[2] = event_list[i-1].displacement[2];
            strcpy(selected_event.fromA, event_list[i - 1].fromA);
            strcpy(selected_event.toB, event_list[i - 1].toB);
            break;
        }
    }

    free(cumulative_function);

    return selected_event;

}

int get_barrier(barrier *barrier_list, double dist, int num_barriers, const char *fromA, const char *toB) {

    for (int i = 0; i < num_barriers; i++) {
        if (fabs(barrier_list[i].distance - dist) < 0.05 &&
            strcmp(barrier_list[i].fromA, fromA) == 0 &&
            strcmp(barrier_list[i].toB, toB) == 0) {
            return i;
        }
    }


    printf("Critical error! Check barrier list for types: %s -> %s! (%lf)\n ", fromA, toB, dist);

    return -1;
}


barrier* build_barrier_list(const char *filename, int *num_barriers) {
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s.\n", filename);
        return NULL;
    }

    barrier *barrier_list = (barrier *)calloc(MAX_BARRIERS,  sizeof(barrier));
    if (barrier_list == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        fclose(file);
        return NULL;
    }

    printf("Building barrier list\n");

    int counter = 0;
    char line[MAX_LINE_LENGTH];
    fgets(line, MAX_LINE_LENGTH, file); // Skip comment

    while (fgets(line, MAX_LINE_LENGTH, file) != NULL) {
        if (counter >= MAX_BARRIERS) {
            fprintf(stderr, "Warning: Maximum number of barriers reached.\n");
            break;
        }
        if (sscanf(line, "%s %s %lf %lf %lf", 
                   barrier_list[counter].fromA,
                   barrier_list[counter].toB,
                   &barrier_list[counter].distance, 
                   &barrier_list[counter].barrier_ev,
                   &barrier_list[counter].frequency
                   ) == 5) {
            printf("%s %s %lf %lf %lf\n", barrier_list[counter].fromA, barrier_list[counter].toB, barrier_list[counter].distance, barrier_list[counter].barrier_ev, barrier_list[counter].frequency);
            counter++;
        }
    }

    fclose(file);
    *num_barriers = counter;
    printf("Total barriers: %d\n", counter);
    return barrier_list;
}

void calculate_D_coef(const double x[], const double y[], int n, double *slope, double *intercept, double *r_value, double *stderr_val)
{
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_yy = 0.0, sum_xy = 0.0;

    for (int i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xx += x[i] * x[i];
        sum_yy += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }

    double mean_x = sum_x / n;
    double mean_y = sum_y / n;

    *slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x) / 6; // 3D-case
    *intercept = mean_y - (*slope) * mean_x;

    double r_num = n * sum_xy - sum_x * sum_y;
    double r_den = sqrt((n * sum_xx - sum_x * sum_x) * (n * sum_yy - sum_y * sum_y));
    *r_value = r_num / r_den;

    double ss_res = 0.0;
    for (int i = 0; i < n; i++) {
        double y_pred = (*slope) * x[i] + (*intercept);
        ss_res += (y[i] - y_pred) * (y[i] - y_pred);
    }
    *stderr_val = sqrt(ss_res / (n - 2)) / sqrt(sum_xx - (sum_x * sum_x) / n);
}