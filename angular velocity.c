// mean_angular_velocity.c
// Compute mean angular velocity from evo/time_*.dat snapshots
// Configure the DIRECTORY_PATH, OUTPUT_FILE_PATH and time-scaling params below.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ----------------- CONFIGURATION -----------------
const char* DIRECTORY_PATH = "/home/somnath2/Codes/Trial_64/TA200_R2.5/evo/";
const char* OUTPUT_FILE_PATH = "/home/somnath2/Codes/polar_order/Trial_96/mean_angular_velocity_200_2.5.dat";
const int TIME_MIN = 0;
const int TIME_MAX = 2000;

// Simulation time scaling (set these to match your sim)
const long OUTPUT_EVERY_N = 100;     // number of integration steps between successive output snapshots
const double INTEGRATION_DT = 1e-5;  // integration timestep (s)
const double DR = 1.0;               // rotational diffusion coefficient (use 1.0 if Dr == 1)

// IMPORTANT: Set this depending on how your saved theta values are stored:
// 1 -> theta in files are UNWRAPPED (accumulated rotations). Use raw difference for omega.
// 0 -> theta in files are WRAPPED (in [-pi,pi) or [0,2pi)). Use minimal-angle difference.
const int DATA_UNWRAPPED = 1;

// -------------------------------------------------

typedef struct {
    int time_num;
    char path[1024];
} FileInfo;

int compare_files(const void* a, const void* b) {
    const FileInfo* A = (const FileInfo*)a;
    const FileInfo* B = (const FileInfo*)b;
    return (A->time_num - B->time_num);
}

// count only lines that contain three doubles
int count_particles_in_file(const char* filepath) {
    FILE* file = fopen(filepath, "r");
    if (!file) {
        perror("Error opening file for counting particles");
        return -1;
    }
    int count = 0;
    double d1, d2, d3;
    while (1) {
        int ret = fscanf(file, " %lf %lf %lf", &d1, &d2, &d3);
        if (ret == 3) {
            count++;
            continue;
        } else if (ret == EOF) {
            break;
        } else {
            int c;
            while ((c = fgetc(file)) != '\n' && c != EOF) {}
            continue;
        }
    }
    fclose(file);
    return count;
}

// read exactly num_particles third-column orientation values; skip bad lines
int read_orientations(const char* filepath, double* orientations, int num_particles) {
    FILE* file = fopen(filepath, "r");
    if (!file) {
        fprintf(stderr, "Error opening file '%s'\n", filepath);
        return -1;
    }
    int read = 0;
    double a,b,c;
    while (read < num_particles) {
        int ret = fscanf(file, " %lf %lf %lf", &a, &b, &c);
        if (ret == 3) {
            orientations[read++] = c;
        } else if (ret == EOF) {
            break;
        } else {
            int ch;
            while ((ch = fgetc(file)) != '\n' && ch != EOF) {}
        }
    }
    fclose(file);
    if (read != num_particles) {
        fprintf(stderr, "Warning: expected %d particles but read %d from '%s'\n", num_particles, read, filepath);
        return -1;
    }
    return 0;
}

// robust minimal-angle difference in (-pi, pi]
static inline double angle_diff_min(double a, double b) {
    double d = a - b;
    // atan2(sin(d), cos(d)) returns minimal signed difference in (-pi, pi]
    return atan2(sin(d), cos(d));
}

int main(void) {
    printf("Scanning directory: %s\n", DIRECTORY_PATH);

    DIR* d = opendir(DIRECTORY_PATH);
    if (!d) {
        fprintf(stderr, "Error opening directory '%s'\n", DIRECTORY_PATH);
        return 1;
    }

    struct dirent* dir;
    FileInfo* file_list = NULL;
    int file_count = 0;
    int file_capacity = 0;

    while ((dir = readdir(d)) != NULL) {
        int time_num;
        if (sscanf(dir->d_name, "time_%d.dat", &time_num) == 1) {
            if (time_num >= TIME_MIN && time_num <= TIME_MAX) {
                if (file_count >= file_capacity) {
                    file_capacity = (file_capacity == 0) ? 16 : file_capacity * 2;
                    FileInfo* new_list = realloc(file_list, file_capacity * sizeof(FileInfo));
                    if (!new_list) {
                        fprintf(stderr, "Memory allocation failed.\n");
                        free(file_list);
                        closedir(d);
                        return 1;
                    }
                    file_list = new_list;
                }
                file_list[file_count].time_num = time_num;
                snprintf(file_list[file_count].path, sizeof(file_list[file_count].path), "%s%s", DIRECTORY_PATH, dir->d_name);
                file_count++;
            }
        }
    }
    closedir(d);

    if (file_count < 2) {
        fprintf(stderr, "Not enough files to compute angular velocity (found %d files).\n", file_count);
        free(file_list);
        return 1;
    }

    qsort(file_list, file_count, sizeof(FileInfo), compare_files);
    printf("Found %d files between %d and %d.\n", file_count, TIME_MIN, TIME_MAX);

    // Determine number of particles from the first file (count only valid data lines)
    int num_particles = count_particles_in_file(file_list[0].path);
    if (num_particles <= 0) {
        fprintf(stderr, "Could not determine particle count (or file empty/non-data).\n");
        free(file_list);
        return 1;
    }
    printf("Detected %d particles per snapshot.\n", num_particles);
    printf("Time scaling: OUTPUT_EVERY_N = %ld, INTEGRATION_DT = %.6e, DR = %.6g\n",
           OUTPUT_EVERY_N, INTEGRATION_DT, DR);
    printf("DATA_UNWRAPPED = %d (1 = unwrapped snapshots; 0 = wrapped snapshots)\n", DATA_UNWRAPPED);

    // Allocate memory for two snapshots
    double* prev_orientations = malloc((size_t)num_particles * sizeof(double));
    double* curr_orientations = malloc((size_t)num_particles * sizeof(double));
    if (!prev_orientations || !curr_orientations) {
        fprintf(stderr, "Memory allocation failed for orientation arrays.\n");
        free(file_list);
        free(prev_orientations);
        free(curr_orientations);
        return 1;
    }

    FILE* out_file = fopen(OUTPUT_FILE_PATH, "w");
    if (!out_file) {
        fprintf(stderr, "Could not open output file '%s'.\n", OUTPUT_FILE_PATH);
        free(file_list);
        free(prev_orientations);
        free(curr_orientations);
        return 1;
    }

    // Header: time_index, cumulative_time_real(s), cumulative_time_nondim, mean_omega (rad/s), mean_omega_nondim, mean_abs_omega, rms_omega
    fprintf(out_file, "# time_index\tcumulative_time_real(s)\tcumulative_time_nondim\tmean_omega_rad_per_s\tmean_omega_nondim\tmean_abs_omega\trms_omega\n");

    // Read the first file (snapshot at file_list[0].time_num)
    if (read_orientations(file_list[0].path, prev_orientations, num_particles) != 0) {
        fprintf(stderr, "Error reading first file '%s'.\n", file_list[0].path);
        free(file_list);
        free(prev_orientations);
        free(curr_orientations);
        fclose(out_file);
        return 1;
    }

    // cumulative real time (seconds) since first snapshot
    double cumulative_real_time = 0.0;
    double cumulative_nondim_time = 0.0;

    for (int t = 1; t < file_count; ++t) {
        int idx_prev = file_list[t-1].time_num;
        int idx_curr = file_list[t].time_num;
        long idx_diff = (long)(idx_curr - idx_prev);
        if (idx_diff <= 0) {
            fprintf(stderr, "Non-positive step difference between indices %d and %d. Skipping.\n", idx_prev, idx_curr);
            if (read_orientations(file_list[t].path, prev_orientations, num_particles) != 0) break;
            continue;
        }

        // Real time between snapshots (seconds)
        double dt_real = (double)idx_diff * (double)OUTPUT_EVERY_N * INTEGRATION_DT;
        if (dt_real <= 0.0) {
            fprintf(stderr, "Warning: dt_real == 0 for indices %d -> %d. Skipping.\n", idx_prev, idx_curr);
            if (read_orientations(file_list[t].path, prev_orientations, num_particles) != 0) break;
            continue;
        }

        // nondimensional dt scaled by Dr
        double dt_nondim = dt_real * DR;

        // Update cumulative times (we report cumulative at current snapshot)
        cumulative_real_time += dt_real;
        cumulative_nondim_time += dt_nondim;

        printf("Processing %s -> %s (idx_diff=%ld, dt_real=%.6e s, cumulative=%.6e s)\n",
               file_list[t-1].path, file_list[t].path, idx_diff, dt_real, cumulative_real_time);

        if (read_orientations(file_list[t].path, curr_orientations, num_particles) != 0) {
            fprintf(stderr, "Error reading file '%s'.\n", file_list[t].path);
            break;
        }

        // Compute mean angular velocity and other stats
        double sum_omega = 0.0;
        double sum_abs_omega = 0.0;
        double sum_sq_omega = 0.0;

        for (int p = 0; p < num_particles; ++p) {
            double dtheta;
            if (DATA_UNWRAPPED) {
                // For unwrapped data, use raw difference (no modulo). This preserves large rotations.
                dtheta = curr_orientations[p] - prev_orientations[p];
                // optional check: if dtheta is extremely large relative to dt_real, warn (possible file mismatch)
                if (!isfinite(dtheta)) {
                    fprintf(stderr, "NaN/Inf encountered in dtheta at particle %d between files %s and %s\n",
                            p, file_list[t-1].path, file_list[t].path);
                    dtheta = 0.0;
                }
            } else {
                // For wrapped data, use minimal-angle difference
                dtheta = angle_diff_min(curr_orientations[p], prev_orientations[p]);
            }

            double omega_real = dtheta / dt_real;   // rad per real time unit (s)
            sum_omega += omega_real;
            sum_abs_omega += fabs(omega_real);
            sum_sq_omega += omega_real * omega_real;
        }

        double mean_omega = sum_omega / (double)num_particles;
        double mean_abs_omega = sum_abs_omega / (double)num_particles;
        double rms_omega = sqrt(sum_sq_omega / (double)num_particles);
        double mean_omega_nondim = mean_omega / DR;

        // Write out: time index (current snapshot), cumulative times, stats
        fprintf(out_file, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                file_list[t].time_num,
                cumulative_real_time,
                cumulative_nondim_time,
                mean_omega,
                mean_omega_nondim,
                mean_abs_omega,
                rms_omega);

        // swap buffers (fast)
        double* tmp = prev_orientations;
        prev_orientations = curr_orientations;
        curr_orientations = tmp;
    }

    printf("Processing complete. Results saved to %s\n", OUTPUT_FILE_PATH);

    fclose(out_file);
    free(prev_orientations);
    free(curr_orientations);
    free(file_list);

    return 0;
}
