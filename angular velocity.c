#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Configuration ---
// Directory containing snapshot files named like "time_%d.dat"
const char* DIRECTORY_PATH = "/home/somnath2/Codes/Trial_64/TA200_R2.5/evo/";
// Output file
const char* OUTPUT_FILE_PATH = "/home/somnath2/Codes/polar_order/Trial_96/mean_angular_velocity_200_2.5.dat";
// Which snapshot indices to consider
const int TIME_MIN = 0;
const int TIME_MAX = 2000;

// ---------------- Time scaling parameters (set these manually) ----------------
// Number of integration steps between successive output snapshots (N)
const long OUTPUT_EVERY_N = 100;          // <---- set this to the value used in your simulation
// Integration timestep used in the simulation
const double INTEGRATION_DT = 1e-5;       // integration time step (e.g. 1e-5)
// Rotational diffusion coefficient (Dr)
const double DR = 1.0;                    // set to 1 if Dr = 1

// If you want mean_omega_nondim computed as omega_real / DR, leave as is.
// If you prefer omega_nondim = omega_real * DR, change below accordingly.

typedef struct {
    int time_num;
    char path[1024];
} FileInfo;

int compare_files(const void* a, const void* b) {
    const FileInfo* fileA = (const FileInfo*)a;
    const FileInfo* fileB = (const FileInfo*)b;
    return (fileA->time_num - fileB->time_num);
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

// compute shortest angular difference wrapped to (-pi, pi]
static inline double angle_diff(double a, double b) {
    double d = a - b;
    d = fmod(d + M_PI, 2.0*M_PI);
    if (d < 0) d += 2.0*M_PI;
    d -= M_PI;
    return d;
}

int main() {
    printf("Scanning directory for data files...\n");

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
        fprintf(stderr, "Not enough files to compute angular velocity.\n");
        free(file_list);
        return 1;
    }

    qsort(file_list, file_count, sizeof(FileInfo), compare_files);
    printf("Found %d files.\n", file_count);

    // Determine number of particles from the first file (count only valid data lines)
    int num_particles = count_particles_in_file(file_list[0].path);
    if (num_particles <= 0) {
        fprintf(stderr, "Could not determine particle count (or file empty/non-data).\n");
        free(file_list);
        return 1;
    }
    printf("Detected %d particles per snapshot.\n", num_particles);
    printf("Time scaling: OUTPUT_EVERY_N = %ld, INTEGRATION_DT = %.6e, DR = %.6g\n", OUTPUT_EVERY_N, INTEGRATION_DT, DR);

    // Allocate memory for two snapshots
    double* prev_orientations = malloc(num_particles * sizeof(double));
    double* curr_orientations = malloc(num_particles * sizeof(double));
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
    // Header: time_index, time_real (s), time_nondim, mean_omega (rad/s), mean_omega_nondim, mean_abs_omega, rms_omega
    fprintf(out_file, "# time_index\ttime_real(s)\ttime_nondim\tmean_omega_rad_per_s\tmean_omega_nondim\tmean_abs_omega\t_rms_omega\n");

    // Read the first file (snapshot at file_list[0].time_num)
    if (read_orientations(file_list[0].path, prev_orientations, num_particles) != 0) {
        fprintf(stderr, "Error reading first file '%s'.\n", file_list[0].path);
        free(file_list);
        free(prev_orientations);
        free(curr_orientations);
        fclose(out_file);
        return 1;
    }

    // Loop over the rest of the files
    for (int t = 1; t < file_count; ++t) {
        int time_prev_idx = file_list[t-1].time_num;
        int time_curr_idx = file_list[t].time_num;
        long step_diff = (long)(time_curr_idx - time_prev_idx);
        if (step_diff <= 0) {
            fprintf(stderr, "Non-positive step difference between indices %d and %d. Skipping.\n", time_prev_idx, time_curr_idx);
            if (read_orientations(file_list[t].path, prev_orientations, num_particles) != 0) break;
            continue;
        }

        // Real time between snapshots in seconds (or simulation time units)
        double dt_real = (double)step_diff * (double)OUTPUT_EVERY_N * INTEGRATION_DT;
        // Nondimensional dt (scaled by Dr): t_nondim = t_real * Dr
        double dt_nondim = dt_real * DR;

        printf("Processing %s -> %s (step_diff=%ld, dt_real=%.6e, dt_nondim=%.6e)\n",
               file_list[t-1].path, file_list[t].path, step_diff, dt_real, dt_nondim);

        if (read_orientations(file_list[t].path, curr_orientations, num_particles) != 0) {
            fprintf(stderr, "Error reading file '%s'\n", file_list[t].path);
            break;
        }

        // Compute mean angular velocity and other stats
        double sum_omega = 0.0;
        double sum_abs_omega = 0.0;
        double sum_sq_omega = 0.0;
        for (int p = 0; p < num_particles; ++p) {
            double dtheta = angle_diff(curr_orientations[p], prev_orientations[p]);
            double omega_real = dtheta / dt_real;       // rad per unit real-time
            double omega_nondim = omega_real / DR;      // nondimensional (divide by Dr)
            sum_omega += omega_real;
            sum_abs_omega += fabs(omega_real);
            sum_sq_omega += omega_real * omega_real;
            // if you want to use nondimensional omega in stats, compute similarly
            (void)omega_nondim; // placeholder to show both are available
        }
        double mean_omega = sum_omega / (double)num_particles;
        double mean_abs_omega = sum_abs_omega / (double)num_particles;
        double rms_omega = sqrt(sum_sq_omega / (double)num_particles);
        double mean_omega_nondim = mean_omega / DR;

        // Write results using current time index and both real and nondimensional times
        fprintf(out_file, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                time_curr_idx, dt_real, dt_nondim, mean_omega, mean_omega_nondim, mean_abs_omega, rms_omega);

        // Swap buffers
        double* tmp = prev_orientations;
        prev_orientations = curr_orientations;
        curr_orientations = tmp;
    }

    printf("Processing complete. Results saved to %s\n", OUTPUT_FILE_PATH);

    // Cleanup
    fclose(out_file);
    free(prev_orientations);
    free(curr_orientations);
    free(file_list);

    return 0;
}
