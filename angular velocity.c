/*
 * mean_angular_velocity.c
 * -----------------------
 *
 * Title: Mean angular velocity (CLI-ready, GitHub-friendly)
 * Author: Somnath Roy
 * Created: 2025-11-19
 
 *
 * Purpose
 * -------
 * Post-process simulation snapshot files named `time_<INDEX>.dat` to compute
 * ensemble angular velocity statistics. The program reads a chosen column
 * (default: third column) from successive snapshot files, computes angular
 * displacements between snapshots (either using raw/unwrapped differences or
 * the minimal wrapped difference), and writes a tab-separated time series with
 * mean angular velocity, nondimensionalized mean, mean absolute value and RMS.
 *

 *
 * Input file requirements
 * -----------------------
 * - Files must be named `time_<index>.dat`, where <index> is an integer.
 * - Each data row should contain whitespace-separated numeric columns.
 * - By default the program reads column index 2 (0-based), i.e. the 3rd
 *   numeric column as the orientation (theta, in radians). Use -c to change.
 * - Non-numeric or incomplete lines are skipped.
 *
 * Output
 * ------
 * A tab-separated file with header:
 * # time_index\tcumulative_time_real(s)\tcumulative_time_nondim\tmean_omega_rad_per_s\tmean_omega_nondim\tmean_abs_omega\trms_omega
 *
 * Improvements in this version
 * ----------------------------
 * - Command-line options (see usage below) — no need to edit source to change
 *   paths or parameters.
 * - Configurable orientation column, wrapped/unwrapped choice, time range.
 * - Robust counting of numeric lines and configurable tolerance for missing
 *   lines (optionally fail on mismatch).
 * - Better error messages and exit codes for CI and automated runs.
 *
 * Usage examples
 * --------------
 * Compile:
 *   gcc -O2 -std=c11 -Wall -Wextra -o mean_angular_velocity mean_angular_velocity.c -lm
 *
 * Basic run (uses defaults inside program):
 *   ./mean_angular_velocity -d ./evo/ -o mean_omega.dat
 *
 * Provide integration scaling and treat snapshots as wrapped angles:
 *   ./mean_angular_velocity -d ./evo/ -o mean_omega.dat --output-every 10000 --dt 1e-5 --dr 1.0 --wrapped
 *
 * Process a subset of indices:
 *   ./mean_angular_velocity -d ./evo/ -o mean_omega.dat --time-min 0 --time-max 2000
 *
 * See --help for full option list.
 *
 * Notes for maintainers
 * ---------------------
 * - Consider adding unit tests in tests/ that contain small synthetic
 *   time_0.dat, time_1.dat files along with expected output.
 * - Add a Makefile and README.md describing repository structure.
 *
 * Changelog
 * ---------
 * 2025-11-19  v1.0  CLI-enabled, header, robust parsing and column selection
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <getopt.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Default configuration (overridden by CLI)
static const char *DEFAULT_DIRECTORY = "./";
static const char *DEFAULT_OUTPUT = "mean_angular_velocity.out";
static const int DEFAULT_TIME_MIN = INT_MIN; // accept any
static const int DEFAULT_TIME_MAX = INT_MAX;
static const long DEFAULT_OUTPUT_EVERY_N = 10000;
static const double DEFAULT_INTEGRATION_DT = 1e-5;
static const double DEFAULT_DR = 1.0;
static const int DEFAULT_DATA_UNWRAPPED = 1; // 1 = unwrapped (raw diffs); 0 = wrapped
static const int DEFAULT_ORIENTATION_COLUMN = 2; // 0-based index (third column by default)
static const int DEFAULT_FAIL_ON_MISMATCH = 0; // 1 = exit when particle counts mismatch

typedef struct {
    int time_num;
    char path[PATH_MAX];
} FileInfo;

int compare_files(const void* a, const void* b) {
    const FileInfo* A = (const FileInfo*)a;
    const FileInfo* B = (const FileInfo*)b;
    return (A->time_num - B->time_num);
}

// Read a single double from a token, returns 1 on success
static int read_double_from_token(const char* token, double* val) {
    char *endptr = NULL;
    errno = 0;
    double d = strtod(token, &endptr);
    if (endptr == token) return 0; // no conversion performed
    if (errno == ERANGE) return 0; // out of range
    // allow trailing whitespace only
    while (*endptr) {
        if (!isspace((unsigned char)*endptr)) return 0;
        endptr++;
    }
    *val = d;
    return 1;
}

// count numeric lines that have at least (orientation_col+1) numeric columns
int count_particles_in_file(const char* filepath, int orientation_col) {
    FILE* file = fopen(filepath, "r");
    if (!file) {
        perror("Error opening file for counting particles");
        return -1;
    }
    int count = 0;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, file)) != -1) {
        // tokenize and try parsing doubles
        int col = 0;
        char *saveptr = NULL;
        char *token = strtok_r(line, " \t\n\r", &saveptr);
        int numeric_columns = 0;
        while (token != NULL) {
            double tmp;
            if (read_double_from_token(token, &tmp)) {
                numeric_columns++;
            } else {
                // non-numeric token; mark as break
                numeric_columns = -1000000; // impossible
                break;
            }
            token = strtok_r(NULL, " \t\n\r", &saveptr);
        }
        if (numeric_columns >= (orientation_col + 1)) count++;
    }
    free(line);
    fclose(file);
    return count;
}

// read exactly num_particles values from orientation_col (0-based); skip bad lines
int read_orientations(const char* filepath, double* orientations, int num_particles, int orientation_col) {
    FILE* file = fopen(filepath, "r");
    if (!file) {
        fprintf(stderr, "Error opening file '%s': %s\n", filepath, strerror(errno));
        return -1;
    }
    int read_count = 0;
    char *line = NULL;
    size_t len = 0;
    ssize_t llen;
    while ((llen = getline(&line, &len, file)) != -1 && read_count < num_particles) {
        // parse tokens
        char *saveptr = NULL;
        char *token = strtok_r(line, " \t\n\r", &saveptr);
        int col = 0;
        int success = 1;
        double val = 0.0;
        while (token != NULL) {
            if (col == orientation_col) {
                if (read_double_from_token(token, &val)) {
                    orientations[read_count++] = val;
                } else {
                    success = 0;
                }
                break;
            }
            col++;
            token = strtok_r(NULL, " \t\n\r", &saveptr);
        }
        if (!token) {
            // line ended before reaching orientation_col -> skip
            continue;
        }
        if (!success) {
            // invalid token -> skip
            continue;
        }
    }
    free(line);
    fclose(file);
    if (read_count != num_particles) {
        fprintf(stderr, "Warning: expected %d particles but read %d from '%s'\n", num_particles, read_count, filepath);
        return -1;
    }
    return 0;
}

// minimal-angle difference in (-pi, pi]
static inline double angle_diff_min(double a, double b) {
    double d = a - b;
    return atan2(sin(d), cos(d));
}

void print_usage(const char *prog) {
    fprintf(stderr,
            "Usage: %s [options]\n"
            "Options:\n"
            "  -d, --dir DIR             directory containing time_<N>.dat files (default ./)\n"
            "  -o, --output FILE         output file path (default mean_angular_velocity.out)\n"
            "  -m, --time-min INT        minimum index to consider (inclusive)\n"
            "  -M, --time-max INT        maximum index to consider (inclusive)\n"
            "  -n, --output-every N      saved frames between outputs (default 10000)\n"
            "  -t, --dt DOUBLE          integration timestep in seconds (default 1e-5)\n"
            "  -r, --dr DOUBLE          rotational diffusion coefficient Dr (default 1.0)\n"
            "  -u, --unwrapped          treat stored angles as unwrapped (default)\n"
            "  -w, --wrapped            treat stored angles as wrapped in [-pi,pi)\n"
            "  -c, --col INT            orientation column index (0-based, default 2)\n"
            "  -f, --fail-on-mismatch   exit on particle-count mismatch (default: warn only)\n"
            "  -h, --help               show this help and exit\n",
            prog);
}

int main(int argc, char *argv[]) {
    const char *directory = DEFAULT_DIRECTORY;
    const char *output_path = DEFAULT_OUTPUT;
    int time_min = DEFAULT_TIME_MIN;
    int time_max = DEFAULT_TIME_MAX;
    long output_every_n = DEFAULT_OUTPUT_EVERY_N;
    double integration_dt = DEFAULT_INTEGRATION_DT;
    double dr = DEFAULT_DR;
    int data_unwrapped = DEFAULT_DATA_UNWRAPPED;
    int orientation_col = DEFAULT_ORIENTATION_COLUMN;
    int fail_on_mismatch = DEFAULT_FAIL_ON_MISMATCH;

    static struct option long_options[] = {
        {"dir", required_argument, 0, 'd'},
        {"output", required_argument, 0, 'o'},
        {"time-min", required_argument, 0, 'm'},
        {"time-max", required_argument, 0, 'M'},
        {"output-every", required_argument, 0, 'n'},
        {"dt", required_argument, 0, 't'},
        {"dr", required_argument, 0, 'r'},
        {"unwrapped", no_argument, 0, 'u'},
        {"wrapped", no_argument, 0, 'w'},
        {"col", required_argument, 0, 'c'},
        {"fail-on-mismatch", no_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "d:o:m:M:n:t:r:uwc:fh", long_options, NULL)) != -1) {
        switch (opt) {
            case 'd': directory = optarg; break;
            case 'o': output_path = optarg; break;
            case 'm': time_min = atoi(optarg); break;
            case 'M': time_max = atoi(optarg); break;
            case 'n': output_every_n = strtol(optarg, NULL, 10); break;
            case 't': integration_dt = strtod(optarg, NULL); break;
            case 'r': dr = strtod(optarg, NULL); break;
            case 'u': data_unwrapped = 1; break;
            case 'w': data_unwrapped = 0; break;
            case 'c': orientation_col = atoi(optarg); break;
            case 'f': fail_on_mismatch = 1; break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    printf("Scanning directory: %s\n", directory);

    DIR* d = opendir(directory);
    if (!d) {
        fprintf(stderr, "Error opening directory '%s': %s\n", directory, strerror(errno));
        return 1;
    }

    struct dirent* direntp;
    FileInfo* file_list = NULL;
    int file_count = 0;
    int file_capacity = 0;

    while ((direntp = readdir(d)) != NULL) {
        int time_num;
        if (sscanf(direntp->d_name, "time_%d.dat", &time_num) == 1) {
            if (time_num >= time_min && time_num <= time_max) {
                if (file_count >= file_capacity) {
                    file_capacity = (file_capacity == 0) ? 32 : file_capacity * 2;
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
                snprintf(file_list[file_count].path, sizeof(file_list[file_count].path), "%s/%s", directory, direntp->d_name);
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
    printf("Found %d files between indices (considering time_min/time_max).\n", file_count);

    // Determine number of particles from the first file
    int num_particles = count_particles_in_file(file_list[0].path, orientation_col);
    if (num_particles <= 0) {
        fprintf(stderr, "Could not determine particle count (or file empty/non-data).\n");
        free(file_list);
        return 1;
    }
    printf("Detected %d particles per snapshot (orientation column %d).\n", num_particles, orientation_col);
    printf("Time scaling: OUTPUT_EVERY_N = %ld, INTEGRATION_DT = %.6e, DR = %.6g\n",
           output_every_n, integration_dt, dr);
    printf("DATA_UNWRAPPED = %d (1 = unwrapped snapshots; 0 = wrapped snapshots)\n", data_unwrapped);

    double* prev_orientations = malloc((size_t)num_particles * sizeof(double));
    double* curr_orientations = malloc((size_t)num_particles * sizeof(double));
    if (!prev_orientations || !curr_orientations) {
        fprintf(stderr, "Memory allocation failed for orientation arrays.\n");
        free(file_list);
        free(prev_orientations);
        free(curr_orientations);
        return 1;
    }

    FILE* out_file = fopen(output_path, "w");
    if (!out_file) {
        fprintf(stderr, "Could not open output file '%s' for writing: %s\n", output_path, strerror(errno));
        free(file_list);
        free(prev_orientations);
        free(curr_orientations);
        return 1;
    }

    fprintf(out_file, "# time_index\tcumulative_time_real(s)\tcumulative_time_nondim\tmean_omega_rad_per_s\tmean_omega_nondim\tmean_abs_omega\trms_omega\n");

    if (read_orientations(file_list[0].path, prev_orientations, num_particles, orientation_col) != 0) {
        fprintf(stderr, "Error reading first file '%s'.\n", file_list[0].path);
        free(file_list);
        free(prev_orientations);
        free(curr_orientations);
        fclose(out_file);
        return 1;
    }

    double cumulative_real_time = 0.0;
    double cumulative_nondim_time = 0.0;

    for (int t = 1; t < file_count; ++t) {
        int idx_prev = file_list[t-1].time_num;
        int idx_curr = file_list[t].time_num;
        long idx_diff = (long)(idx_curr - idx_prev);
        if (idx_diff <= 0) {
            fprintf(stderr, "Non-positive step difference between indices %d and %d. Skipping.\n", idx_prev, idx_curr);
            if (read_orientations(file_list[t].path, prev_orientations, num_particles, orientation_col) != 0) break;
            continue;
        }

        double dt_real = (double)idx_diff * (double)output_every_n * integration_dt;
        if (dt_real <= 0.0) {
            fprintf(stderr, "Warning: dt_real == 0 for indices %d -> %d. Skipping.\n", idx_prev, idx_curr);
            if (read_orientations(file_list[t].path, prev_orientations, num_particles, orientation_col) != 0) break;
            continue;
        }

        double dt_nondim = dt_real * dr;
        cumulative_real_time += dt_real;
        cumulative_nondim_time += dt_nondim;

        printf("Processing %s -> %s (idx_diff=%ld, dt_real=%.6e s, cumulative=%.6e s)\n",
               file_list[t-1].path, file_list[t].path, idx_diff, dt_real, cumulative_real_time);

        if (read_orientations(file_list[t].path, curr_orientations, num_particles, orientation_col) != 0) {
            fprintf(stderr, "Error reading file '%s'.\n", file_list[t].path);
            if (fail_on_mismatch) break;
            // attempt to continue by copying prev into curr (no motion)
            for (int p = 0; p < num_particles; ++p) curr_orientations[p] = prev_orientations[p];
        }

        double sum_omega = 0.0;
        double sum_abs_omega = 0.0;
        double sum_sq_omega = 0.0;

        for (int p = 0; p < num_particles; ++p) {
            double dtheta;
            if (data_unwrapped) {
                dtheta = curr_orientations[p] - prev_orientations[p];
                if (!isfinite(dtheta)) {
                    fprintf(stderr, "NaN/Inf encountered in dtheta at particle %d between files %s and %s\n",
                            p, file_list[t-1].path, file_list[t].path);
                    dtheta = 0.0;
                }
            } else {
                dtheta = angle_diff_min(curr_orientations[p], prev_orientations[p]);
            }

            double omega_real = dtheta / dt_real;
            sum_omega += omega_real;
            sum_abs_omega += fabs(omega_real);
            sum_sq_omega += omega_real * omega_real;
        }

        double mean_omega = sum_omega / (double)num_particles;
        double mean_abs_omega = sum_abs_omega / (double)num_particles;
        double rms_omega = sqrt(sum_sq_omega / (double)num_particles);
        double mean_omega_nondim = mean_omega / dr;

        fprintf(out_file, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                file_list[t].time_num,
                cumulative_real_time,
                cumulative_nondim_time,
                mean_omega,
                mean_omega_nondim,
                mean_abs_omega,
                rms_omega);

        double* tmp = prev_orientations;
        prev_orientations = curr_orientations;
        curr_orientations = tmp;
    }

    printf("Processing complete. Results saved to %s\n", output_path);

    fclose(out_file);
    free(prev_orientations);
    free(curr_orientations);
    free(file_list);

    return 0;
}
