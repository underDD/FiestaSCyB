#include "head.h"

int main(){

    int N;
    int Ntries;
    int MCsteps;

    double L;
    double *x, *y, *MSDx, *MSDy; 
    double *x0, *y0;
    double delta, sigma;
    double phi;
    double cpu_time_used;
    double gravity;

    char filename[MAX_STRING_LENGTH];

    FILE *f; 
    FILE *fpos;
    FILE *initial_pos;
    FILE *energy_file;

    clock_t start, end;

    bool hard_core, grav;

    // #####Initialization#####

    ini_ran(time(NULL)); // Seed the random number generator
    N = 1000; // Number of particles
    phi = 0.05; // Packing fraction
    sigma = 1; // Particle diameter
    delta = sigma/10; // Maximum displacement
    Ntries = 1000; // Number of tries
    MCsteps = 50000; // Number of Monte Carlo steps
    gravity = 10; // Gravity strength
    // L = sqrt((PI*sigma*sigma*N)/(4*phi)); // Box size
    L = 12.5;
    phi = (PI*sigma*sigma*N)/(4*L*10*L);
    printf("phi: %.5f\n", phi); 
    printf("Box size L: %lf\n", L);

    // sprintf(filename, "results/2_MSD_triang_phi%.3f_delta%.3f.txt", phi, delta);
    // if ((f = fopen(filename, "w")) == NULL) {
    //     printf("Error opening file %s!\n", filename);
    //     return 1;
    // }
    // sprintf(filename, "results/2_positions_triang_phi%.3f_delta%.3f.txt", phi, delta);
    // if ((fpos = fopen(filename, "w")) == NULL) {
    //     printf("Error opening file %s!\n", filename);
    //     return 1;
    // }
    sprintf(filename, "results/5_final_positions_grav_g%.2f_Lx%.1f_delta%.3f.txt", gravity, L, delta);
    if ((fpos = fopen(filename, "w")) == NULL) {
        printf("Error opening file %s!\n", filename);
        return 1;
    }
    sprintf(filename, "results/5_energy_grav_g%.2f_Lx%.1f_delta%.3f.txt", gravity, L, delta);
    if ((energy_file = fopen(filename, "w")) == NULL) {
        printf("Error opening file %s!\n", filename);
        return 1;
    }
    // sprintf(filename, "2_initial_positions_triang_phi0.5.txt");
    // initial_pos = fopen(filename, "w");
    // if (initial_pos == NULL) {
    //     printf("Error opening file %s!\n", filename);
    //     return 1;
    // }
    x = (double*) malloc(N * sizeof(double));
    y = (double*) malloc(N * sizeof(double));
    x0 = (double*) malloc(N * sizeof(double));
    y0 = (double*) malloc(N * sizeof(double));
    MSDx = (double*) malloc(N * sizeof(double));
    MSDy = (double*) malloc(N * sizeof(double));

    hard_core = true;
    grav = true;

    start = clock();
    // ########################

    initialize_system(x, y, N, L, sigma, grav);
    // triangular_lattice(x, y, N, L);

    for(int i = 0; i<N; i++){
        x0[i] = x[i];
        y0[i] = y[i];
        MSDx[i] = x[i];
        MSDy[i] = y[i];
        // fprintf(initial_pos, "%lf\t%lf\n", x[i], y[i]);
    }

    fprintf(energy_file, "0\t%lf\n", energy(y, gravity, N, 10*L));

    printf("Running Monte Carlo simulation with N=%d, phi=%.3f, sigma=%.3f, delta=%.5f, Ntries=%d, MCsteps=%d ...\n", N, phi, sigma, delta, Ntries, MCsteps);
    
    for(int mc = 0; mc < MCsteps; mc++){
        
        for(int i = 0; i<Ntries; i++){
            metropolis(x, y, L, sigma, delta, N, MSDx, MSDy, hard_core, gravity, grav);
        }
        if (mc == MCsteps-1){
            for(int i = 0; i<N; i++){
                fprintf(fpos, "%lf\t%lf\t", x[i], y[i]);
                fprintf(fpos, "\n");
            }
        }
        if (gravity != 0.0) {
            fprintf(energy_file, "%d\t%lf\n", mc, energy(y, gravity, N, 10*L));
        }
    }
    
        
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("CPU time used: %lf seconds\n", cpu_time_used);
    
    free(x);
    free(y);
    free(x0);
    free(y0);
    free(MSDx);
    free(MSDy);
    // fclose(f);
    fclose(fpos);
    fclose(energy_file);
    // fclose(initial_pos);
    return 0;
    
    }