#include "head.h"

int main(){

    int N;
    int Ntries;
    int MCsteps;
    int MCterm;
    int nbins, n_samples;
    int Nsim;

    double L;
    double *x, *y, *MSDx, *MSDy; 
    double *x0, *y0;
    double delta, sigma;
    double phi;
    double cpu_time_used;
    double gravity;
    double b;

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
    MCsteps = 20000; // Number of Monte Carlo steps
    MCterm = 30000; // Number of MC steps for equilibration
    gravity = 1; // Gravity strength
    // L = sqrt((PI*sigma*sigma*N)/(4*phi)); // Box size
    L = 39.633;
    phi = (PI*sigma*sigma*N)/(4*L*10*L);
    nbins = 75;
    b = (10*L)/nbins;
    n_samples = 0;
    Nsim = 10;
    double density_sims[nbins];
    for(int i = 0; i<nbins; i++){
        density_sims[i] = 0;
    }
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
    // sprintf(filename, "results/5_final_positions_grav_g%.2f_Lx%.1f_delta%.3f.txt", gravity, L, delta);
    // if ((fpos = fopen(filename, "w")) == NULL) {
    //     printf("Error opening file %s!\n", filename);
    //     return 1;
    // }
    // sprintf(filename, "results/5_energy_grav_g%.2f_Lx%.1f_delta%.3f.txt", gravity, L, delta);
    // if ((energy_file = fopen(filename, "w")) == NULL) {
    //     printf("Error opening file %s!\n", filename);
    //     return 1;
    // }
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

    // initialize_system(x, y, N, L, sigma, grav);
    // triangular_lattice(x, y, N, L);

    for(int i = 0; i<N; i++){
        x0[i] = x[i];
        y0[i] = y[i];
        MSDx[i] = x[i];
        MSDy[i] = y[i];
        // fprintf(initial_pos, "%lf\t%lf\n", x[i], y[i]);
    }

    // fprintf(energy_file, "0\t%lf\n", energy(y, gravity, N, 10*L));

    printf("Running Monte Carlo simulation with N=%d, phi=%.3f, sigma=%.3f, delta=%.5f, Ntries=%d, MCsteps=%d ...\n", N, phi, sigma, delta, Ntries, MCsteps);
    
    // for(int mc = 0; mc < MCsteps; mc++){
        
    //     for(int i = 0; i<Ntries; i++){
    //         metropolis(x, y, L, sigma, delta, N, MSDx, MSDy, hard_core, gravity, grav);
    //     }
    //     if (mc > (MCsteps-20000)){
    //         if (mc % 100 == 0){
    //             for(int i = 0; i<N; i++){
    //                 int bin = (int)(y[i]/b);
    //                 if (bin >= 0 && bin < nbins){
    //                     count[bin]++;
    //                 }
    //             }
    //             n_samples++;
    //         }
    //     }
    //     // if (mc == MCsteps-1){
    //     //     for(int i = 0; i<N; i++){
    //     //         fprintf(fpos, "%lf\t%lf\t", x[i], y[i]);
    //     //         fprintf(fpos, "\n");
    //     //     }
    //     // }
    //     // if (gravity != 0.0) {
    //     //     fprintf(energy_file, "%d\t%lf\n", mc, energy(y, gravity, N, 10*L));
    //     // }
    // }

    for(int sim = 0; sim<Nsim; sim++){

        initialize_system(x, y, N, L, sigma, grav);

        for(int equil = 0; equil<MCterm; equil++){
            for(int i = 0; i<Ntries; i++){
                metropolis(x, y, L, sigma, delta, N, MSDx, MSDy, hard_core, gravity, grav);
            }
        }

        int count[nbins];
        for(int i = 0; i<nbins; i++){
            count[i] = 0;
        }
        n_samples = 0;

        for(int mc = 0; mc < MCsteps; mc++){
            for(int i = 0; i<Ntries; i++){
                metropolis(x, y, L, sigma, delta, N, MSDx, MSDy, hard_core, gravity, grav);
            }

            for(int i = 0; i<N; i++){
                int bin = (int)((y[i]+10*L/2)/b);
                if (bin < 0) bin =0;
                if (bin >= nbins-1) bin = nbins-1;
                count[bin]++;
            }
            n_samples++;

        }

        double area = L * b;
        for(int i = 0; i<nbins; i++){
            double density = count[i] / ((double)n_samples * area);
            // printf("Sim %d, Bin %d, Count: %d, Density: %lf\n", sim, i, count[i], density);
            density_sims[i] += density;

            // printf("Simulation %d, Bin %d, Density: %lf\n", sim, i, density);
        }
    }

    
    FILE *density_file;
    sprintf(filename, "results/bonus_density_profile_grav_g%.2f_Lx%.1f_delta%.3f_nbin%d_T1.0.txt", gravity, L, delta, nbins);
    if ((density_file = fopen(filename, "w")) == NULL) {
        printf("Error opening file %s!\n", filename);
        return 1;
    }
    for(int i = 0; i<nbins; i++){
        double y_pos = -10*L/2.0 + (i + 0.5)*b;
        double density_sim = density_sims[i] / (double)Nsim;
        fprintf(density_file, "%lf\t%lf\n", y_pos, density_sim);
    }
    fclose(density_file);
        
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