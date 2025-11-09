#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
#include <stdbool.h>

#define MAX_STRING_LENGTH 100000
#define fran rand()/((double)RAND_MAX+1)
#define PI 2*asin(1)
#define EPS 1E-5
#define NMAX 10000

//**********PARISI RAPUANO*************
#define NormRANu (2.3283063671E-10F)

extern unsigned int irr[256];
extern unsigned int ir1;
extern unsigned char ind_ran,ig1,ig2,ig3;

//************************************

//**********functions.c*************  

float randomIn(double min, double max);
float Random(void);
void ini_ran(int SEMILLA);
float randomInPR(double min, double max);
double box_muller(void);
double sqrt_distance(double x, double y);
void initialize_system(double *x, double *y, int N, double L, double sigma, bool grav);
bool touches(int i, double *x, double *y, double sigma, double L, int N, bool grav);
void metropolis(double *x, double *y, double L, double sigma, double delta, int N, double *MSDx, double *MSDy, bool hard_core, double gravity, bool grav);
double MSD(double *x0, double *y0, double *x, double *y, int N);
void triangular_lattice(double *x, double *y, int N, double L);
void print_positions_time(double *x, double *y, int N, FILE *f, double time);
bool touches_wall(int i, double *x, double *y, double sigma, double Lx, double Ly);
void PBC(double *r, double L);
double mean(double *data, int size);
double energy(double *y, double gravity, int N, double Ly);
//************************************