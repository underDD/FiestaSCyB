#include "head.h"
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;

/*
    Generate a random number between min and max
*/
float randomIn(double min, double max) {
    return min + (max - min) * fran;
}

/*
    Initialize the random number generator
    SEMILLA: seed
*/
void  ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;

    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }

    ind_ran=ig1=ig2=ig3=0;
}

/*
    Generate a random number between 0 and 1 using the Parisi-Rapuano method
*/
float Random(void)
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;

    return r;
}

/*
    Generate a random number between min and max
*/
float randomInPR(double min, double max) {
    return min + (max - min) * Random();
}


/*
    Generate a random number using the Box-Muller method
*/
double box_muller(){

    double r1, r2, s, u1;
    
    // Generate two uniform random numbers between 0 and 1
    r1 = Random();
    r2 = Random();

    // Apply the Box-Muller transform
    s = sqrt(-2.0 * log(r1));
    u1 = s * cos(2.0 * PI * r2);
    //*u2 = s * sin(2.0 * M_PI * r2);
    return u1;

}

double sqrt_distance(double x, double y){
    return sqrt(x*x + y*y);
}

void PBC(double *x, double L){
    if (L <= 0.0) return;
    double r = fmod(*x + L/2.0, L);   // desplazar para obtener resto en torno a 0..L
    if (r < 0.0) r += L;              // fmod puede devolver negativo
    *x = r - L/2.0;                   // volver al intervalo [-L/2, L/2)
}

bool touches(int i, double *x, double *y, double sigma, double L, int N, bool grav){

    double dx,dy;
    double Lx = L;
    double Ly = 10*Lx;

    if (grav){
        for(int j = 0; j<N; j++){
            if (j == i) continue;
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            if (fabs(dx) > Lx/2) dx = dx - copysign(Lx, dx);
            if (fabs(dy) > Ly/2) dy = dy - copysign(Ly, dy);
            if(sqrt_distance(dx, dy) < sigma){
                return true;
            }
        }
    }
    else{
        for(int j = 0; j<N; j++){
            if(j == i) continue;
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            if (fabs(dx) > L/2) dx = dx - copysign(L, dx);
            if (fabs(dy) > L/2) dy = dy - copysign(L, dy);
            if(sqrt_distance(dx, dy) < sigma){
                return true;
            }
        }
    }
    
    return false;

}

bool touches_wall(int i, double *x, double *y, double sigma, double Lx, double Ly){

    double radius = sigma/2.0;
    if(x[i] < -Lx/2 + radius || x[i] > Lx/2 - radius) return true;
    if(y[i] < -Ly/2 + radius || y[i] > Ly/2 - radius) return true;
    return false;
}

void initialize_system(double *x, double *y, int N, double L, double sigma, bool grav){

    double Lx = L;
    double Ly = 10*Lx;

    if (grav){
        
        for(int i = 0; i<N; i++){
            do{
                x[i] = randomIn(-Lx/2, Lx/2);
                y[i] = randomIn(-Ly/2, Ly/2);
            }while(touches(i, x, y, sigma, Lx, N, grav) || touches_wall(i, x, y, sigma, Lx, Ly));
            // printf("Initialized particle %d at (%lf, %lf)\n", i, x[i], y[i]);
        }
    }
    else{
        for(int i = 0; i<N; i++){
            do{
                x[i] = randomIn(-L/2, L/2);
                y[i] = randomIn(-L/2, L/2);
            }while(touches(i, x, y, sigma, L, N, grav));
        }
    }

}

void metropolis(double *x, double *y, double L, double sigma, double delta, int N, double *MSDx, double *MSDy, bool hard_core, double gravity, bool grav){
    
    int n;
    double dx, dy;
    double old_x, old_y, old_MSDx, old_MSDy;
    double Lx = L;
    double Ly = 10*Lx;

    double kBT = 1;

    n = (int)randomIn(0, N);
    dx = delta*randomIn(-0.5, 0.5);
    dy = delta*randomIn(-0.5, 0.5);
    old_x = x[n];
    old_y = y[n];
    old_MSDx = MSDx[n];
    old_MSDy = MSDy[n];
    x[n] = x[n] + dx;
    y[n] = y[n] + dy;

    
    if(hard_core && grav){
        double old_energy = gravity * old_y;
        double new_energy = gravity * y[n];
        double delta_E = new_energy - old_energy;
        
        if(touches(n, x, y, sigma, L, N, grav) || touches_wall(n, x, y, sigma, Lx, Ly)){
            x[n] = old_x;
            y[n] = old_y;
            return;
        }
        
        if(delta_E > 0.0){
            double prob = exp(-delta_E/kBT);
            if(randomIn(0.0, 1.0) > prob){
                x[n] = old_x;
                y[n] = old_y;
                return;
            }
        }
        
    }else{

        PBC(&x[n], L);
        PBC(&y[n], L);

        if(touches(n, x, y, sigma, L, N, grav)){
            x[n] = old_x;
            y[n] = old_y;
            MSDx[n] = old_MSDx;
            MSDy[n] = old_MSDy;
            return;
        }
        else{
            MSDx[n] = old_MSDx + dx;
            MSDy[n] = old_MSDy + dy;
        }
    }

}

double MSD(double *x0, double *y0, double *x, double *y, int N){
    double msd = 0.0;
    double dx, dy;

    for(int i = 0; i<N; i++){
        dx = x[i] - x0[i];
        dy = y[i] - y0[i];
        msd += dx*dx + dy*dy;
    }

    return msd/N;
}

void triangular_lattice(double *x, double *y, int N, double L){
    
    double a = L*sqrt(2.0/sqrt(3.0))/sqrt(N); // Lattice constant
    double ax = a;
    double ay = a*sqrt(3.0)/2.0;
    int ny = (int)(L/ay)+1;
    int nx = (int)(L/ax)+1;

    int part = 0; 
    
    for(int i = 0; i<ny && part<N ; i++){
        double yi = -L/2.0+ay*i;
        for(int j = 0; j<nx && part<N ; j++){
            double xi = -L/2.0 + ax*j + (i%2)*ax/2.0;
            x[part] = xi;
            y[part] = yi;
            part++;
        }
    }

}

void print_positions_time(double *x, double *y, int N, FILE *f, double time){

    fprintf(f, "%lf\t", time);
    for(int i = 0; i<N; i++){
        fprintf(f, "%lf\t%lf\t", x[i], y[i]);
    }
    fprintf(f, "\n");
}

double mean(double *data, int size){
    double sum = 0.0;
    for(int i = 0; i<size; i++){
        sum += data[i];
    }
    return sum/size;
}

double energy(double *y, double gravity, int N, double Ly){
    double E = 0.0;
    for(int i = 0; i<N; i++){
        E += gravity * (y[i] + Ly/2.0);
    }
    return E;
}