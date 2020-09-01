#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include <time.h>

const int N = 5000;             // Iterations
const int frame_interval = 500; // Iterations between exporting data
const int NX = 520;             // x resolution
const int NY = 180;             // y resolution

const double uw = 0.1;          // Initial fluid speed in x
const double Re = 55.;          // Reynolds number
const double rho_0 = 1.;        // Initial pressure

const char file_path[128] = "omp_img/"; // Path where files are saved

const int ex[] = {0, 1, -1, 0, 0, 1, -1, -1, 1};    // Vector directions of distribution func.
const int ey[] = {0, 0, 0, 1, -1, 1, -1, 1, -1};    // Vector directions of distribution func.

const double w[] = {4. / 9., 1. / 9., 1. / 9.,      // Weights of vectors
                    1. / 9., 1. / 9., 1. / 36.,
                    1. / 36., 1. / 36., 1. / 36.};
const int opp[] = {0, 2, 1, 4, 3, 6, 5, 8, 7};      // Gives the vector pointing opposite to i

// DEFINE FUNCTIONS

// Functions to return the location of data in arrays:
// (this was done to allow timing experiments by reversing locations)
unsigned int scalar(unsigned int x, unsigned int y) { return NY * x + y; }
unsigned int scalar3(unsigned int x, unsigned int y) { return 3 * x + y; }
unsigned int vector(unsigned int x, unsigned int y, unsigned int z) { return 9 * (NY * x + y) + z; }

// Function to assign equilibrium values eq to the given array:
void equilibrium(double eq[], double vx[], double vy[], double r[], bool bound[])
{
#pragma omp parallel for
    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            if (!bound[scalar(x, y)])
            {
                // These three do no depend upon i and so are raised
                //  to the loop above to prevent excess computation
                double dot2 = pow(vx[scalar(x, y)], 2) + pow(vy[scalar(x, y)], 2);
                double vx_ = vx[scalar(x, y)];
                double vy_ = vy[scalar(x, y)];

                for (int i = 0; i < 9; i++)
                {
                    double dot1 = ex[i] * vx_ + ey[i] * vy_;
                    eq[vector(x, y, i)] = r[scalar(x, y)] * w[i] * (1 + (3. * dot1) + (9. * pow(dot1, 2.) / 2.) - (3. * dot2 / 2.));
                }
            }
        }
    }
}

// Function to output data to binary files:
int image(double image[], char name[], unsigned int num)
{
    // Allocate array to memory
    const size_t mem_val = sizeof(double) * (NX * NY + 3);
    double *val = (double *)malloc(mem_val);

    val[0] = (double)NX;
    val[1] = (double)NY;
    val[2] = (double)num;

    #pragma omp parralel for
    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            val[scalar(x, y) + 3] = image[scalar(x, y)];
        }
    }

    // Open, write to and close the file
    FILE *f;
    char str[50];
    sprintf(str, "%somp_%s_%05d.bin", file_path, name, num);

    f = fopen(str, "w");
    if (!f)
        return 1;
    fwrite(val, mem_val, 1, f);
    fclose(f);

    // Remove the array from memory
    free(val);
}

// Function to make a modulus operatorwhich works
//  for negative numbers (turns it into python modulus)
int mod(int x, int m) { return (x % m + m) % m; }

//-----------------------------------------------------\\
//                  MAIN FUNCTION                      \\
//-----------------------------------------------------\\

int main(int argc, char *argv[])
{
    omp_set_num_threads(atoi(argv[1])); // Enters thread number when calling

    const double centre[] = {NX / 4, NY / 2}; // Location of circle
    const double radius = NY / 8;             // Radius of circle

    const double tau = 0.5 + 3 * uw * radius / Re; // Collision parameter

    // Define sizes for memory allocation
    const size_t mem_scalar = sizeof(double) * NX * NY;
    const size_t mem_vector = sizeof(double) * NX * NY * 9;
    const size_t mem_init_vel = sizeof(double) * NY;
    const size_t mem_init_eq = sizeof(double) * NY * 3;
    const size_t mem_bound = sizeof(bool) * NX * NY;

    // Allocate memory to variables
    double *rho = (double *)malloc(mem_scalar);
    double *ux = (double *)malloc(mem_scalar);
    double *uy = (double *)malloc(mem_scalar);
    double *vel = (double *)malloc(mem_scalar);

    double *f = (double *)malloc(mem_vector);
    double *f1 = (double *)malloc(mem_vector);
    double *f_eq = (double *)malloc(mem_vector);

    double *ux0 = (double *)malloc(mem_init_vel);

    bool *boundary = (bool *)malloc(mem_bound);

    // Initialise velocity, pressure and boundary
    #pragma omp parralel for
    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            // Set to ux to uw and add a small amount of noise
            //  (noise helps vortices to form)
            ux[scalar(x, y)] = uw + ((double) rand() / pow(2, 31) - 0.5) * uw / 100;
            uy[scalar(x, y)] = 0;

            rho[scalar(x, y)] = rho_0;

            // Boolean array, 1 within circle, 0 outside
            boundary[scalar(x, y)] = pow((centre[0] - x), 2) + pow((centre[1] - y), 2) <= pow(radius, 2);

            // Save initial velocity for later
            if (x == 0)
                ux0[y] = ux[scalar(0, y)];
        }
    }

    // Initialise f_eq, and f to the equilibrium value
    equilibrium(f, ux, uy, rho, boundary);
    equilibrium(f_eq, ux, uy, rho, boundary);

    // Initialise various timers
    double amdahl_s_time = 0;                // Counts serial times
    double amdahl_p_time = -omp_get_wtime(); // Counts parallel times
    double T = -omp_get_wtime();             // Counts total loop time
    double t = - omp_get_wtime();     // Counts code time

    //-----------------------------------------------------\\
    //                      MAIN LOOP                      \\
    //-----------------------------------------------------\\

    printf("Starting simulation...\n");

    for (int step = 0; step < N + 1; step++)
    {
        // Print progress to command line
        if (step % 10 == 0 & step > 0)
        {
            t =  (t + omp_get_wtime()) / 10.;
            printf("\rstep: %d  loop_time: %.05f", step, t);
            fflush(stdout);
            t = - omp_get_wtime();
        }
        

        // Data/Image output
        if (step % frame_interval == 0)
        {
            #pragma omp parralel for
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    // Make boundary clear to output with NAN
                    if (boundary[scalar(x, y)])
                    {
                        vel[scalar(x, y)] = NAN;
                        rho[scalar(x, y)] = NAN;
                    }
                    
                    // Calculate fluid speed
                    else
                        vel[scalar(x, y)] = sqrt(pow(ux[scalar(x, y)], 2) + pow(uy[scalar(x, y)], 2));
                }
            }

            amdahl_p_time += omp_get_wtime();
            amdahl_s_time -= omp_get_wtime();

            // Data saving is the only serial part to ensure 
            //  no data is corrupted
            image(rho, "rho", (step / frame_interval));
            image(vel, "vel", (step / frame_interval));
            
            amdahl_s_time += omp_get_wtime();
            amdahl_p_time -= omp_get_wtime();
        }

        //Calculate new velocity and pressure
        #pragma omp parallel for
        for (int x = 1; x < NX; x++)
        {
            for (int y = 0; y < NY; y++)
            {
                if (!boundary[scalar(x, y)])
                {
                    double sum_rho = 0;
                    double sum_ux = 0;
                    double sum_uy = 0;

                    for (int i = 0; i < 9; i++)
                    {
                        sum_rho += f[vector(x, y, i)];
                        sum_ux += f[vector(x, y, i)] * ex[i];
                        sum_uy += f[vector(x, y, i)] * ey[i];
                    }

                    rho[scalar(x, y)] = sum_rho;
                    ux[scalar(x, y)] = sum_ux / sum_rho;
                    uy[scalar(x, y)] = sum_uy / sum_rho;

                    // This is talked about in the report.
                    //  It saves time but only when optimisation is not used when compiling.
                    /*
                    rho[scalar(x, y)] = f[vector(x, y, 0)] + f[vector(x, y, 1)] + f[vector(x, y, 2)] +
                                        f[vector(x, y, 3)] + f[vector(x, y, 4)] + f[vector(x, y, 5)] +
                                        f[vector(x, y, 6)] + f[vector(x, y, 7)] + f[vector(x, y, 8)];
                                        
                    ux[scalar(x, y)] = ((f[vector(x, y, 1)] + f[vector(x, y, 5)] + f[vector(x, y, 8)])  -
                                        (f[vector(x, y, 2)] + f[vector(x, y, 6)] + f[vector(x, y, 7)])) / 
                                        rho[scalar(x, y)];
                    
                    uy[scalar(x, y)] = ((f[vector(x, y, 3)] + f[vector(x, y, 5)] + f[vector(x, y, 7)])  -
                                        (f[vector(x, y, 4)] + f[vector(x, y, 6)] + f[vector(x, y, 8)])) / 
                                        rho[scalar(x, y)];

                    */
                }
            }
        }

        #pragma omp parallel for
        for (int y = 0; y < NY; y++)
        {
            // Set left boundary velocity
            ux[scalar(0, y)] = ux0[y];
            uy[scalar(0, y)] = 0.;

            // Calculate left density
            double col1 = f[vector(0, y, 2)] + f[vector(0, y, 6)] + f[vector(0, y, 7)];
            double col2 = f[vector(0, y, 0)] + f[vector(0, y, 3)] + f[vector(0, y, 4)];
            rho[scalar(0, y)] = (col2 + 2 * col1) / (1 - ux[scalar(0, y)]);

            // Right pressure and velocity are not prescribed and are left open
        }

        // Calculate equilibrium distribution functions
        equilibrium(f_eq, ux, uy, rho, boundary);

        // Apply boundary coundions to inlet and outlet
        #pragma omp parallel for
        for (int y = 0; y < NY; y++)
        {
            // Apply full LBM to left boundary - replaces collision and streaming
            f1[vector(0, y, 1)] = f_eq[vector(0, y, 1)] + (f[vector(0, y, opp[1])] - f_eq[vector(0, y, opp[1])]);
            f1[vector(0, y, 5)] = f_eq[vector(0, y, 5)] + (f[vector(0, y, opp[5])] - f_eq[vector(0, y, opp[5])]);
            f1[vector(0, y, 8)] = f_eq[vector(0, y, 8)] + (f[vector(0, y, opp[8])] - f_eq[vector(0, y, opp[8])]);
        
            // Apply full LBM to right boundary - replaces collision and streaming
            f1[vector(NX - 1, y, 2)] = f_eq[vector(NX - 1, y, 2)] + (f[vector(NX - 1, y, opp[2])] - f_eq[vector(NX - 1, y, opp[2])]);
            f1[vector(NX - 1, y, 7)] = f_eq[vector(NX - 1, y, 7)] + (f[vector(NX - 1, y, opp[7])] - f_eq[vector(NX - 1, y, opp[7])]);
            f1[vector(NX - 1, y, 6)] = f_eq[vector(NX - 1, y, 6)] + (f[vector(NX - 1, y, opp[6])] - f_eq[vector(NX - 1, y, opp[6])]);
        }

        // Collision step
        #pragma omp parallel for
        for (int x = 0; x < NX; x++)
            for (int y = 0; y < NY; y++)
                if (!boundary[scalar(x, y)])
                    for (int i = 0; i < 9; i++)
                        f[vector(x, y, i)] = (1 - (1 / tau)) * f[vector(x, y, i)] + (1 / tau) * f_eq[vector(x, y, i)];

        // Stream step
        #pragma omp parallel for
        for (int x = 0; x < NX; x++)
        {
            for (int y = 0; y < NY; y++)
            {
                if (!boundary[scalar(x, y)])        // Don't waste time computing in boundary circle
                {
                    for (int i = 0; i < 9; i++)
                    {
                        // True if distribution function will stream to the boundary circle in next step
                        bool going_to_boundary = boundary[scalar(mod(x + ex[i], NX), mod(y + ey[i], NY))];
                        // True if distribution function will stream to the inlet from outlet
                        //  because of the periodic boundary -- will overwrite loop line ~270
                        bool going_to_inlet = (x == NX - 1 & (i == 1 | i == 5 | i == 8));
                        // True if distribution function will stream to the outlet from inlet
                        //  because of the periodic boundary -- will overwrite loop line ~270
                        bool going_to_outlet = (x == 0 & (i == 2 | i == 7 | i == 6));

                        if (!(going_to_inlet | going_to_outlet))
                        {
                            if (going_to_boundary)
                            {
                                // Apply bounce-back boundary from circle
                                f1[vector(x, y, opp[i])] = f[vector(x, y, i)];
                            }
                            else
                            {
                                // Apply periodic boundary along top and bottom, and stream in centre
                                f1[vector(mod(x + ex[i], NX), mod(y + ey[i], NY), i)] = f[vector(x, y, i)];
                            }
                        }
                    }
                }
            }
        }

        // Swap f1 to f
        double *a = f1;
        f1 = f;
        f = a;
    }

    // Finish timing
    amdahl_p_time += omp_get_wtime();
    T += omp_get_wtime();

    printf("\nSimulation complete\n");
    printf("Total time: %f\n", T);
    printf("Files saved to $%s\n\n", file_path);

    return 0;
}
