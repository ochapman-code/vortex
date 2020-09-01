#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

const int N = 5000;             // Iterations
const int frame_interval = 500; // Iterations between exporting data
const int NX = 520;             // x resolution
const int NY = 180;             // y resolution

const int gate_start = 0;   // First interation to start time gates
const int TC = 50;         // Amount of time gates to pass
                            //  must finish before end of script!

// Exactly ONE of each of the following must be defined:
// MPI task distribution - (see report)
//const char processor_split[15] = "EQUAL_WIDTH";       // Equal MPI computation widths (equal no. of columns)
const char processor_split[15] = "EQUAL_AREA";          // Equal area of computation (adaptive no. of columns)

const char file_path[128] = "mpi_img/"; // Path where files are saved

const double uw = 0.1;          // Initial fluid speed in x
const double Re = 55.;          // Reynolds number
const double rho_0 = 1.;        // Initial pressure

const int ex[] = {0, 1, -1, 0, 0, 1, -1, -1, 1};    // Vector directions of distribution func.
const int ey[] = {0, 0, 0, 1, -1, 1, -1, 1, -1};    // Vector directions of distribution func.

const double w[] = {4. / 9., 1. / 9., 1. / 9.,      // Weights of vectors
                    1. / 9., 1. / 9., 1. / 36.,
                    1. / 36., 1. / 36., 1. / 36.};
const int opp[] = {0, 2, 1, 4, 3, 6, 5, 8, 7};      // Gives the vector pointing opposite to i

// Initialise MPI parameters
unsigned int rank_xstart, rank_xend, rank_xnum;     
int col_start, col_num;
int size, rank;

// DEFINE FUNCTIONS

// Functions to return the location of data in arrays:
// (this was done to allow timing experiments by reversing locations)
unsigned int scalar(unsigned int x, unsigned int y) { return NY * x + y; }
unsigned int scalar3(unsigned int x, unsigned int y) { return 3 * x + y; }
unsigned int vector(unsigned int x, unsigned int y, unsigned int z) { return 9 * (NY * x + y) + z; }

// Function to assign equilibrium values eq to the given array:
void equilibrium(double eq[], double vx[], double vy[], double r[], bool bound[])
{
    for (int x = 0; x < rank_xnum; x++)
    {
        int b_x = rank_xstart + x;
        for (int y = 0; y < NY; y++)
        {
            if (!bound[scalar(b_x, y)])
            {
                // These three do no depend upon i and so are raised
                //  to the loop above to prevent excess computation
                double vx_ = vx[scalar(x, y)];
                double vy_ = vy[scalar(x, y)];
                double dot2 = vx_ * vx_ + vy_ * vy_;

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

    val[0] = (double) NX;
    val[1] = (double) NY;
    val[2] = (double) num;

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
    sprintf(str, "%smpi_%s_%05d.bin", file_path, name, num);

    f = fopen(str, "w");
    if (!f)
        return 1;
    fwrite(val, mem_val, 1, f);
    fclose(f);

    // Remove the array from memory
    free(val);
}

// Function to time task distribution between MPI ranks
void time_gate(int loop_step, int gate_type, double time_log[], int type_log[])
{
    // Types:
    // 0 - Computation
    // 1 - Idle
    // 2 - Internal Communication
    // 3 - External Saving

    static int counts = 0;
    static double counter;

    // If tripped, save time and reset
    if (counts < TC & loop_step >= gate_start)
    {
        time_log[counts] = MPI_Wtime() - counter;
        type_log[counts] = gate_type;
        counts++;
    }

    // Align timing on each processor on loop before using gates
    if (loop_step - gate_start == -1)
        //Only interrupts flow on one loop - minimal disturbance
        MPI_Barrier(MPI_COMM_WORLD);

    counter = MPI_Wtime();
}

// Function to make a modulus operatorwhich works
//  for negative numbers (turns it into python modulus)
int mod(int x, int m) { return (x % m + m) % m; }

//-----------------------------------------------------\\
//                  MAIN FUNCTION                      \\
//-----------------------------------------------------\\

int main()
{
    const double centre[] = {NX / 4, NY / 2};
    const double radius = NY / 8;

    // Initialise boundary circle array
    const size_t mem_bound = sizeof(bool) * NX * NY;
    bool *boundary = (bool *)malloc(mem_bound);

    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            // Boolean array, 1 within circle, 0 outside
            boundary[scalar(x, y)] = pow((centre[0] - x), 2) + pow((centre[1] - y), 2) <= pow(radius, 2);
        }
    }

    // Initialise MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size < 2 & rank == 0)
    {
        printf("The MPI Lattice Boltzmann Simulation must have size >= 2\n");
        MPI_Finalize();
        return 1;
    }

    // Calculate area covered by boundary
    unsigned int sum = 0;
    unsigned int sum_min[NX];
    unsigned int sum_max[NX];

    // Decide on the columns each rank will operate on
    for (int x = 0; x < NX; x++)
    {
        sum_min[x] = sum;
        for (int y = 0; y < NY; y++)
        {
            if (strstr(processor_split, "EQUAL_AREA") != NULL)
            {
                // Count NON-BOUNDARY points
                if (!boundary[scalar(x, y)])
                    sum++;
            }

            else if (strstr(processor_split, "EQUAL_WIDTH") != NULL)
            {
                // Count ALL points
                sum++;
            }

            else if (rank == 0)
            {
                printf("Variable processor_split is not set correctly\n");
                printf("\t- must be set to \"EQUAL_AREA\" or \"EQUAL_WIDTH\" \n");
                printf("\t- currently \"%s\"\n", processor_split);
                MPI_Finalize();
                return 1;
            }
        }
        sum_max[x] = sum;       // Save column totals
    }

    double ideal_sum = (double)sum / (double)size;          // Ideal column width       
    double ideal_min = (double)rank * ideal_sum;            // Ideal 1st column for rank
    double ideal_max = ((double)rank + 1.0) * ideal_sum;    // Ideal 2nd column for rank
                                                            //  - (often fractions) 
    // Assign the range that each rank covers
    for (int x = 0; x < NX; x++)
    {
        if (sum_min[x] <= ideal_min & ideal_min < sum_max[x])
        {
            rank_xstart = x;                    // Actual first columns for rank
            if (ideal_min - (int)ideal_min > 0)
                rank_xstart++;
        }
        if (sum_min[x] < ideal_max & ideal_max <= sum_max[x])
        {
            rank_xend = x;                      // Actual last columns for rank
        }
    }

    rank_xnum = rank_xend - rank_xstart + 1;    // Actual number of columns for rank

    // Share rank info with all other ranks
    int col_starts[size];       // Amount of cells before rank_xstart for each rank
    int col_nums[size];         // Amount of cells between rank_xstart and rank_xend for each rank

    col_start = rank_xstart * NY;
    col_num = rank_xnum * NY;

    MPI_Allgather(&col_start, 1, MPI_INT, col_starts, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&col_num, 1, MPI_INT, col_nums, 1, MPI_INT, MPI_COMM_WORLD);
    // End of MPI initialisation

    const double tau = 0.5 + 3 * uw * radius / Re;      // Collision parameter

    // Define sizes for memory allocation
    const size_t mem_scalar = sizeof(double) * rank_xnum * NY;
    const size_t mem_scalar_full = sizeof(double) * NX * NY;
    const size_t mem_vector = sizeof(double) * rank_xnum * NY * 9;
    const size_t mem_init_vel = sizeof(double) * NY;
    const size_t mem_init_eq = sizeof(double) * NY * 3;

    // Allocate memory to variables
    double *rho = (double *)malloc(mem_scalar);
    double *ux = (double *)malloc(mem_scalar);
    double *uy = (double *)malloc(mem_scalar);
    double *vel = (double *)malloc(mem_scalar);

    double *img_rho = (double *)malloc(mem_scalar_full);
    double *img_vel = (double *)malloc(mem_scalar_full);

    double *f = (double *)malloc(mem_vector);
    double *f1 = (double *)malloc(mem_vector);
    double *f_eq = (double *)malloc(mem_vector);

    double *ux0 = (double *)malloc(mem_init_vel);

    // Initialise velocity and pressure
    for (int x = 0; x < rank_xnum; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            // Set to ux to uw and add a small amount of noise
            //  (noise helps vortices to form)
            ux[scalar(x, y)] = uw + ((double) rand() / pow(2, 31) - 0.5) * uw / 100;
            uy[scalar(x, y)] = 0;
            rho[scalar(x, y)] = 1;

            // Save initial ux for future boundaries
            if (x == 0)
                ux0[y] = ux[scalar(0, y)];
        }
    }

    // Initialise f_eq, and f to the equilibrium value
    equilibrium(f, ux, uy, rho, boundary);
    equilibrium(f_eq, ux, uy, rho, boundary);

    
    // Initialise various timers
    double amdahl_s_time = 0;
    double amdahl_p_time = - MPI_Wtime();

    double T = - MPI_Wtime();     // Counts code time
    double t = - MPI_Wtime();     // Counts code time

    double count_log[TC];       // Logs time gate times
    int count_type[TC];         // Logs time gate time types
    
    // Initialises timer if starting on step = 0
    time_gate(-1, 0, count_log, count_type);

    //-----------------------------------------------------\\
    //                      MAIN LOOP                      \\
    //-----------------------------------------------------\\

    if (rank == 0)
        printf("Starting lattice Boltzman model simulation...\n");

    for (int step = 0; step < N + 1; step++)
    {
        // Print progress to command line
        if (rank == 0)
        {
            if (step % 10 == 0 & step > 0)
            {
                t =  (t + MPI_Wtime()) / 10.;
                printf("\rstep: %d  loop_time: %.05f", step, t);
                fflush(stdout);
                t = - MPI_Wtime();
            }
        }

        // Start timer at start of first applicable loop
        // TIME GATE 0 - Computation
        time_gate(step, 0, count_log, count_type);

        // Data/Image output
        if (step % frame_interval == 0)
        {
            // Create velocity image
            for (int x = 0; x < rank_xnum; x++)
            {
                int b_x = rank_xstart + x;
                for (int y = 0; y < NY; y++)
                {
                    // Make boundary clear to output with NAN
                    if (boundary[scalar(b_x, y)])
                    {
                        rho[scalar(x, y)] = NAN;
                        vel[scalar(x, y)] = NAN;
                    }
                    
                    // Calculate fluid speed
                    else        
                        vel[scalar(x, y)] = sqrt(pow(ux[scalar(x, y)], 2) + pow(uy[scalar(x, y)], 2));
                }
            }

            // TIME GATE 1 - Computation
            time_gate(step, 0, count_log, count_type);
            // Barrier doesn't interrupt flow due to blocking communications, but allows idle timing
            MPI_Barrier(MPI_COMM_WORLD);
            // TIME GATE 2 - Idle
            time_gate(step, 1, count_log, count_type);

            // Gather rhos into img_rho on rank 0
            MPI_Gatherv(rho, rank_xnum * NY, MPI_DOUBLE, img_rho, col_nums, col_starts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // Gather vels into img_vel on rank 1
            MPI_Gatherv(vel, rank_xnum * NY, MPI_DOUBLE, img_vel, col_nums, col_starts, MPI_DOUBLE, 1, MPI_COMM_WORLD);

            // TIME GATE 3 - Communication
            time_gate(step, 2, count_log, count_type);

            if (rank == 0) amdahl_p_time += MPI_Wtime();
            if (rank == 0) amdahl_s_time -= MPI_Wtime();

            // Data saving is the only serial part to ensure 
            //  no data is corrupted
            // Separate files split between two ranks to increae speed
            if (rank == 0)
                image(img_rho, "rho", (step / frame_interval));
            if (rank == 1)
                image(img_vel, "vel", (step / frame_interval));

            // TIME GATE 4 - Saving
            time_gate(step, 3, count_log, count_type);

            if (rank == 0) amdahl_s_time += MPI_Wtime();
            if (rank == 0) amdahl_p_time -= MPI_Wtime();
        }


        // Calculate new velocity and rho
        for (int x = 0; x < rank_xnum; x++)
        {
            int b_x = rank_xstart + x;
            for (int y = 0; y < NY; y++)
            {
                if (!boundary[scalar(b_x, y)])
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
                }
            } 
        }

        // Apply left boundary
        if (rank == 0)
        {
            for (int y = 0; y < NY; y++)
            {
                // Apply left boundary velocity
                ux[scalar(0, y)] = ux0[y];
                uy[scalar(0, y)] = 0.;

                // Calculate left density
                double col1 = f[vector(0, y, 2)] + f[vector(0, y, 6)] + f[vector(0, y, 7)];
                double col2 = f[vector(0, y, 0)] + f[vector(0, y, 3)] + f[vector(0, y, 4)];
                rho[scalar(0, y)] = (col2 + 2 * col1) / (1 - ux[scalar(0, y)]);

                // Right pressure and velocity are not prescribed and are left open
            }
        }

        // Calculate equilibrium distribution functions
        equilibrium(f_eq, ux, uy, rho, boundary);

        for (int y = 0; y < NY; y++)
        {
            // Apply full LBM to left boundary - replaces collision and streaming
            if (rank == 0)
            {
                f1[vector(0, y, 1)] = f_eq[vector(0, y, 1)] + (f[vector(0, y, opp[1])] - f_eq[vector(0, y, opp[1])]);
                f1[vector(0, y, 5)] = f_eq[vector(0, y, 5)] + (f[vector(0, y, opp[5])] - f_eq[vector(0, y, opp[5])]);
                f1[vector(0, y, 8)] = f_eq[vector(0, y, 8)] + (f[vector(0, y, opp[8])] - f_eq[vector(0, y, opp[8])]);
            }

            // Apply full LBM to left boundary - replaces collision and streaming
            if (rank == size - 1)
            {
                f1[vector(rank_xnum - 1, y, 2)] = f_eq[vector(rank_xnum - 1, y, 2)] + (f[vector(rank_xnum - 1, y, opp[2])] - f_eq[vector(rank_xnum - 1, y, opp[2])]);
                f1[vector(rank_xnum - 1, y, 7)] = f_eq[vector(rank_xnum - 1, y, 7)] + (f[vector(rank_xnum - 1, y, opp[7])] - f_eq[vector(rank_xnum - 1, y, opp[7])]);
                f1[vector(rank_xnum - 1, y, 6)] = f_eq[vector(rank_xnum - 1, y, 6)] + (f[vector(rank_xnum - 1, y, opp[6])] - f_eq[vector(rank_xnum - 1, y, opp[6])]);
            }
        }

        // Collision step
        for (int x = 0; x < rank_xnum; x++)
        {
            int b_x = rank_xstart + x;      // Gives global position on boundary array
            for (int y = 0; y < NY; y++)
                if (!boundary[scalar(b_x, y)])
                    for (int i = 0; i < 9; i++)
                        f[vector(x, y, i)] = (1 - (1 / tau)) * f[vector(x, y, i)] + (1 / tau) * f_eq[vector(x, y, i)];
        }

        // Share data between ranks
        MPI_Status status;
        const size_t mem_send = sizeof(double) * 3 * NY;
        double *rightward = (double *)malloc(mem_send);
        double *leftward = (double *)malloc(mem_send);

        if (rank > 0)
        {
            // Fill array with elements heading to rank to the left
            for (int y = 0; y < NY; y++)
            {
                leftward[scalar3(y, 0)] = f[vector(0, y, 2)];
                leftward[scalar3(y, 1)] = f[vector(0, y, 7)];
                leftward[scalar3(y, 2)] = f[vector(0, y, 6)];
            }
        }

        if (rank < size - 1)
        {
            // Fill array with elements heading to rank to the right
            for (int y = 0; y < NY; y++)
            {
                rightward[scalar3(y, 0)] = f[vector(rank_xnum - 1, y, 1)];
                rightward[scalar3(y, 1)] = f[vector(rank_xnum - 1, y, 5)];
                rightward[scalar3(y, 2)] = f[vector(rank_xnum - 1, y, 8)];
            }
        }

        // TIME GATE 5 - Computation
        time_gate(step, 0, count_log, count_type);
        // Barrier doesn't interrupt flow due to blocking communications, but allows idle timing
        MPI_Barrier(MPI_COMM_WORLD);
        // TIME GATE 6 - Idle
        time_gate(step, 1, count_log, count_type);

        //Send rightward elements to rank to the right
        if (rank < size - 1)
            MPI_Send(rightward, NY * 3, MPI_DOUBLE, rank + 1, 999, MPI_COMM_WORLD);

        if (rank > 0)
        {
            //Receive rightward elements from block to the left
            MPI_Recv(rightward, NY * 3, MPI_DOUBLE, rank - 1, 999, MPI_COMM_WORLD, &status);

            //Send leftward elements to rank to the left
            MPI_Send(leftward, NY * 3, MPI_DOUBLE, rank - 1, 111, MPI_COMM_WORLD);
        }

        //Receive leftward elements from block to the right
        if (rank < size - 1)
            MPI_Recv(leftward, NY * 3, MPI_DOUBLE, rank + 1, 111, MPI_COMM_WORLD, &status);

        // TIME GATE 7 - Communication
        time_gate(step, 2, count_log, count_type);

        int mem[9] = {-1, 0, 0, -1, -1, 1, 2, 1, 2}; // -1 in unused positions

        // Streaming
        for (int x = 0; x < rank_xnum; x++)
        {
            int b_x = rank_xstart + x;      // x position in full grid
            for (int y = 0; y < NY; y++)
            {
                if (!boundary[scalar(b_x, y)])      // Don't waste time computing in boundary circle
                {
                    for (int i = 0; i < 9; i++)
                    {
                        // True if distribution function will stream to the boundary circle in next step
                        bool going_to_boundary = boundary[scalar(mod(b_x + ex[i], NX), mod(y + ey[i], NY))];

                        // Bounce-back boundary: from circle
                        if (going_to_boundary)
                            f1[vector(x, y, opp[i])] = f[vector(x, y, i)];

                        // Periodic boundaries: away from border between ranks
                        else if (x + ex[i] >= 0 & x + ex[i] < rank_xnum)
                            f1[vector(x + ex[i], mod(y + ey[i], NY), i)] = f[vector(x, y, i)];

                        // Periodic boundaries: incoming data from process to the left
                        else if (x + ex[i] == -1 & rank > 0)
                            f1[vector(0, y, opp[i])] = rightward[scalar3(mod(y + ey[i], NY), mem[opp[i]])];

                        // Periodic boundaries: incoming data from process to the right
                        else if (x + ex[i] == rank_xnum & rank < size - 1)
                            f1[vector(rank_xnum - 1, y, opp[i])] = leftward[scalar3(mod(y + ey[i], NY), mem[opp[i]])];

                        // Inlet and outlet boundaries are covered around line 416
                    }
                }
            }
        }

        // Swap f1 to f
        double *a = f1;
        f1 = f;
        f = a;
    }

    if (rank == 0) amdahl_p_time += MPI_Wtime();

    // Send timing data to rank 0
    if (rank > 0)
        MPI_Send(count_log, TC, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);

    // Save timing data
    if (rank == 0)
    {
        FILE *f;

        char str[256];
        sprintf(str, "%sgantt_chart_data.txt", file_path);
        f = fopen(str, "w");

        // Check file is open
        if (!f)
        {
            printf("\nCould not open file: %s\n", str);
            printf("Aborting file write");
            MPI_Finalize();
            return 1;
        }

        // Write time types e.g. idle, communicating, computation...
        char log_types[5] = "grbc";
        char letter[5];
        for (int i = 0; i < TC; i++)
        {
            sprintf(letter, "%c", log_types[count_type[i]]);
            fwrite(letter, sizeof(char), 1, f);
        }
        fwrite("\n", sizeof(char), 1, f);

        // Write times from rank 0
        for (int i = 0; i < TC; i++)
        {
            sprintf(str, "%.08f\t", count_log[i]);
            fwrite(str, sizeof(char), 11, f);
        }

        // Write times from other ranks
        for (int i = 1; i < size; i++)
        {
            MPI_Status status;
            MPI_Recv(count_log, TC, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);

            fwrite("\n", sizeof(char), 1, f);
            for (int i = 0; i < TC; i++)
            {
                sprintf(str, "%.08f\t", count_log[i]);
                fwrite(str, sizeof(char) * 11, 1, f);
            }
        }

        T += MPI_Wtime();

        if (rank == 0)
        {
            printf("\nSimulation complete\n");
            printf("Total time: %f\n", T);
            printf("Files saved to $%s\n\n", file_path);
        }

        // Close text file
        fclose(f);
    }

    MPI_Finalize();

    return 0;
}
