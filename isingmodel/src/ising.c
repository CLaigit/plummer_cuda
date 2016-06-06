/*
Ising model: Halmitonian H = /sum_ij J(sigma_i)(sigma_j)
We set J = 1 first
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>       /* time */

#define  LATTICE_LENGTH 200
#define  BOLTZMANN_CONST 1
#define  WARMSTEPS 1e3
#define  NSWEEPS 1e3

// Calculate the energy of the (up, center) (down, center) (left, center) ( right, center)
double energy(int up, int down, int left, int right, int center){
    return -center * (up + down + left + right);
}

int main (int argc, char *argv[]){

    static int lattice[LATTICE_LENGTH][LATTICE_LENGTH] = {};

    double T = 2;
    int col, row;
    int new;
    double deltaE = 0.0;
    double tmpE = 0.0, tmpE2 = 0.0, averE = 0.0, averE2 = 0.0;
    double tmpmag = 0.0, tmpmag2 = 0.0, avermag = 0.0, avermag2 = 0.0;
    double siteE = 0.0;

    T = argc > 1 ? atof(argv[1]) : 2;
    col = argc > 2 ? atoi(argv[2]) : 20;
    row = col;

    double beta = 1.0 / BOLTZMANN_CONST / T;

    // Tempurature
    srand (time(NULL));
    // Initialize every grid point
    for (int i = 0; i < col; i++){
        for (int j = 0; j < row; j++){
            lattice[i][j] = 2 * (rand() % 2) - 1;
            // lattice[i][j] = 1;
        }
    }
    // Warmup process
    for (int nstep = 0; nstep < WARMSTEPS; nstep++){
        for (int i = 0; i < col; i++){
            for (int j = 0; j < row; j++){
                // flip a spin
                // If the energy becomes small, accept it
                // Else accept w.r.t the probability e^-beta*deltaE
                new = -lattice[i][j];
                deltaE = energy(lattice[ (i - 1 + row) % row][j], lattice[(i + 1 + row) % row][j], lattice[i][(j - 1 + row) % row], lattice[i][(j + 1 + row) % row], new);
                deltaE -= energy(lattice[ (i - 1 + row) % row][j], lattice[(i + 1 + row) % row][j], lattice[i][(j - 1 + row) % row], lattice[i][(j + 1 + row) % row], lattice[i][j]);
                if (deltaE < 0 || (double)rand() / (double)RAND_MAX <= exp(- beta * deltaE)){
                    lattice[i][j] = new;
                }
            }
        }
    }

    // Measure steps
    for (int nstep = 0; nstep < NSWEEPS; nstep++){
        for (int i = 0; i < col; i++){
            for (int j = 0; j < row; j++){
                new = -lattice[i][j];
                deltaE = energy(lattice[ (i - 1 + row) % row][j], lattice[(i + 1 + row) % row][j], lattice[i][(j - 1 + row) % row], lattice[i][(j + 1 + row) % row], new);
                deltaE -= energy(lattice[ (i - 1 + row) % row][j], lattice[(i + 1 + row) % row][j], lattice[i][(j - 1 + row) % row], lattice[i][(j + 1 + row) % row], lattice[i][j]);
                if (deltaE < 0 || (double)rand() / (double)RAND_MAX <= exp(- beta * deltaE)){
                    lattice[i][j] = new;
                }
            }
        }
        tmpE = 0, tmpE2 = 0, tmpmag = 0, tmpmag2 = 0;
        for (int i = 0; i < col; i++){
            for (int j = 0; j < row; j++){
                siteE = energy(lattice[ (i - 1 + row) % row][j], lattice[(i + 1 + row) % row][j], lattice[i][(j - 1 + row) % row], lattice[i][(j + 1 + row) % row], lattice[i][j])/2.0;
                tmpE += siteE;
                tmpE2 += siteE * siteE;
                tmpmag += lattice[i][j];
                tmpmag2 += lattice[i][j] * lattice[i][j];
            }
        }
        averE += (1.0 * tmpE / col / row / NSWEEPS);
        averE2 += (1.0 * tmpE2 / col / row / NSWEEPS);
        avermag += (1.0 * tmpmag / col / row / NSWEEPS);
        avermag2 += (1.0 * tmpmag2 /col / row / NSWEEPS);
        // Output data every NOUT
#ifdef PICORE
        for (int i = 0; i < col; i++){
            for (int j = 0; j < col-1; j++){
                printf("%d,", lattice[i][j]);
            }
            printf("%d\n", lattice[i][col-1]);
        }
    }
#else
    }
    printf("%f\n", T);
    printf("%d\n", col);
    printf("%f\n", averE);
    printf("%f\n", 1.0*(averE2 - averE * averE) / T / T);
    printf("%f\n", fabs(avermag));
    printf("%f\n", 1.0*(avermag2 - avermag * avermag) / T );
#endif
}
