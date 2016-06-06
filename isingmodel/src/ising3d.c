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
double energy(int up, int down, int left, int right, int front, int back, int center){
    return -center * (up + down + left + right + front + back);
}

int main (int argc, char *argv[]){

    static int lattice[LATTICE_LENGTH][LATTICE_LENGTH][LATTICE_LENGTH] = {};

    double T = 2;
    int col, row, zaxis;
    int new;
    double deltaE = 0.0;
    double tmpE = 0.0, tmpE2 = 0.0, averE = 0.0, averE2 = 0.0;
    double tmpmag = 0.0, tmpmag2 = 0.0, avermag = 0.0, avermag2 = 0.0;
    double siteE = 0.0;

    T = argc > 1 ? atof(argv[1]) : 2;
    col = argc > 2 ? atoi(argv[2]) : 20;
    row = col;
    zaxis = col;
    double beta = 1.0 / BOLTZMANN_CONST / T;
    // Tempurature
    srand (time(NULL));
    // Initialize every grid point
    for (int i = 0; i < col; i++){
        for (int j = 0; j < row; j++){
            for (int k = 0; k < zaxis; k++){
                lattice[i][j][k] = 2 * (rand() % 2) - 1;
                // lattice[i][j][k] = 1;
            }
        }
    }
    // Warmup process
    for (int nstep = 0; nstep < WARMSTEPS; nstep++){
        for (int i = 0; i < col; i++){
            for (int j = 0; j < row; j++){
                for (int k = 0; k < zaxis; k++){
                    // flip a spin
                    // If the energy becomes small, accept it
                    // Else accept w.r.t the probability e^-beta*deltaE
                    new = -lattice[i][j][k];
                    deltaE = energy(lattice[ (i - 1 + row) % row][j][k], lattice[(i + 1 + row) % row][j][k], lattice[i][(j - 1 + row) % row][k], lattice[i][(j + 1 + row) % row][k], lattice[i][j][(k - 1 + zaxis) % zaxis], lattice[i][j][(k + 1 + zaxis) % zaxis],  new);
                    deltaE -= energy(lattice[ (i - 1 + row) % row][j][k], lattice[(i + 1 + row) % row][j][k], lattice[i][(j - 1 + row) % row][k], lattice[i][(j + 1 + row) % row][k], lattice[i][j][(k - 1 + zaxis) % zaxis], lattice[i][j][(k + 1 + zaxis) % zaxis], lattice[i][j][k]);
                    if ((double)rand() / (double)RAND_MAX <= exp(- beta * deltaE)){
                        lattice[i][j][k] = new;
                    }
                }
            }
        }
    }

    // Measure steps
    for (int nstep = 0; nstep < NSWEEPS; nstep++){
        for (int i = 0; i < col; i++){
            for (int j = 0; j < row; j++){
                for (int k = 0; k < zaxis; k++){
                    new = -lattice[i][j][k];
                    deltaE = energy(lattice[ (i - 1 + row) % row][j][k], lattice[(i + 1 + row) % row][j][k], lattice[i][(j - 1 + row) % row][k], lattice[i][(j + 1 + row) % row][k], lattice[i][j][(k - 1 + zaxis) % zaxis], lattice[i][j][(k + 1 + zaxis) % zaxis],  new);
                    deltaE -= energy(lattice[ (i - 1 + row) % row][j][k], lattice[(i + 1 + row) % row][j][k], lattice[i][(j - 1 + row) % row][k], lattice[i][(j + 1 + row) % row][k], lattice[i][j][(k - 1 + zaxis) % zaxis], lattice[i][j][(k + 1 + zaxis) % zaxis], lattice[i][j][k]);
                    if (deltaE < 0 || (double)rand() / (double)RAND_MAX <= exp(- beta * deltaE)){
                        lattice[i][j][k] = new;
                    }
                }
            }
        }
        tmpE = 0, tmpE2 = 0, tmpmag = 0, tmpmag2 = 0;
        for (int i = 0; i < col; i++){
            for (int j = 0; j < row; j++){
                for (int k = 0; k < zaxis; k++){
                    siteE = energy(lattice[ (i - 1 + row) % row][j][k], lattice[(i + 1 + row) % row][j][k], lattice[i][(j - 1 + row) % row][k], lattice[i][(j + 1 + row) % row][k], lattice[i][j][(k - 1 + zaxis) % zaxis], lattice[i][j][(k + 1 + zaxis) % zaxis], lattice[i][j][k])/3;
                    tmpE += siteE;
                    tmpE2 += siteE * siteE;
                    tmpmag += lattice[i][j][k];
                    tmpmag2 += lattice[i][j][k] * lattice[i][j][k];
                }
            }
        }
        averE += (1.0 * tmpE / col / row / zaxis / NSWEEPS);
        averE2 += (1.0 * tmpE2 / col / row / zaxis / NSWEEPS);
        avermag += (1.0 * tmpmag / col / row / zaxis / NSWEEPS);
        avermag2 += (1.0 * tmpmag2 /col / row / zaxis / NSWEEPS);
        // Output data every NOUT
#ifdef PICORE
        for (int i = 0; i < col; i++){
            for (int j = 0; j < row; j++){
                for (int k = 0; k < zaxis; k++){
                    printf("%d,", lattice[i][j][k]);
                }
                printf("%d\n", lattice[i][col-1][zaxis-1]);
            }
            // printf("%d\n", lattice[i][col-1][k]);
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
