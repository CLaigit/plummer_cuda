/*
 * LEAPINT.C: program to integrate hamiltonian system using leapfrog.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // rand(), srand()
#include <time.h>       /* time */

// #define GM 0.00029632889       /*4pi/365^2*/
// #define Gmearth 8.8987655E-10
// #define Gmsaturn 8.46797626E-8
// #define Gmpartical 8.8987655E-11
// #define PI 3.1415

#define BLOCK_SIZE 256
#define GM 1
#define PI 3.14159265359
#define NUM_PLANET 1000
#define DT 0.1
#define SOFTPARAMETER 0.00001

typedef struct Vec_3{
    double x;
    double y;
    double z;
}Vec_3;

typedef struct Planet{          /*define a structure to store the position, velocity and dt for a planet*/
  Vec_3* pos;            // WHY USE POINTER INSTEAD OF ACTUAL DATA?? Because I am using
  Vec_3* vel;            // I THINK IF U CHANGE THIS, IT MIGHT WORK
} Planet;

void initialize(Planet *planet);					/* number of points         */
void printstate();				/* print out system state   */



void initialize(Planet *planet)
{
    int i;
    double x1, x2, x3, x4, x5, x6, x7;
    double radius, vra;
    double q;

    srand (time(NULL));
                                /* set number of planet  */
/*Mercury*/
    for (i = 0; i < NUM_PLANET; i++){
        x1 = (double)rand() / (double)RAND_MAX;
        x2 = (double)rand() / (double)RAND_MAX;					/* set initial position */
        x3 = (double)rand() / (double)RAND_MAX;
        x4 = (double)rand() / (double)RAND_MAX;
        x5 = (double)rand() / (double)RAND_MAX;
        x6 = (double)rand() / (double)RAND_MAX;
        x7 = (double)rand() / (double)RAND_MAX;
        radius =  pow( (pow(x1, (-2.0/3.0)) - 1), -0.5 );
        planet[i].pos->z =  (1.0 - 2.0*x2) * radius;
        planet[i].pos->x =  pow( radius*radius - planet[i].pos->z*planet[i].pos->z, 0.5 ) * cos(2.0 * PI * x3);
        planet[i].pos->y =  pow( radius*radius - planet[i].pos->z*planet[i].pos->z, 0.5 ) * sin(2.0 * PI * x3);

        planet[i].vel->x =  0;
        planet[i].vel->y =  0;
        planet[i].vel->z =  0;
    }
}

int main(int argc, char **argv)
{
    int i, n, mstep, nout, nstep;
    double q;
    double x1, x2, x3, x4, x5, x6, x7;

    mstep = 2;			/* number of steps to take  */
    nout = 1;
    	    /* steps between outputs    */
    const int bytes = 2 * NUM_PLANET * sizeof(Vec_3);
    Planet planet[NUM_PLANET];
    for(i = 0; i < NUM_PLANET; i++){
        planet[i].pos = (Vec_3 *)malloc(sizeof(Vec_3));
        planet[i].vel = (Vec_3 *)malloc(sizeof(Vec_3));
    }

    initialize(planet);

    printstate(planet);		/* then output last step    */

}


/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

/*
 * ACCEL: compute accelerations for harmonic oscillator(s).
 */



/*
 * PRINTSTATE: output system state variables.
 */

void printstate(Planet *planet)					/* number of points         */
{
    int i;

    for (i = 0; i < NUM_PLANET; i++)			/* loop over all points...  */
      	printf("%d,%12.6f,%12.6f,%12.6f\n", i, planet[i].pos->x, planet[i].pos->y,  planet[i].pos->z);
}
