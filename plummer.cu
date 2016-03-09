/*
 * LEAPINT.C: program to integrate hamiltonian system using leapfrog.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // rand(), srand()
#include <time.h>       /* time */

#define MAXPNT  100				/* maximum number of points */
// #define GM 0.00029632889       /*4pi/365^2*/
// #define Gmearth 8.8987655E-10
// #define Gmsaturn 8.46797626E-8
// #define Gmpartical 8.8987655E-11
// #define PI 3.1415
#define BLOCK_SIZE 256
#define GM 1
#define PI 3.14159265359
#define NUM_PLANET 1024
#define DT 0.1
#define SOFTPARAMETER 0.00001



typedef struct Planet{          /*define a structure to store the position, velocity and dt for a planet*/
  double* pos;
  double* vel;
} Planet;

void initialize(Planet *planet);
__global__ void leap_step(double *pos, double *vel);
__global__ void accel(double *pos, double *vel);
__global__ void printstate(double *pos);					/* number of points         */


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
        planet[i].pos[2] =  (1.0 - 2.0*x2) * radius;
        planet[i].pos[0] =  pow( radius*radius - planet[i].pos[2]*planet[i].pos[2], 0.5 ) * cos(2.0 * PI * x3);
        planet[i].pos[1] =  pow( radius*radius - planet[i].pos[2]*planet[i].pos[2], 0.5 ) * sin(2.0 * PI * x3);


        // while(0.1 * x5 >= ( x4 * x4 * pow((1 - x4 * x4), 3.5) )  ){
        //     x4 = (double)rand() / (double)RAND_MAX;
        //     x5 = (double)rand() / (double)RAND_MAX;
        // }
        //
        // q = x4;
        //
        //
        // vra =  q * pow(2.0, 0.5) * pow(1 + planet[i]->r * planet[i]->r, -0.25);					/* set initial position */
        // planet[i]->vel[2] =  (1.0 - 2.0 * x6) * vra;
        // planet[i]->vel[0] =  pow( vra * vra - planet[i]->vel[2] * planet[i]->vel[2], 0.5) * cos(2.0 * PI * x7);
        // planet[i]->vel[1] =  pow( vra * vra - planet[i]->vel[2] * planet[i]->vel[2], 0.5) * sin(2.0 * PI * x7);

        planet[i].vel[0] =  0;
        planet[i].vel[1] =  0;
        planet[i].vel[2] =  0;
    }
}

int main(int argc, char **argv)
{
    int mstep, nout, nstep;

    mstep = 100;			/* number of steps to take  */
    nout = 1;		    /* steps between outputs    */

    /* first, set up initial conditions */

    int bytes = 2 * NUM_PLANET * sizeof(double);

    double *buf = (double*)malloc(bytes);
    double *d_buf;
    cudaMalloc(&d_buf, bytes);

    Planet planet = (Planet){(double*)buf, ((double*)buf) + NUM_PLANET};
    Planet d_planet = (Planet){(double*)d_buf, ((double*)d_buf) + NUM_PLANET};
    const unsigned long nBlocks = (NUM_PLANET + BLOCK_SIZE - 1)/BLOCK_SIZE;

    initialize(&planet);

    cudaDeviceSetLimit(cudaLimitPrintfFifoSize, NUM_PLANET * 512);

    cudaMemcpy(d_buf, buf, bytes, cudaMemcpyHostToDevice);

    for (nstep = 0; nstep < mstep; nstep++){	/* loop mstep times in all  */
	       if (nstep % nout == 0)			/* if time to output state  */
	           printstate<<<nBlocks, BLOCK_SIZE>>>(d_planet.pos);		/* then call output routine */
           cudaDeviceSynchronize();
           accel<<<nBlocks, BLOCK_SIZE>>>(d_planet.pos, d_planet.vel);
           cudaDeviceSynchronize();
           leap_step<<<nBlocks, BLOCK_SIZE>>>(d_planet.pos, d_planet.vel);
           cudaDeviceSynchronize();
           accel<<<nBlocks, BLOCK_SIZE>>>(d_planet.pos, d_planet.vel);
           cudaDeviceSynchronize();
    }
    if (mstep % nout == 0)			/* if last output wanted    */
        printstate<<<nBlocks, BLOCK_SIZE>>>(d_planet.pos);

    free(buf);
    cudaFree(d_buf);
}


/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */
__global__ void accel(double *pos, double *vel){
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int tdx = threadIdx;

    double ax = 0.0, ay = 0.0, az = 0.0;
    double d_x = *pos[i].pos, d_y = *(pos[i].pos + 1), d_z = *(pos[i].pos + 2);
    double norm;
    int j, k;

    __shared__ double sx[BLOCK_SIZE];
    __shared__ double sy[BLOCK_SIZE];
    __shared__ double sz[BLOCK_SIZE];

    for(j = 0; j < gridDim.x; j++){
        sx[tdx] = pos[j * BLOCK_SIZE + tdx].pos[0];
        sy[tdx] = pos[j * BLOCK_SIZE + tdx].pos[1];
        sz[tdx] = pos[j * BLOCK_SIZE + tdx].pos[2];
        __syncthreads();

        for(k = 0; k < BLOCK_SIZE; k++){
            norm = pow(SOFTPARAMETER + (d_x - sx[k]) * (d_x - sx[k]) + (d_y - sy[k]) * (d_y - sy[k]) + (d_z - sz[k]) * (d_z - sz[k]), 1.5 );
            ax -= (1.0 * GM/NUM_PLANET) * (d_x - sx[k]) / norm;
            ay -= (1.0 * GM/NUM_PLANET) * (d_y - sy[k]) / norm;
            az -= (1.0 * GM/NUM_PLANET) * (d_z - sz[k]) / norm;
        }
        __syncthreads();
    }

    if(i < NUM_PLANET){
        vel[i].vel[0] += 0.5 * DT * ax;
        vel[i].vel[1] += 0.5 * DT * ay;
        vel[i].vel[2] += 0.5 * DT * az;
    }
}

__global__ void leap_step(double *pos, double *vel){
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < NUM_PLANET){
        pos[i].pos[0] += DT * vel[i].vel[0];
        pos[i].pos[1] += DT * vel[i].vel[1];
        pos[i].pos[2] += DT * vel[i].vel[2];
    }
}

/*
 * PRINTSTATE: output system state variables.
 */

__global__ void printstate(double *pos)					/* number of points         */
{
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < NUM_PLANET){		/* loop over all points...  */
      	printf("%lu,%12.6f,%12.6f,%12.6f\n", i, pos[i].pos[0], pos[i].pos[1], pos[i].pos[2]);
    }
}
