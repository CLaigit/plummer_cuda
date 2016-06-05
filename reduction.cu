/*
Ising model: Halmitonian H = /sum_ij J(sigma_i)(sigma_j)
*/

/*
* TODO:
*   1. Calculate the energy in the program
*   2. Calculate the heat capacity in the program
*   3. Add more inputs to adjust the length of input
*   4. A matlab code to plot data.
*       data format example:
*                    position.x  position.y   spin(-1, 1)
*       Iteattion 1:    1           4               -1
*                       *           *                *
*                       *           *                *
*       Iteattion 2:    4           3                1
*                       *           *                *
*                       *           *                *
*       Iteattion N:    35          76               1
*                       *           *                *
*                       *           *                *
*   5. Compare the numerical value with the analytic value
*   6. Move to 3D
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>       /* time */
#include <curand.h>
#include <curand_kernel.h>

#define  LATTICE_LENGTH 256
#define  LATTICE_2 (LATTICE_LENGTH * LATTICE_LENGTH)
#define  BOLTZMANN_CONST 1
#define  N LATTICE_LENGTH

__global__ void reduce0(int *g_idata, int *g_odata)
{
    extern __shared__ int sdata[];
    // each thread loads one element from global to shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    sdata[tid] = g_idata[i];
    __syncthreads();
    // do reduction in shared mem
    for(unsigned int s=1; s < blockDim.x; s *= 2) {
        if (tid % (2*s) == 0) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

__global__ void printstate(int* output) {
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

    if (idx < N && idy < N){
        printf("%d, %d, %d\n", idx, idy, lattice[idx + idy * N]);
    }
}

int main (int argc, char *argv[]){

    int *input;
    int *d_input;

    int *output;
    int *d_output;

    double T = 2;
    int warmsteps = 1e3;
    int nout;
    nout = 1e5;
    int warp = 1e3;
    int numthreadx = 16;
    int numthready = 16;
    int numblocksX = LATTICE_LENGTH / numthreadx;
    int numblocksY = LATTICE_LENGTH / numthready;

    input = (int*)malloc(LATTICE_2 * sizeof(int));
    output = (int*)malloc(LATTICE_2 * sizeof(int));

    for(int i = 0; i < LATTICE_2; i++){
        input[i] = i;
        output[i] = 0;
    }

    // Set dimensions of block and grid
    dim3 grid(numblocksX, numblocksY, 1);
    dim3 thread(numthreadx, numthready,1);

    // beta is a parameter in the probability
    double beta = 1.0 / BOLTZMANN_CONST / T;

    // Allocate memoery in device and copy from host to device
    cudaMalloc((void **)&d_input, bytes_input);
    cudaMalloc((void **)&d_energy, bytes_input);

    cudaMemcpy(d_input, input, bytes_input, cudaMemcpyHostToDevice);

    cudaDeviceSetLimit(cudaLimitPrintfFifoSize, N * N * sizeof(int) * N);

    reduce0<<<grid, thread>>>(d_input, d_output);
    printstate<<<grid, thread>>>(d_output);

    free(input);
    cudaFree(d_input);
}
