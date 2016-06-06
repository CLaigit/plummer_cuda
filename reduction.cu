#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>       /* time */
#include <curand.h>
#include <curand_kernel.h>

#define  LATTICE_LENGTH 256
#define  LATTICE_2 (LATTICE_LENGTH * LATTICE_LENGTH)
#define  N LATTICE_LENGTH
//
// __global__ void reduce0(int *g_idata, int *g_odata)
// {
//     extern __shared__ int sdata[];
//     // each thread loads one element from global to shared mem
//     unsigned int tidx = threadIdx.x;
//     unsigned int tidy = threadIdx.y;
//
//     unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
//     unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
//
//     sdata[tidx] = g_idata[i + j * N];
//     __syncthreads();
//     // do reduction in shared mem
//     for(unsigned int s=1; s < blockDim.x; s *= 2) {
//         if (tid % (2*s) == 0) {
//             sdata[tid] += sdata[tid + s];
//         }
//         __syncthreads();
//     }
//
//     // write result for this block to global mem
//     if (tid == 0) g_odata[blockIdx.x] = sdata[0];
//
//     printf("die\n");
// }

__global__ void test(int *input, int *output){
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

    if (idx < N && idy < N){
        output[idx + idy * N] = input[idx + idy * N];
    }
}


__global__ void printstate(int* output) {
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

    if (idx < N && idy < N){
        printf("%d, %d, %d\n", idx, idy, output[idx + idy * N]);
    }
}

int main (int argc, char *argv[]){

    int *input;
    int *d_input;

    int *output;
    int *d_output;

    int numthreadx = 16;
    int numthready = 16;
    int numblocksX = LATTICE_LENGTH / numthreadx;
    int numblocksY = LATTICE_LENGTH / numthready;

    input = (int*)malloc(LATTICE_2 * sizeof(int));
    output = (int*)malloc(LATTICE_2 * sizeof(int));

    for(int i = 0; i < LATTICE_2; i++){
        input[i] = 1;
        output[i] = 0;
    }

    const size_t bytes_input = LATTICE_2 * sizeof(int);
    const size_t bytes_output = LATTICE_2 * sizeof(int);

    // Set dimensions of block and grid
    dim3 grid(numblocksX, numblocksY, 1);
    dim3 thread(numthreadx, numthready,1);

    // Allocate memoery in device and copy from host to device
    cudaMalloc((void **)&d_input, bytes_input);
    cudaMalloc((void **)&d_output, bytes_output);

    cudaMemcpy(d_input, input, bytes_input, cudaMemcpyHostToDevice);
    cudaMemcpy(d_output, output, bytes_output, cudaMemcpyHostToDevice);

    cudaDeviceSetLimit(cudaLimitPrintfFifoSize, N * N * sizeof(int) * N);

    cudaMemcpy(output, d_output, bytes_output, cudaMemcpyDeviceToHost);

    int sum = 0;
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N; j++){
            sum += output[i + j * N];
        }
    }
    printf("%d\n", sum);
    // printstate<<<grid, thread>>>(d_output);
    cudaDeviceSynchronize();

    free(input);
    cudaFree(d_input);

    free(output);
    cudaFree(d_output);
}
