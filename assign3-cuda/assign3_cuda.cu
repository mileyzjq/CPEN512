#include "ctools.h"

__global__ void gaussianNormalization(float *A, int N, int pivot) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    float tmp = A[pivot*N+tid];
    if(tid > pivot && tid < N) {
        tmp /= A[pivot*N+pivot];
    }
    A[pivot*N+tid] = tmp;
}

__global__ void setToOne(float* A, int N, int pivot) {
    A[pivot*N+pivot] = 1.0;
}

__global__ void gaussianElimination(float *A, int N, int pivot) {
    int rowId = blockIdx.x * blockDim.x + threadIdx.x;
    if(rowId > pivot && rowId < N) {
        float scale = A[rowId*N+pivot];
        for(int i=pivot+1; i<N; i++) {
            A[rowId*N+i] -= scale*A[pivot*N+i];
        }
        A[rowId*N+pivot] = 0;
    }
}

__global__ void gaussianElimination2(float *A, int N, int pivot) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    float scale = A[row*N+pivot];
    if(row > pivot && row < N) {
        if(col>=pivot && col<N) {
            A[row*N+col] -= scale*A[pivot*N+col];
        }
    }
}

int main(int argc, char *argv[]) {
    int N = atoi(argv[1]);
    struct timeval begin, end;
    gettimeofday(&begin, 0);
    int NUM_THREADS = atoi(argv[2]);
    int NUM_BLOCKS = (N + NUM_THREADS - 1) / NUM_THREADS;
    dim3 threads(NUM_THREADS, NUM_THREADS);
    dim3 blocks(NUM_BLOCKS, NUM_BLOCKS);
    size_t bytes = N * N * sizeof(float);

    float *h_a;
    cudaMallocHost(&h_a, bytes);

    h_a = generateMatrix(N);

    float *d_a;
    cudaMalloc(&d_a, bytes);
    cudaMemcpy(d_a, h_a, bytes, cudaMemcpyHostToDevice);

    for(int i=0; i<N; i++) {
        gaussianNormalization<<<NUM_BLOCKS, NUM_THREADS>>>(d_a, N, i);
        cudaDeviceSynchronize();
        setToOne<<<1, 1>>>(d_a, N, i);
        gaussianElimination2<<<blocks, threads>>>(d_a, N, i);
        cudaDeviceSynchronize();
    }

    gettimeofday(&end, 0);
    double duration = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec) * 1e-6;
    printf("Cuda Gaussian Elimination takes %lf seconds\n", duration);
    cudaMemcpy(h_a, d_a, bytes, cudaMemcpyDeviceToHost);

    //printMatrix(h_a, N, N);

    cudaFree(d_a);

    return 0;
}