#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include "barrier.h"
#include <time.h>
#include <sys/time.h>
#include "ptools.h"

struct Matrix_t {
    int threadId;
    pthread_barrier_t *barrier;
};

int N;
int NUM_THREADS;
pthread_mutex_t lock;
float *A, *B;

int main(int argc, char *argv[]){
    double duration;
    N = atoi(argv[1]);
    NUM_THREADS = atoi(argv[2]);
    struct timeval start, stop;
    B = (float*)malloc(N*N*sizeof(float));
    A = generateMatrix(N);
    memcpy(B, A, N*N*sizeof(float));

    gettimeofday(&start, 0);
    setThreads(A);
    gettimeofday(&stop, 0);
    duration = (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) * 1e-6;
    printf("Pthread Gaussian Elimination takes %lf seconds \n", duration);

    gettimeofday(&start, 0);
    serialGaussianEliminationOptimized(N, B);
    gettimeofday(&stop, 0);
    duration = (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) * 1e-6;
    printf("Serial Gaussian Elimination with tile optimization takes %lf seconds \n", duration);

    validateMatrix(A, B, N);

    free(A);
    free(B);
    return 0;
}

void setThreads(float* A) {
    pthread_t threads[NUM_THREADS];
    pthread_barrier_t barrier;
    struct Matrix_t matrixThreads[NUM_THREADS];
    pthread_barrier_init(&barrier, NULL, NUM_THREADS);

    for(int i=0; i<NUM_THREADS; i++) {
        matrixThreads[i].threadId = i;
        matrixThreads[i].barrier = &barrier;
        if(pthread_create(&threads[i], NULL, gaussPthreads, (void*)&matrixThreads[i]) != 0) {
            perror("ERROR: fail to create thread...\n");
        }
    }
    for(int i=0; i<NUM_THREADS; i++) {
        if(pthread_join( threads[i], NULL) != 0) {
            perror("ERROR: fail to join thread ...\n");
        }
    }
    pthread_barrier_destroy(&barrier);
}

void *gaussPthreads(void *args){
    struct Matrix_t *localArgs = (struct Matrix_t*)args;
    int threadId = localArgs->threadId, localId;
    pthread_barrier_t *barrier = localArgs->barrier;

    for(int i=0; i<N; i++) {
        localId = i % NUM_THREADS;
        if(localId == threadId) {
            normaliseSingleRow(i, i, N, A);
        }
        pthread_barrier_wait(barrier);

        for(int j=i+1; j<N; j++) {
            if((j % NUM_THREADS) == threadId) {
                eliminateSingleRow(i, j, N, A);
            }
        }
    }
    pthread_barrier_wait(barrier);

    return NULL;
}
