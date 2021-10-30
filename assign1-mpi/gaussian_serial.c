#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "tools.h"

#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define TILE_SIZE 64

void serialGaussianElimination1(int n, float *matrix);
void serialGaussianEliminationOptimized(int n, float *matrix);
void serialGaussianEliminationGroup(int n, float *matrix, int group);
void serialGaussianEliminationBetween(int n, float *matrix, int i, int group);

int main(int argc, char *argv[]) {
    clock_t start, stop;
    double duration;
    int n = atoi(argv[1]);
    float* matrix = generateMatrix(n);
    float* matrix2 = (float*)malloc(n*n*sizeof(float));
    memcpy(matrix2, matrix, n*n*sizeof(float));

    start = clock();
    serialGaussianElimination(n, matrix);
    stop = clock();
    duration = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("Serial Gaussian Elimination time: %lf seconds\n", duration);

    start = clock();
    serialGaussianEliminationOptimized(n, matrix2);
    stop = clock();
    duration = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("Serial Gaussian Elimination with Tile Optimization time: %lf seconds\n", duration);

    validateMatrix(matrix, matrix2, n);
    return 0;
}

void serialGaussianEliminationOptimized(int n, float* matrix) {
    int groups = n / TILE_SIZE;
    for (int i = 0; i < groups; i++) {
        serialGaussianEliminationGroup(n, matrix, i);
    }
}

void serialGaussianEliminationGroup(int n, float* matrix, int group) {
    for (int i = 0; i < group; i++) {
        serialGaussianEliminationBetween(n, matrix, i, group);
    }

    for (int pivot = group*TILE_SIZE; pivot < (group+1)*TILE_SIZE; pivot++) {
        float scale = matrix[pivot * n + pivot];
        for (int i = pivot + 1; i < n; i++) {
            matrix[pivot * n + i] /= scale;
        }
        matrix[pivot * n + pivot] = 1;
        for (int row = pivot + 1; row < (group+1)*TILE_SIZE; row++) {
            float ratio = matrix[row * n + pivot];
            for (int col = pivot + 1; col < n; col++) {
                matrix[row * n + col] -= matrix[pivot * n + col] * ratio;
            }
            matrix[row * n + pivot] = 0;
        }
    }
}

void serialGaussianEliminationBetween(int n, float *matrix, int i, int group) {
    for (int pivot = i*TILE_SIZE; pivot < (i+1)*TILE_SIZE; pivot++) {
        for (int row = group*TILE_SIZE; row < (group+1)*TILE_SIZE; row++) {
            float ratio = matrix[row * n + pivot];
            for (int col = pivot + 1; col < n; col++) {
                matrix[row * n + col] -= matrix[pivot * n + col] * ratio;
            }
            matrix[row * n + pivot] = 0;
        }
    }
}
