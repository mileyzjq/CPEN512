#include <errno.h>
#include <math.h>

#define TILE_SIZE 64

void setThreads(float* A);
void *gaussPthreads(void *args);
void serialGaussianEliminationGroup(int n, float* matrix, int group);
void serialGaussianEliminationBetween(int n, float *matrix, int i, int group);

void normaliseSingleRow(int startPoint, int pivot, int n, float* matrix) {
    float ratio = matrix[pivot*n+startPoint];
    for(int i=startPoint+1; i<n; i++) {
        matrix[pivot*n+i] /= ratio;
    }
    matrix[pivot*n+startPoint] = 1;
}

void eliminateSingleRow(int baseRow, int targetRow, int N, float *matrix) {
    float scale = matrix[targetRow*N+baseRow];
    for(int k=baseRow+1; k<N; k++) {
        matrix[targetRow*N+k] -= matrix[baseRow*N+k] * scale;
    }
    matrix[targetRow*N+baseRow] = 0;
}

void serialGaussianElimination(int n, float* matrix) {
    float ratio;
    for(int pivot=0; pivot<n; pivot++) {
        normaliseSingleRow(pivot, pivot, n, matrix);
        for(int row=pivot+1; row<n; row++) {
            float ratio = matrix[row*n+pivot];
            for(int col = pivot+1; col<n; col++) {
                matrix[row*n+col] -= matrix[pivot*n+col]*ratio;
            }
            matrix[row*n+pivot] = 0;
        }
    }
}

float* getTestMatrix() {
    int n = 4;
    float *matrix = (float*)malloc((n*n)*sizeof(float));
    float ori[16] = {2,4,6,8,2,3,1,2,3,3,2,1,1,3,2,5};
    for(int i=0; i<n*n; i++) {
        matrix[i] = ori[i];
    }
    return matrix;
}

float* generateMatrix(int n) {
    float *arr = (float*)malloc((n*n)*sizeof(float));
    srand(time(0));
    for(int i=0; i<n*n; i++) {
        arr[i] = (float)(rand()%100+1)/10;
    }
    return arr;
}

void printMatrix(float *matrix, int m, int n) {
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            printf("%f ", matrix[i*n+j]);
        }
        printf("\n");
    }
}

void validateMatrix(float *A, float *B, int n) {
    float epsilon = 0.0005;
    for(int i=0; i<n*n; i++) {
        if(fabsf(A[i]-B[i]) > epsilon) {
            perror("Validation Fails! \n");
            return;
        }
    }
    printf("Validation is successful!\n");
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
