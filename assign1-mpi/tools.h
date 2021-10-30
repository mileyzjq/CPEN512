#include <math.h>

float* generateMatrix(int n) {
    float *arr = (float*)malloc((n*n)*sizeof(float));
    srand(time(0));
    for(int i=0; i<n*n; i++) {
        arr[i] = (float)(rand()%10)+1;
    }
    return arr;
}

void normaliseSingleRow(int start_point, int pivot, int n, float* matrix) {
    float ratio = matrix[pivot*n+start_point];
    for(int i=start_point+1; i<n; i++) {
        matrix[pivot*n+i] /= ratio;
    }
    matrix[pivot*n+start_point] = 1;
}

void serialGaussianElimination(int n, float* matrix) {
    float ratio;
    for(int pivot=0; pivot<n; pivot++) {
        // normalize first row
        normaliseSingleRow(pivot, pivot, n, matrix);
        //update all other rows
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
