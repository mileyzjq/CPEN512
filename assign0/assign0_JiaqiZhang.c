/
// Created by Jiaqi Zhang on 2021/9/19.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TILE_SIZE 512
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#if defined(SET_INT)
typedef int new_t;
#elif defined(SET_FLOAT)
typedef float new_t;
#elif defined(SET_DOUBLE)
typedef double new_t;
#endif

// declare functions
new_t** initializeMatrix(int n);
new_t** generateSquareMatrix(int n);
void freeMatrix(int n, new_t** arr);
void setZero(int n, new_t **arr);
new_t** transformMatrix(int n, new_t **arr);
void multiplyVanillaMatrix(int n, new_t** arr1, new_t** arr2, new_t** res);
void transformOptimization(int n, new_t** arr1, new_t** arr2, new_t** res);
void loopOptimization(int n, new_t** arr1, new_t** arr2, new_t** res);
void tileOptimization(int n, new_t** arr1, new_t** arr2, new_t** res);
void combinedOptimization(int n, new_t** arr1, new_t** arr2, new_t** res);
void matrixAddtion(int n, new_t** arr1, new_t** arr2, new_t** res);
void matrixSubstraction(int n, new_t** arr1, new_t** arr2, new_t** res);
new_t** generateSubMatrix(int n, new_t** matrix, int x, int y);

int main(int argc, char *argv[]) {
    clock_t start, stop;
    double duration;
    //int type = atoi(argv[1]);
    // n is matrix size
    int n = atoi(argv[1]);
    
    new_t **matrix1 = generateSquareMatrix(n);
    new_t **matrix2 = generateSquareMatrix(n);
    new_t **res = initializeMatrix(n);

    start = clock();
    multiplyVanillaMatrix(n, matrix1, matrix2, res);
    stop = clock();
    duration = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("Vanilla Matrix time is %lf seconds\n", duration);
    setZero(n, res);
    
    //transform Optimization
    start = clock();
    transformOptimization(n, matrix1, matrix2, res);
    stop = clock();
    duration = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("Transformation optimization Square Matrix time is %lf seconds\n", duration);
    setZero(n, res);

    // loop Optimization
    start = clock();
    loopOptimization(n, matrix1, matrix2, res);
    stop = clock();
    duration = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("Loop optimization Square Matrix time is %lf seconds\n", duration);
    setZero(n, res);

    // tile optimization
    start = clock();
    tileOptimization(n, matrix1, matrix2, res);
    stop = clock();
    duration = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("Tile optimization Square Matrix time is %lf seconds\n", duration);
    setZero(n, res);

    // Combine Optimization: combine Strassen and Tile algorithm
    start = clock();
    combinedOptimization(n, matrix1, matrix2, res);
    stop = clock();
    duration = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("Combined optimization Square Matrix time is %lf seconds\n", duration);
    setZero(n, res);

    freeMatrix(n, matrix1);
    freeMatrix(n, matrix2);
    freeMatrix(n, res);

    return 0;
}

// generate A[n][n] matrix with random numbers between 0 and 10
new_t** generateSquareMatrix(int n) {
    int i, j;
    new_t** arr;
    srand(time(0));
    arr = (new_t**)malloc(n*sizeof(new_t*));
    for(i=0; i<n; i++) {
        arr[i] = (new_t*)malloc(n*sizeof(new_t));
    }
    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            #if defined(SET_INT)
            arr[i][j] = rand()%10;
            #else
            arr[i][j] = (new_t)rand()/RAND_MAX;
            #endif
        }
    }

    return arr;
}

// initialize n*n matrix
new_t** initializeMatrix(int n) {
    new_t** res;
    res = (new_t**)malloc(n*sizeof(new_t*));
    for(int i=0; i<n; i++) {
        res[i] = (new_t*)malloc(n*sizeof(new_t));
    }
    return res;
}

// free malloc
void freeMatrix(int n, new_t** arr) {
    for(int i=0; i<n; i++) {
        free(arr[i]);
    }
    free(arr);
}

void setZero(int n, new_t **arr) {
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            arr[i][j] = (new_t)0.0;
        }
    }
}

// multiply Vanilla matrix arr1 and arr2
void multiplyVanillaMatrix(int n, new_t** arr1, new_t** arr2, new_t** res) {
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            res[i][j] = 0;
            for(int k=0; k<n; k++) {
                res[i][j] += arr1[i][k] * arr2[k][j];
            }
        }
    }
}

/*** Manual Optimization 1: Matrix transformation ***/
void transformOptimization(int n, new_t** arr1, new_t** arr2, new_t** res) {
    new_t** transformed = transformMatrix(n, arr2);
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            for(int k=0; k<n; k++) {
                res[i][j] += arr1[i][k]*transformed[j][k];
            }
        }
    }
}

new_t** transformMatrix(int n, new_t **arr) {
    new_t** transformed = initializeMatrix(n);
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            transformed[i][j] = arr[j][i];
        }
    }
    return transformed;
}

/*** Manual Optimization 2: loop Optimization, change the order of loops ***/
void loopOptimization(int n, new_t** arr1, new_t** arr2, new_t** res) {
    for(int i=0; i<n; i++) {
        for(int k=0; k<n; k++) {
            for(int j=0; j<n; j++) {
                res[i][j] += arr1[i][k]*arr2[k][j];
            }
        }
    }
}

/*** Manual Optimization 3: Tile Matrix multiplication ***/
void tileOptimization(int n, new_t** arr1, new_t** arr2, new_t** res) {
    int minI2, minI3;
    for (int II2 = 0; II2 < n; II2 += TILE_SIZE) {
        for(int II3 = 0; II3<n; II3 += TILE_SIZE) {
            minI2 = min(II2+TILE_SIZE, n);
            minI3 = min(II3+TILE_SIZE, n);
            for(int I1 = 0; I1<n; I1++) {
                for(int I2 = II2; I2 < minI2; I2++) {
                    for(int I3 = II3; I3 < minI3; I3++) {
                        res[I1][I3] += arr1[I1][I2] * arr2[I2][I3];
                    }
                }
            }
        }
    }
}


/*** Manual Optimization 3: Combined algorithm: Tile Optimization plus Strassen Optimization ***/
void combinedOptimization(int n, new_t** arr1, new_t** arr2, new_t** res) {
    int half = n/2;
    // get sub matrix
    new_t** a11 = generateSubMatrix(half, arr1, 0, 0);
    new_t** a12 = generateSubMatrix(half, arr1, 0, half);
    new_t** a21 = generateSubMatrix(half, arr1, half, 0);
    new_t** a22 = generateSubMatrix(half, arr1, half, half);
    new_t** b11 = generateSubMatrix(half, arr2, 0, 0);
    new_t** b12 = generateSubMatrix(half, arr2, 0, half);
    new_t** b21 = generateSubMatrix(half, arr2, half, 0);
    new_t** b22 = generateSubMatrix(half, arr2, half, half);
    // addtion and substraction part
    new_t** s1 = initializeMatrix(half);
    new_t** s2 = initializeMatrix(half);
    new_t** s3 = initializeMatrix(half);
    new_t** s4 = initializeMatrix(half);
    new_t** s5 = initializeMatrix(half);
    new_t** s6 = initializeMatrix(half);
    new_t** s7 = initializeMatrix(half);
    new_t** s8 = initializeMatrix(half);
    new_t** s9 = initializeMatrix(half);
    new_t** s10 = initializeMatrix(half);
    matrixSubstraction(half, b12, b22, s1);
    matrixAddtion(half, a11, a12, s2);
    matrixAddtion(half, a21, a22, s3);
    matrixSubstraction(half, b21, b11, s4);
    matrixAddtion(half, a11, a22, s5);
    matrixAddtion(half, b11, b22, s6);
    matrixSubstraction(half, a12, a22, s7);
    matrixAddtion(half, b21, b22, s8);
    matrixSubstraction(half, a11, a21, s9);
    matrixAddtion(half, b11, b12, s10);

    // multiplication part
    new_t** p1 = initializeMatrix(half);
    new_t** p2 = initializeMatrix(half);
    new_t** p3 = initializeMatrix(half);
    new_t** p4 = initializeMatrix(half);
    new_t** p5 = initializeMatrix(half);
    new_t** p6 = initializeMatrix(half);
    new_t** p7 = initializeMatrix(half);
    tileOptimization(half, a11, s1, p1);
    tileOptimization(half, s2, b22, p2);
    tileOptimization(half, s3, b11, p3);
    tileOptimization(half, a22, s4, p4);
    tileOptimization(half, s5, s6, p5);
    tileOptimization(half, s7, s8, p6);
    tileOptimization(half, s9, s10, p7);

    new_t** c11 = initializeMatrix(half);
    new_t** c22 = initializeMatrix(half);
    new_t** c12 = initializeMatrix(half);
    new_t** c21 = initializeMatrix(half);
    matrixAddtion(half, p5, p4, c11);
    matrixSubstraction(half, c11, p2, c11);
    matrixAddtion(half, c11, p6, c11);
    matrixAddtion(half, p1, p2, c12);
    matrixAddtion(half, p3, p4, c21);
    matrixAddtion(half, p5, p1, c22);
    matrixAddtion(half, p5, p1, c22);
    matrixSubstraction(half, c22, p3, c22);
    matrixSubstraction(half, c22, p7, c22);
    for(int i=0; i<half; i++) {
        for(int j=0; j<half; j++) {
            res[i][j] = c11[i][j];
            res[i][j+half] = c12[i][j];
            res[i+half][j] = c21[i][j];
            res[i+half][j+half] = c22[i][j];
        }
    }
    // free matrix
    freeMatrix(half, s1);
    freeMatrix(half, s2);
    freeMatrix(half, s3);
    freeMatrix(half, s4);
    freeMatrix(half, s5);
    freeMatrix(half, s6);
    freeMatrix(half, s7);
    freeMatrix(half, s8);
    freeMatrix(half, s9);
    freeMatrix(half, s10);
    freeMatrix(half, p1);
    freeMatrix(half, p2);
    freeMatrix(half, p3);
    freeMatrix(half, p4);
    freeMatrix(half, p5);
    freeMatrix(half, p6);
    freeMatrix(half, p7);
    freeMatrix(half, a11);
    freeMatrix(half, a12);
    freeMatrix(half, a21);
    freeMatrix(half, a22);
    freeMatrix(half, b11);
    freeMatrix(half, b12);
    freeMatrix(half, b21);
    freeMatrix(half, b22);
}

// generate size n Sub Matrix of Original Matrix from start poition(x,y)
new_t** generateSubMatrix(int n, new_t** matrix, int x, int y) {
    new_t** res = initializeMatrix(n);
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            res[i][j] = matrix[x+i][y+j];
        }
    }
    return res;
}

void matrixAddtion(int n, new_t** arr1, new_t** arr2, new_t** res) {
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            res[i][j] = arr1[i][j] + arr2[i][j];
        }
    }
}

void matrixSubstraction(int n, new_t** arr1, new_t** arr2, new_t** res) {
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            res[i][j] = arr1[i][j] - arr2[i][j];
        }
    }
}
