
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include "tools.h"

#define ROOT 0
#define MPI_SEND_TAG 5000
#define MPI_RECEIVE_TAG 5001

void gaussian(float *matrix, float *subMatrix);
void distributeRows(float *matrix, float *subMatrix);
void distributeRowsV2(float *matrix, float *subMatrix);
void gatherRows(float *matrix, float *subMatrix);
void receiveFromOne(int targetRank, float *row);
void sendToManyIsend(float *row, float *subMatrix, int sub_row, int main_row);
void awaitRequest(MPI_Request *request);
void sendToManyWithElimination(float *row, float *subMatrix, int sub_row, int m_row);

int rank;
int size;
int N;
int numRows;
MPI_Status status;
int debug = 1;

int main(int argc, char *argv[]){
    double start, stop, duration;
    N = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    numRows = N / size;
    float *matrix;
    float *subMatrix = (float*)malloc((N*numRows)*sizeof(float));

    if(rank == ROOT){
        matrix = generateMatrix(N);
        //matrix = getTestMatrix();
    }

    if(size == 1){
        start = MPI_Wtime();
        serialGaussianElimination(N, matrix);
        stop = MPI_Wtime();
        printf("1 processor: total time: %lf\n", (stop-start));
        MPI_Finalize();
        return 0;
    }

    distributeRows(matrix, subMatrix);
    //distributeRowsV2(matrix, subMatrix);

    if(rank == ROOT){
        start = MPI_Wtime();
    }

    gaussian(matrix, subMatrix);

    if(rank == ROOT){
        stop = MPI_Wtime();
        duration = stop - start;
        printf("%d processors cost total time: %lf\n", size, duration);
    }

    gatherRows(matrix, subMatrix);

    MPI_Finalize();

    if(rank == ROOT){
        //printMatrix(matrix, N, N);
        free(matrix);
    }
    free(subMatrix);

    return 0;
}

void distributeRows(float *matrix, float *subMatrix) {
    for(int i = 0; i < numRows; i++){
        MPI_Scatter(&matrix[i*N*size], N, MPI_FLOAT, &subMatrix[i*N], N, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
    }
}

void distributeRowsV2(float *matrix, float *subMatrix) {
    MPI_Request request;
    int index;
    float *row = (float*)malloc(N*sizeof(float));
    int tSize;
    MPI_Type_size(MPI_FLOAT, &tSize);
    for(int i = 0; i < numRows; i++){
        if(rank == ROOT) {
            for(int j=0; j<size; j++) {
                MPI_Isend(matrix+(i*size+j)*N, N, MPI_FLOAT, j, MPI_SEND_TAG, MPI_COMM_WORLD, &request);
            }
        }
        MPI_Recv(&subMatrix[i*N], N, MPI_FLOAT, ROOT, MPI_SEND_TAG, MPI_COMM_WORLD, &status);
    }
}

void gatherRows(float *matrix, float *subMatrix) {
    for(int i = 0; i < numRows; i++){
        MPI_Gather(&subMatrix[i*N], N, MPI_FLOAT, &matrix[i*size*N], N, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
    }
}

void gaussian(float *matrix, float *subMatrix) {
    float *row = (float*)malloc(N*sizeof(float));
    int localRow, localRank;
    float scale, ratio;

    for(int i = 0; i < N; i++){
        localRow = i / size;
        localRank = i % size;

        if(rank == localRank){
            scale = subMatrix[localRow*N+i];
            for(int j=i+1; j<N; j++){
                subMatrix[localRow*N+j] /= scale;
            }
            subMatrix[localRow*N+i] = 1;
            memcpy(row, &subMatrix[localRow*N], N*sizeof(float));

            //sendToManyWithElimination(row, subMatrix, localRow+1, i);
            sendToManyIsend(row, subMatrix, localRow+1, i);
        }else{
            //receiveFromOne(localRank, row);
            MPI_Recv(row, N, MPI_FLOAT, localRank, MPI_RECEIVE_TAG, MPI_COMM_WORLD, &status);

            for(int j=localRow; j<numRows; j++){
                if((j>localRow)|| (rank>localRank)){
                    ratio = subMatrix[j*N+i];
                    for(int k=i + 1; k<N; k++){
                        subMatrix[j*N+k] -= ratio * row[k];
                    }
                    subMatrix[j*N+i] = 0;
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(row);
}

void sendToManyWithElimination(float *row, float *subMatrix, int sub_row, int m_row) {
    float ratio;
    for(int i=0; i<size; i++) {
        if(rank != i) {
            MPI_Send(row, N, MPI_FLOAT, i, MPI_RECEIVE_TAG, MPI_COMM_WORLD);
        }
    }
    for(int j=sub_row; j<numRows; j++){
        ratio = subMatrix[j*N+m_row];
        for(int k=m_row+1; k<N; k++){
            subMatrix[j*N+k] -= ratio * row[k];
        }
        subMatrix[j*N+m_row] = 0;
    }
}

void sendToManyIsend(float *row, float *subMatrix, int sub_row, int main_row) {
    MPI_Request send_request;
    int s_row = sub_row;
    int n_updates = (numRows-s_row)/(size-1);
    float ratio;
    float *base = (float*)malloc(N*sizeof(float));
    memcpy(base, row, N * sizeof(float));
    for(int i=0; i<size; i++) {
        if(rank != i) {
            MPI_Isend(row, N, MPI_FLOAT, i, MPI_RECEIVE_TAG, MPI_COMM_WORLD, &send_request);
            for(int j=0; j<n_updates; j++) {
                ratio = subMatrix[s_row*N+main_row];
                for(int k = main_row+1; k<N; k++){
                    subMatrix[s_row*N+k] -= ratio*base[k];
                }
                subMatrix[s_row*N+main_row] = 0;
                s_row++;
            }
        }
    }
    while(s_row<numRows) {
        ratio = subMatrix[s_row*N+main_row];
        for(int k = main_row+1; k < N; k++){
            subMatrix[s_row*N+k] -= ratio * row[k];
        }
        subMatrix[s_row*N+main_row] = 0;
        s_row++;
    }
}

void awaitRequest(MPI_Request *request) {
    int count = 0;
    int flag = 0;
    while(flag == 0) {
        MPI_Test(request, &flag, &status);
        count++;
    }
}

void receiveFromOne(int targetRank, float *row) {
    MPI_Request recv_request;
    MPI_Irecv(row, N, MPI_FLOAT, targetRank, MPI_RECEIVE_TAG, MPI_COMM_WORLD, &recv_request);
    awaitRequest(&recv_request);
}
