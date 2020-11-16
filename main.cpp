#include<stdlib.h>
#include<mpi.h>
#define ROWSA 10
#define COLSA 32
#define ROWSB 32
#define COLSB 10
#define ROWSC 10
#define COLSC 10

void makeArrA(int arr[ROWSA][COLSA]);
void makeArrB(int arr[ROWSB][COLSB]);
void printArrA(int arr[ROWSA][COLSA]);
void printArrB(int arr[ROWSB][COLSB]);
void printArrC(int arr[ROWSC][COLSC]);

int matA[ROWSA][COLSA];
int matB[ROWSB][COLSB];
int matC[ROWSC][COLSC];

int main(int argc, char* argv[]) {

    makeArrA(matA);
    makeArrB(matB);
    printArrA(matA);
    printArrB(matB);

    MPI_Init(&argc, &argv);
    MPI 
}

void makeArrA(int arr[ROWSA][COLSA]) {

    for(int i = 0; i < ROWSA; i++) {
        for(int j = 0; j < COLSA; j++){
            arr[i][j] = i+j;
        }
    }
}

void makeArrB(int arr[ROWSB][COLSB]) {

    for(int i = 0; i < ROWSB; i++) {
        for(int j = 0; j < COLSB; j++){
            arr[i][j] = i + 2*j;
        }
    }
}

void printArrA(int arr[ROWSA][COLSA]) {

    printf("Matrix A\n");

    for(int i = 0; i < ROWSA; i++) {
        for(int j = 0; j < COLSA; j++){
            printf("%d\t", arr[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void printArrB(int arr[ROWSB][COLSB]) {

    printf("Matrix B\n");

    for(int i = 0; i < ROWSB; i++) {
        for(int j = 0; j < COLSB; j++){
            printf("%d\t", arr[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void printArrC(int arr[ROWSC][COLSC]) {

    printf("Matrix C\n");

    for(int i = 0; i < ROWSC; i++) {
        for(int j = 0; j < COLSC; j++){
            printf("%d\t", arr[i][j]);
        }
        printf("\n");
    }

    printf("\n\n");
}
