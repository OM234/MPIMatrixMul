#include<stdlib.h>
#include<time.h>
#include<mpi.h>
#define ROWSA 4
#define COLSA 3 //TODO: Change to 32
#define ROWSB 3 //TODO: Change to 32
#define COLSB 4
#define ROWSC 4
#define COLSC 4

void makeArrA(int arr[ROWSA][COLSA]);
void makeArrB(int arr[ROWSB][COLSB]);
void printArrA(int arr[ROWSA][COLSA]);
void printArrB(int arr[ROWSB][COLSB]);
void printArrC(int arr[ROWSC][COLSC]);
void nonParallelExec();
void makeAndPrintMatrices();

void initMPI(int &argc, char **&argv);

int matA[ROWSA][COLSA];
int matB[ROWSB][COLSB];
int matC[ROWSC][COLSC];
int numTasks, taskID;;

int main(int argc, char* argv[]) {

    makeAndPrintMatrices();
    initMPI(argc, argv);

    if(numTasks == 1) {
        nonParallelExec();
    }


    MPI_Finalize();
}

void initMPI(int &argc, char **&argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
}

void makeAndPrintMatrices() {

    makeArrA(matA);
    makeArrB(matB);
    printArrA(matA);
    printArrB(matB);
}

void nonParallelExec() {

    clock_t begin = clock();

    for(int i = 0; i < ROWSC; i++) {
        for(int j = 0; j < COLSC; j++){

            matC[i][j] = 0;

            for(int k = 0; k < COLSA; k++){
                matC[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }

    clock_t end = clock();
    double execTime = (double)(end-begin)/CLOCKS_PER_SEC * 1000000;
    printArrC(matC);
    printf("Non-parellel execution done in %f us\n", execTime);

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

void makeArrA(int arr[ROWSA][COLSA]) {

    for(int i = 0; i < ROWSA; i++) {
        for(int j = 0; j < COLSA; j++){
            arr[i][j] = i;
        }
    }
}

void makeArrB(int arr[ROWSB][COLSB]) {

    for(int i = 0; i < ROWSB; i++) {
        for(int j = 0; j < COLSB; j++){
            arr[i][j] = j;
        }
    }
}

void printArrA(int arr[ROWSA][COLSA]) {

    printf("Matrix A\n");

    for(int i = 0; i < ROWSA; i++) {
        for(int j = 0; j < COLSA; j++){
            printf("%3d", arr[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void printArrB(int arr[ROWSB][COLSB]) {

    printf("Matrix B\n");

    for(int i = 0; i < ROWSB; i++) {
        for(int j = 0; j < COLSB; j++){
            printf("%3d", arr[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void printArrC(int arr[ROWSC][COLSC]) {

    printf("Matrix C\n");

    for(int i = 0; i < ROWSC; i++) {
        for(int j = 0; j < COLSC; j++){
            printf("%6d", arr[i][j]);
        }
        printf("\n");
    }

    printf("\n\n");
}