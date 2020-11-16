#include<stdlib.h>
#include<time.h>
#include<mpi.h>

#define ROWSA 10
#define COLSA 10 //TODO: Change to 32
#define ROWSB 10 //TODO: Change to 32
#define COLSB 10
#define ROWSC 10
#define COLSC 10
#define MASTER 0
#define MASTERTAG 1
#define WORKERTAG 2

void makeArrA(int arr[ROWSA][COLSA]);
void makeArrB(int arr[ROWSB][COLSB]);
void printArrA(int arr[ROWSA][COLSA]);
void printArrB(int arr[ROWSB][COLSB]);
void printArrC(int arr[ROWSC][COLSC]);
void nonParallelExec();
void makeAndPrintMatrices();
void initMPI(int &argc, char **&argv);

void masterTask();

void workerTask();

int findRowsWithOneExtra();

int matA[ROWSA][COLSA];
int matB[ROWSB][COLSB];
int matC[ROWSC][COLSC];
int numTasks, taskID, rowsPerTask = 0, rowsWithExtra = 0, rowStart = 0, rowsToDo = 0;
MPI_Status status;

int main(int argc, char* argv[]) {

    makeAndPrintMatrices();
    initMPI(argc, argv);

    if(numTasks == 1) {
        nonParallelExec();
    }

    rowsPerTask = ROWSA / (numTasks-1);
    rowsWithExtra = findRowsWithOneExtra();

    if(taskID == MASTER) {
        masterTask();
    } else {
        workerTask();
    }

    MPI_Finalize();
}

void masterTask() {

    int MPICountA, MPICountB;

    for(int i = 1 ; i < numTasks ; i++) {

        //MPICountA = i != numTasks - 1 ? rowsPerTask * COLSA : (rowsPerTask+remainderRows) * COLSA;
        MPICountA = i <= rowsWithExtra ? (rowsPerTask + 1) * COLSA : rowsPerTask * COLSA;
        MPICountB = ROWSB * COLSB;
        rowStart = (i-1) * rowsPerTask;
        rowsToDo = MPICountA / COLSA;

        MPI_Send(&rowsToDo, 1, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&rowStart, 1, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&matA, MPICountA, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&matB, MPICountB, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
    }
}

void workerTask() {

    int numRows;
    int count = 0;

    MPI_Recv(&rowsToDo, 1, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
//    MPI_Recv(&matA, MPICountA,  MPI_INT, MASTER, MASTERTAG, MPI_COMM_WORLD, &status);

    printf("%d\n", rowsToDo);
//    for(int i = 0; i < ROWSC; i++) {
//        for(int j = 0; j < COLSC; j++){
//
//            matC[i][j] = 0;
//
//            for(int k = 0; k < COLSA; k++){
//                matC[i][j] += matA[i][k] * matB[k][j];
//            }
//        }
//    }
}

int findRowsWithOneExtra() {

    (int)(((ROWSA + 0.0) / (numTasks-1) - 1) * (numTasks-1));
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
    printf("Non-parallel execution done in %f us\n", execTime);

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