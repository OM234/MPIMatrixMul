#include<stdlib.h>
#include<time.h>
#include<mpi.h>

#define ROWSA 5
#define COLSA 5 //TODO: Change to 32
#define ROWSB 5 //TODO: Change to 32
#define COLSB 5
#define ROWSC 5
#define COLSC 5
#define MASTER 0
#define MASTERTAG 1
#define WORKERTAG 2

void makeArrA(int arr[ROWSA][COLSA]);
void makeArrB(int arr[ROWSB][COLSB]);
void printArrA(int arr[ROWSA][COLSA]);
void printArrB(int arr[ROWSB][COLSB]);
void printArrC(int arr[ROWSC][COLSC], int);
void nonParallelExec();
void makeMatrices();
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

    makeMatrices();
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

    printArrA(matA);
    printArrB(matB);

    for(int i = 1 ; i < numTasks ; i++) {

        //MPICountA = i != numTasks - 1 ? rowsPerTask * COLSA : (rowsPerTask+remainderRows) * COLSA;
        MPICountA = i <= rowsWithExtra ? (rowsPerTask + 1) * COLSA : rowsPerTask * COLSA;
        MPICountB = ROWSB * COLSB;
        rowStart += rowsToDo;
        rowsToDo = MPICountA / COLSA;

        MPI_Send(&rowsToDo, 1, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&rowStart, 1, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&matA, MPICountA, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&matB, MPICountB, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
    }

    for(int i = 1 ; i < numTasks ; i++) {

        MPI_Recv(&rowStart, 1, MPI_INT, i, WORKERTAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&rowsToDo, 1, MPI_INT, i, WORKERTAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&matC, rowsToDo*COLSC, MPI_INT, i, WORKERTAG, MPI_COMM_WORLD, &status);

        //printf("Master %d received start at row: %d and has %d rows to do\n", taskID, rowStart, rowsToDo);

        printArrC(matC, taskID);
    }
}

void workerTask() {

    int numRows;
    int count = 0;

    MPI_Recv(&rowsToDo, 1, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&rowStart, 1, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&matA, rowsToDo*COLSA, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&matB, ROWSB*COLSB, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);

    printf("rank %d starts at row: %d and has %d rows to do\n", taskID, rowStart, rowsToDo);

    for(int i = rowStart; i < rowStart + rowsToDo; i++) {
        for(int j = 0; j < COLSC; j++){

            matC[i][j] = 0;

            for(int k = 0; k < COLSA; k++){
                matC[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }

    printArrC(matC, taskID);

    MPI_Send(&rowStart, 1, MPI_INT, 0, WORKERTAG, MPI_COMM_WORLD);
    MPI_Send(&rowsToDo, 1, MPI_INT, 0, WORKERTAG, MPI_COMM_WORLD);
    MPI_Send(&matC[rowStart], rowsToDo*COLSC, MPI_INT, 0, WORKERTAG, MPI_COMM_WORLD);
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
    printArrC(matC, taskID);
    printf("Non-parallel execution done in %f us\n", execTime);

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

int findRowsWithOneExtra() {

    double rowsWithExtra = (ROWSA + 0.0) / (numTasks-1);
    int intPart = (int) rowsWithExtra;
    rowsWithExtra -= intPart;
    return (int)(rowsWithExtra * (numTasks-1));
}

void initMPI(int &argc, char **&argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
}

void makeMatrices() {

    makeArrA(matA);
    makeArrB(matB);
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

void printArrC(int arr[ROWSC][COLSC], int theTaskID) {

    printf("Matrix C from task %d\n", theTaskID);

    for(int i = 0; i < ROWSC; i++) {
        for(int j = 0; j < COLSC; j++){
            printf("%6d", arr[i][j]);
        }
        printf("\n");
    }

    printf("\n\n");
}