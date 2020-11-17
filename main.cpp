#include<stdlib.h>
#include<time.h>
#include<mpi.h>

#define ROWSA 2000      //Rows in left-sided matrix A
#define COLSA 32        //Columns in left-sided matrix A
#define ROWSB 32        //Rows in right-sided matrix B
#define COLSB 2000      //Columns in right sided matrix B
#define ROWSC 2000      //Rows in result matrix C
#define COLSC 2000      //Columns in result matrix C
#define MASTER 0        //ID of master task
#define MASTERTAG 1     //Tag of master task
#define WORKERTAG 2     //Tag of worker task

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
void calculateRowDistributions();

int matA[ROWSA][COLSA];
int matB[ROWSB][COLSB];
int matC[ROWSC][COLSC];
int result[ROWSC][COLSC];
int numTasks, taskID, rowsPerTask = 0, rowsWithExtra = 0, rowStart = 0, rowsToDo = 0;
clock_t start, end;
MPI_Status status;

int main(int argc, char* argv[]) {

    initMPI(argc, argv);
    makeMatrices();

    if(numTasks == 1) {
        nonParallelExec();
    }

    calculateRowDistributions();

    if(taskID == MASTER) {
        masterTask();
    } else {
        workerTask();
    }

    MPI_Finalize();
}

void masterTask() {

    //printArrA(matA);
    //printArrB(matB);
    start = clock();

    int MPICountA, MPICountB;

    for(int i = 1 ; i < numTasks ; i++) {

        MPICountA = i <= rowsWithExtra ? (rowsPerTask + 1) * COLSA : rowsPerTask * COLSA;
        MPICountB = ROWSB * COLSB;
        rowStart += rowsToDo;
        rowsToDo = MPICountA / COLSA;

        MPI_Send(&rowsToDo, 1, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&rowStart, 1, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&matA[rowStart], MPICountA, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        MPI_Send(&matB, MPICountB, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
    }

    for(int i = 1 ; i < numTasks ; i++) {

        MPI_Recv(&rowStart, 1, MPI_INT, i, WORKERTAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&rowsToDo, 1, MPI_INT, i, WORKERTAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&matC[rowStart], rowsToDo*COLSC, MPI_INT, i, WORKERTAG, MPI_COMM_WORLD, &status);


        for(int y = rowStart; y < rowStart + rowsToDo; y++) {
            for(int x = 0; x < ROWSC; x++) {
                result[y][x] = matC[y][x];
            }
        }
        //printf("Master %d received start at row: %d and has %d rows to do\n", taskID, rowStart, rowsToDo);
    }

    end = clock();
    double execTime = (double)(end-start)/CLOCKS_PER_SEC * 1000000;
    //printArrC(result, taskID);
    printf("Parallel execution done in %f us\n", execTime);
}

void workerTask() {

    MPI_Recv(&rowsToDo, 1, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&rowStart, 1, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&matA[rowStart], rowsToDo*COLSA, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&matB, ROWSB*COLSB, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);

    //printf("rank %d starts at row: %d and has %d rows to do\n", taskID, rowStart, rowsToDo);

    for(int i = rowStart; i < rowStart + rowsToDo; i++) {
        for(int j = 0; j < COLSC; j++){

            matC[i][j] = 0;

            for(int k = 0; k < COLSA; k++){
                matC[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }

    //printArrC(matC, taskID);

    MPI_Send(&rowStart, 1, MPI_INT, 0, WORKERTAG, MPI_COMM_WORLD);
    MPI_Send(&rowsToDo, 1, MPI_INT, 0, WORKERTAG, MPI_COMM_WORLD);
    MPI_Send(&matC[rowStart], rowsToDo*COLSC, MPI_INT, 0, WORKERTAG, MPI_COMM_WORLD);
}


void nonParallelExec() {

    start = clock();

    for(int i = 0; i < ROWSC; i++) {
        for(int j = 0; j < COLSC; j++){

            matC[i][j] = 0;

            for(int k = 0; k < COLSA; k++){
                matC[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }

    end = clock();
    double execTime = (double)(end-start)/CLOCKS_PER_SEC * 1000000;
    //printArrC(matC, taskID);
    printf("Non-parallel execution done in %f us\n", execTime);

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

void calculateRowDistributions() {

    rowsPerTask = ROWSA / (numTasks-1);

    double rowsWithExtraDouble = (ROWSA + 0.0) / (numTasks - 1);
    int intPart = (int) rowsWithExtraDouble;
    rowsWithExtraDouble -= intPart;
    rowsWithExtra = (int)(rowsWithExtraDouble * (numTasks - 1));
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