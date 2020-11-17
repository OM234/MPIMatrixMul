/*
 * Osman Momoh
 * COMP 6331 : Distributed Systems
 * Assignment 4 : Parallel Matrix Multiplication using MPI
 * Fall Semester 2020
 */

#include<stdlib.h>
#include<time.h>
#include<mpi.h>

#define ROWSA 500       //Rows in left-sided matrix A
#define COLSA 32        //Columns in left-sided matrix A
#define ROWSB 32        //Rows in right-sided matrix B
#define COLSB 500       //Columns in right sided matrix B
#define ROWSC 500       //Rows in result matrix C
#define COLSC 500       //Columns in result matrix C
#define MASTER 0        //ID of master task
#define MASTERTAG 1     //Tag of master task
#define WORKERTAG 2     //Tag of worker task

void makeMatrices();                                //Make array A and B
void makeArrA(int arr[ROWSA][COLSA]);               //Make array A
void makeArrB(int arr[ROWSB][COLSB]);               //Make array B
void printArrA(int arr[ROWSA][COLSA]);              //Print array A
void printArrB(int arr[ROWSB][COLSB]);              //Print array B
void printArrC(int arr[ROWSC][COLSC], int);         //Print array C
void nonParallelExec();                             //When number of tasks is 1, perform multiplication sequentially
void initMPI(int &argc, char **&argv);              //Start the MPI
void masterTask();                                  //Perform the master task
void workerTask();                                  //Perform a worker task
void calculateRowDistributions();                   //Calculate how rows are distributed to worker tasks
void checkInput();                                  //Checks the definitions are valid

int matA[ROWSA][COLSA];                 //Matrix A
int matB[ROWSB][COLSB];                 //Matrix B
int matC[ROWSC][COLSC];                 //Matrix C
int result[ROWSC][COLSC];               //The final resulting matrix
int numTasks,                           //number of tasks
    taskID,                             //ID of task
    rowsPerTask = 0,                    //Rows per each task
    tasksWithExtraRow = 0,              //Number of tasks with extra row
    rowStart = 0,                       //Row in which a worker task should begin
    rowsToDo = 0;                       //Number of rows a task should process
clock_t start, end;                     //Start and end times
MPI_Status status;                      //Status of MPI receive

int main(int argc, char* argv[]) {

    checkInput();
    initMPI(argc, argv);
    makeMatrices();

    //If number of tasks is 1, perform sequential multiplication and terminate
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

    //Number of integers to send to worker tasks
    int MPICountA, MPICountB;
    //Arrays for holding which row each task started on and how many rows were processed
    int rowStartArr[numTasks-1];
    int rowsToDoArr[numTasks-1];

    //Send messages to the workers
    for(int i = 1 ; i < numTasks ; i++) {

        //Checks if this is a task which should have an extra row. It then multiplies the number of rows *
        //the number of columns of A, to get the number of integers to send in matrix A
        MPICountA = i <= tasksWithExtraRow ? (rowsPerTask + 1) * COLSA : rowsPerTask * COLSA;
        //The entire matrix B needs to be sent to each worker
        MPICountB = ROWSB * COLSB;
        //Which row the worker should start with in A
        rowStart += rowsToDo;
        //How many rows of A the worker should process
        rowsToDo = MPICountA / COLSA;
        //Store these values in an array
        rowStartArr[i-1] = rowStart;
        rowsToDoArr[i-1] = rowsToDo;

        //Send the number of rows in A to do in worker i
        MPI_Send(&rowsToDo, 1, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        //Which row worker i should start on in A
        MPI_Send(&rowStart, 1, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        //Send matrix A, beginning at the appropriate row, and with the right number of integers
        MPI_Send(&matA[rowStart], MPICountA, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
        //Send the entirety of matrix B
        MPI_Send(&matB, MPICountB, MPI_INT, i, MASTERTAG, MPI_COMM_WORLD);
    }

    //Receive messages from the workers
    for(int i = 1 ; i < numTasks ; i++) {

        //Receive the array results of the worker
        MPI_Recv(&matC[rowStartArr[i-1]], rowsToDoArr[i-1]*COLSC, MPI_INT, i, WORKERTAG, MPI_COMM_WORLD, &status);

        //Take results from worker arrays, add them to result array
        for(int y = rowStartArr[i-1]; y < rowStartArr[i-1] + rowsToDoArr[i-1]; y++) {
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

    //Receive from master number of rows to do
    MPI_Recv(&rowsToDo, 1, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    //Receive from master which row to start on
    MPI_Recv(&rowStart, 1, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    //Receive from master matrix A at appropriate start row
    MPI_Recv(&matA[rowStart], rowsToDo*COLSA, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);
    //Receive from master all of matrix B
    MPI_Recv(&matB, ROWSB*COLSB, MPI_INT, 0, MASTERTAG, MPI_COMM_WORLD, &status);

    //printf("rank %d starts at row: %d and has %d rows to do\n", taskID, rowStart, rowsToDo);

    //Perform the matrix multiplication of the appropriate rows
    for(int i = rowStart; i < rowStart + rowsToDo; i++) {
        for(int j = 0; j < COLSC; j++){

            matC[i][j] = 0;

            for(int k = 0; k < COLSA; k++){
                matC[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }

    //printArrC(matC, taskID);

    //return this work back to the master
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
    tasksWithExtraRow = (int)(rowsWithExtraDouble * (numTasks - 1));
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

void checkInput() {

    if((COLSA != ROWSB) || (ROWSA != ROWSC) || (COLSB != COLSC)) {
        printf("invalid input");
        exit(0);
    }
}