/***************************************************************************
 *
 * MPI-version of row-wise Matrix-Matrix multiplication
 *
 *             File : matmul_mpi.c
 *        Author(s) : Håkan Grahn
 *          Created : 2009-01-30
 *    Last Modified : 2009-01-30
 * Last Modified by : Håkan Grahn
 *
 * © 2009 by Håkan Grahn, Blekinge Institute of Technology.
 * All Rights Reserved
 ***************************************************************************/

/*
 * Compile with:
 * mpicc -o mm matmul_mpi.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 1024   /* assumption: SIZE a multiple of number of nodes */
            /* Hint: use small sizes when testing, e.g., SIZE 8 */
#define FROM_MASTER 1   /* setting a message type */
#define FROM_WORKER 2   /* setting a message type */
#define DEBUG 0     /* 1 = debug on, 0 = debug off */

MPI_Status status;

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];
static double newb[SIZE][SIZE]; //store the column first b.

static void
init_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++) {
        /* Simple initialization, which enables us to easily check
         * the correct answer. Each element in c will have the same 
         * value as SIZE after the matmul operation.
         */
        a[i][j] = 1.0;
        b[i][j] = 1.0;
        }
}

static void
print_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++)
        printf(" %7.2f", c[i][j]);
    printf("\n");
    }
}

int
main(int argc, char **argv)
{
    int myrank, nproc;
    int rows, columns; /* amount of work per node (rows per worker) */
    int mtype; /* message type: send/recv between master and workers */
    int dest, src, offset;
    double start_time, end_time;
    int i, j, k;
    int numofRowBlock, numofColumnBlock;
    int sizeofRowBlock, sizeofColumnBlock;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0) {
    /* Master task */

    /* Initialization */
    printf("SIZE = %d, number of nodes = %d\n", SIZE, nproc);
    init_matrix();
    start_time = MPI_Wtime();

    /* Calculate the number of block of matrix in row and column*/
    switch(nproc){
        case 1:{
            numofRowBlock = 1;
            numofColumnBlock = 1;
            break;
        }
        case 2:{
            numofRowBlock = 1;
            numofColumnBlock = 2;
            break;
        }
        case 4:{
            numofRowBlock = 2;
            numofColumnBlock = 2;
            break;
        }
        case 8:{
            numofRowBlock = 4;
            numofColumnBlock = 2;
            break;
        }
        default:
        printf("the number of processes is wrong");
        break;
    }

    /* Calculate the size of each block in row and column*/
    sizeofRowBlock = SIZE / numofRowBlock ;
    sizeofColumnBlock = SIZE / numofColumnBlock;


    /* Transforming Matrix b into column first, store in newB. */
    for (i = 0; i < SIZE; i++){
            for (j = 0; j < SIZE; j++){
                newb[j][i] = b[i][j];
            }
        }

    /* Send a block of matrix and a block of matrix newb to workers */
    mtype = FROM_MASTER;
    rows = sizeofRowBlock;
    columns = sizeofColumnBlock;
    rowsOffset = rows;
    columnsOffset = columns;
    for (dest = 1; dest < nproc; dest++) {
        if (DEBUG)
        printf("sending %d rows and %d columns to task %d\n",rows,columns,dest);
        MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
        MPI_Send(&columns, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
        MPI_Send(&rowsOffset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
        MPI_Send(&columnsOffset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
        MPI_Send(&a[rowsOffset][columnsOffset], sizeofRowBlock*sizeofColumnBlock, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
        MPI_Send(&newb[columnsOffset][rowsOffset], sizeofRowBlock*sizeofColumnBlock, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
        rowsOffset += sizeofRowBlock;
        columnsOffset += sizeofColumnBlock;
    }

    /* let master do its part of the work */
    for (i = 0; i < sizeofRowBlock; i++) {
        for (j = 0; j < sizeofColumnBlock; j++) {
        c[i][j] = 0.0;
        for (k = 0; k < size; k++)
            c[i][j] = c[i][j] + a[i][k] * newb[j][k];
        }
    }

    /* collect the results from all the workers */
    mtype = FROM_WORKER;
    for (src = 1; src < nproc; src++) {
        MPI_Recv(&sizeofRowBlock, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&sizeofColumnBlock, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rowsOffset, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&columnsOffset, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&c[columnsOffset][rowsOffset], sizeofRowBlock*sizeofColumnBlock, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status);
        if (DEBUG)
        printf("recvd %d rows and %d columns from task %d, rowsOffset = %d, columnsOffset = %d\n",
               rows, columns, src, rowsOffset, columnsOffset);
    }

    end_time = MPI_Wtime();
    if (DEBUG)
        /* Prints the resulting matrix c */
        print_matrix();
    printf("Execution time on %2d nodes: %f\n", nproc, end_time-start_time);

    } else {
    /* Worker tasks */

    /* Receive data from master */
    mtype = FROM_MASTER;
    MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&rowsOffset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&columnsOffset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&a[rowsOffset][columnsOffset], sizeofRowBlock*sizeofColumnBlock, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&newb[columnsOffset][rowsOffset], sizeofRowBlock*sizeofColumnBlock, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
    if (DEBUG)
        printf ("Rank=%d, offset=%d, row =%d, a[offset][0]=%e, b[0][0]=%e\n",
            myrank, offset, rows, a[offset][0], b[0][0]);

    /* do the workers part of the calculation */
    for (i=offset; i<offset+rows; i++)
        for (j=0; j<SIZE; j++) {
        c[i][j] = 0.0;
        for (k=0; k<SIZE; k++)
            c[i][j] = c[i][j] + a[i][k] * newb[j][k];
        }
    if (DEBUG)
        printf ("Rank=%d, offset=%d, row =%d, c[offset][0]=%e\n",
            myrank, offset, rows, a[offset][0]);

    /* send the results to the master */
    mtype = FROM_WORKER;
    MPI_Send(&sizeofRowBlock, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
    MPI_Send(&sizeofColumnBlock, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
    MPI_Send(&rowsOffset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
    MPI_Send(&columnsOffset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
    MPI_Send(c[columnsOffset][rowsOffset], rows*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

