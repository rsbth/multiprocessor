/*****************************************************
 *
 * Gaussian elimination
 *
 * sequential version
 *
 *****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define MAX_SIZE 4096

typedef double matrix[MAX_SIZE][MAX_SIZE];

int	N;		/* matrix size		*/
int	maxnum;		/* max number of element*/
char	*Init;		/* matrix init type	*/
int	PRINT;		/* print switch		*/
matrix	A;		/* matrix A		*/
double	b[MAX_SIZE];	/* vector b             */
double	y[MAX_SIZE];	/* vector y             */

int thread_count; /* Global variable:  accessible to all threads */

/* forward declarations */
void* work(void* rank);
void Init_Matrix(void);
void Print_Matrix(void);
void Init_Default(void);
int Read_Options(int, char **);
pthread_barrier_t barrier;


int
main(int argc, char **argv)
{
	clock_t begin, end;
  	double time_spent;

  	begin = clock();
  	/* here, do your time-consuming job */

	long thread; /* Use long in case of a 64-bit system */
	pthread_t* thread_handles;

    thread_count = strtol(argv[1], NULL, 10); /*Get number of threads from command line */

	thread_handles = malloc (thread_count* sizeof(pthread_t));
	/*allocate storage for one pthread_t object for each thread.
	The pthread_t data structure is used for storing thread-speciﬁc information*/

    Init_Default();		/* Init default values	*/
    Read_Options(argc,argv);	/* Read arguments	*/
    Init_Matrix();		/* Init the matrix	*/
    pthread_barrier_init(&barrier, NULL, thread_count);

    for (thread = 0; thread < thread_count; thread++)
    	pthread_create(&thread_handles[thread], NULL, work, (void*) thread);
    //work();
    //call pthread_barrier_wait(phase1_barrier);
    
    if (PRINT == 1)
	Print_Matrix();

	for (thread = 0; thread < thread_count; thread++)
		pthread_join(thread_handles[thread], NULL);

	pthread_barrier_destroy(&barrier);
	free(thread_handles);

	end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nthe execution time of the Pthread version of the application is %12.4e seconds\n", time_spent);
	return 0;
}

void* work(void* rank)
{
	long my_rank = (long) rank;
    int i, j, k;

    /* Gaussian elimination algorithm, Algo 8.4 from Grama */
    /*
    for (k = 0; k < N; k++) { // Outer loop
		for (j = k+1; j < N; j++)
	    	A[k][j] = A[k][j] / A[k][k]; // Division step

	y[k] = b[k] / A[k][k];
	A[k][k] = 1.0;

	for (i = k+1; i < N; i++) {
	    for (j = k+1; j < N; j++)
			A[i][j] = A[i][j] - A[i][k]*A[k][j]; // Elimination step

	    b[i] = b[i] - A[i][k]*y[k];
	    A[i][k] = 0.0;
	}
    }*/
     for (int k = 0; k < N ;k++) /* Outer loop, second pass */
    {
        y[k] = b[k]/A[k][k];
        A[k][k] = 1.0;
        for (int i = k + 1 ; i < N ;i++)
        {
            b[i] = b[i] − A[i][k] * y[k];
            A[i][k] = 0;
        }
    }

    /* rearrange the algorithm to make this easier, so that there's only one loop over j*/
    for (k = 0; k < N; k++ ) /* Outer loop */
    {
    	// call pthread_barrier_wait(row_barrier);
    	pthread_barrier_wait(&barrier);
    	// for j := k + 1 + thread_number to n − 1 step n_threads do
        //for (j = k + 1; j< n;j++) /*this is the loop that can be done in parallel*/
        for (j = k+1+(my_rank); j < N; j += thread_count)
        {
            A[k][j] = A[k][j]/A[k][k]; /* Division step */
            for (i = k + 1; i < N ; i++ )
                A[i][j] = A[i][j] - A[i][k]*A[k][j];/* Elimination step */
        }
    }

}

void
Init_Matrix()
{
    int i, j;
 
    printf("\nsize      = %dx%d ", N, N);
    printf("\nmaxnum    = %d \n", maxnum);
    printf("Init	  = %s \n", Init);
    printf("Initializing matrix...");
 
    if (strcmp(Init,"rand") == 0) {
	for (i = 0; i < N; i++){
	    for (j = 0; j < N; j++) {
		if (i == j) /* diagonal dominance */
		    A[i][j] = (double)(rand() % maxnum) + 5.0;
		else
		    A[i][j] = (double)(rand() % maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(Init,"fast") == 0) {
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		if (i == j) /* diagonal dominance */
		    A[i][j] = 5.0;
		else
		    A[i][j] = 2.0;
	    }
	}
    }

    /* Initialize vectors b and y */
    for (i = 0; i < N; i++) {
	b[i] = 2.0;
	y[i] = 1.0;
    }

    printf("done \n\n");
    if (PRINT == 1)
	Print_Matrix();
}

void
Print_Matrix()
{
    int i, j;
 
    printf("Matrix A:\n");
    for (i = 0; i < N; i++) {
	printf("[");
	for (j = 0; j < N; j++)
	    printf(" %5.2f,", A[i][j]);
	printf("]\n");
    }
    printf("Vector b:\n[");
    for (j = 0; j < N; j++)
	printf(" %5.2f,", b[j]);
    printf("]\n");
    printf("Vector y:\n[");
    for (j = 0; j < N; j++)
	printf(" %5.2f,", y[j]);
    printf("]\n");
    printf("\n\n");
}

void 
Init_Default()
{
    N = 2048;
    Init = "rand";
    maxnum = 15.0;
    PRINT = 0;
}
 
int
Read_Options(int argc, char **argv)
{
    char    *prog;
 
    prog = *argv;
    while (++argv, --argc > 0)
	if (**argv == '-')
	    switch ( *++*argv ) {
	    case 'n':
		--argc;
		N = atoi(*++argv);
		break;
	    case 'h':
		printf("\nHELP: try sor -u \n\n");
		exit(0);
		break;
	    case 'u':
		printf("\nUsage: sor [-n problemsize]\n");
		printf("           [-D] show default values \n");
		printf("           [-h] help \n");
		printf("           [-I init_type] fast/rand \n");
		printf("           [-m maxnum] max random no \n");
		printf("           [-P print_switch] 0/1 \n");
		exit(0);
		break;
	    case 'D':
		printf("\nDefault:  n         = %d ", N);
		printf("\n          Init      = rand" );
		printf("\n          maxnum    = 5 ");
		printf("\n          P         = 0 \n\n");
		exit(0);
		break;
	    case 'I':
		--argc;
		Init = *++argv;
		break;
	    case 'm':
		--argc;
		maxnum = atoi(*++argv);
		break;
	    case 'P':
		--argc;
		PRINT = atoi(*++argv);
		break;
	    default:
		printf("%s: ignored option: -%s\n", prog, *argv);
		printf("HELP: try %s -u \n\n", prog);
		break;
	    } 
}
