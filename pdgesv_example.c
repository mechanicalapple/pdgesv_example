//Artem Anikeev
//BSD license
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int blacs_exitcode = 0;
//int DLEN_ = 9;//Array descriptor length in ScaLAPACK
int N, NB, M, MB, NPROW, NPCOL, ICTXT, MYROW, MYCOL, i, j, k, l, o, p, DESCA[9], DESCB[9], CSRC, NRHS, RSRC, NBRHS, MXLOCC, MXLOCR, MXRHSC, MXLLDA, MXLLDB, nblockc, nblockr, nbigc, nbigr, nbpc, nbpr, ntailc, ntailr, INFO, LDA, IA, IB, JA, JB, *IPIV;
double *A, *B, *SUB_A, *SUB_B;

//Declaration of Fortran subroutines
void sl_init_(int *, int *, int *);

void blacs_gridinfo_(int *, int *, int *, int *, int *);

void blacs_exit_(int *);

void blacs_gridexit_(int *);

void descinit_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

void dgesd2d_(int *, int *, int *, double *, int *, int *, int *);

void dgerv2d_(int *, int *, int *, double *, int *, int *, int *);

void pdgesv_(int *, int *, double *, int *, int *, int*, int*, double*, int *, int *, int*, int *);
//

//Declaration of Fortran functions
//

int main(int argc, char * argv[])
{
	//Get command line arguments
	if(argc != 5)
	{
		exit(1);
	}

	N = strtol(argv[1], NULL, 10);//Number of columns in the global matrix A
	M = N;//Number of rows in the global matrix A
	NB = strtol(argv[2], NULL, 10);//Column block size for the matrix A
	MB = NB;//Row block size for the matrix A
	NPROW = strtol (argv[3], NULL, 10);//Number of rows in the process grid
	NPCOL = strtol (argv[4], NULL, 10);//Number of columns in the process grid

	//printf("N = %i\nNPROW = %i\nNPCOL = %i\n", N, NPROW, NPCOL);

	//Initialize the process grid
	sl_init_(&ICTXT, &NPROW, &NPCOL);
	blacs_gridinfo_(&ICTXT, &NPROW, &NPCOL, &MYROW, &MYCOL);

	//If I'm not in the process grid, go to the end of the program
	if (MYROW == -1)
	{
		blacs_exit_(&blacs_exitcode);
		exit(0);
	}

	//Generate matrices A and B
	if ((MYROW == 0)&&(MYCOL == 0))
	{
		A = (double*)malloc(sizeof(double) * M * N);
		B = (double*)malloc(sizeof(double) * N);

		//Do not ever do like this in real world. Use high performance random array generators
		srand(time(NULL));

		printf("Matrix A:\n");
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < M; j++)
			{
				//ROW-MAJOR matrix
				A[i * N + j] = (double)rand() / RAND_MAX;
				printf("%e ", A[i * N + j]);
			}
			printf("\n");
		}

		printf("Matrix B:\n");
		for (i = 0; i < N; i++)
		{
			B[i] = (double)rand() / RAND_MAX;
			printf("%e ", B[i]);
		}
		printf("\n");
	}

	//*****
	//Distribute matrices to the process grid
	//*****

	LDA = 1;//The distance between two elements in matrix row
	IA = 1;
	JA = 1;
	IB = 1;
	JB = 1;

	CSRC = 0;//Process column over which the first column of the matrix is distributed
	NRHS = 1;//Number of columns in the global solution matrix B
	RSRC = 0;//Process row over which the first row of the matrix is distributed
	NBRHS = 1;//Column block size for the global solution matrix B
	MXRHSC = 1;//Maximum number of columns of the matrix B owned by any process column

	ntailc = N % NB;//Length of column tail in last block
	ntailr = M % MB;//Length of row tail in last block
	if (ntailc == 0)
	{
		nblockc = N / NB;//Number of column blocks in matrix A
	}
	else
	{
		nblockc = N / NB + 1;
	}
	if (ntailr == 0)
	{
		nblockr = M / MB;
	}
	else
	{
		nblockr = M / MB + 1;
	}

	nbigc = nblockc % NPCOL;//Number of big process columns
	nbpc = nblockc / NPCOL;//Number of blocks per small process column
	if (nbigc == 0)
	{
		MXLOCC = NB * nbpc;//Maximum number of columns of the matrix A owned by any process column
	}
	else
	{
		if (nbigc == 1)
		{
			MXLOCC = NB * nbpc + N % NB;
		}
		else
		{
			MXLOCC = NB * (nbpc + 1);
		}
	}

	nbigr = nblockr % NPROW;//Number of big process rows
	nbpr = nblockr / NPROW;//Number of blocks per small process row
	if (nbigr == 0)
	{
		MXLOCR = MB * nbpr;//Maximum number of rows of the matrix A owned by any process row
	}
	else
	{
		if (nbigr == 1)
		{
			MXLOCR = MB * nbpr + N % NB;
		}
		else
		{
			MXLOCR = MB * (nbpr + 1);
		}
	}

	if (MXLOCC > MXLOCR)
	{
		MXLLDA = MXLOCC;//Maximum local leading dimension of the array SUB_A
	}
	else
	{
		MXLLDA = MXLOCR;
	}

	MXLLDB = MXLOCR;//Maximum local leading dimension of the array SUB_B

	//Initialize array descriptors
	descinit_(DESCA, &M, &N, &MB, &NB, &RSRC, &CSRC, &ICTXT, &MXLLDA, &INFO);
	descinit_(DESCB, &N, &NRHS, &NB, &NBRHS, &RSRC, &CSRC, &ICTXT, &MXLLDB, &INFO);
	IPIV = (int*)malloc(sizeof(int) * (MXLOCR + NB));

	//printf("%i %i %i %i %i %i %i %i %i\n", DESCA[0], DESCA[1], DESCA[2], DESCA[3], DESCA[4], DESCA[5], DESCA[6], DESCA[7], DESCA[8]);
	//printf("%i %i %i %i %i %i %i %i %i\n", DESCB[0], DESCB[1], DESCB[2], DESCB[3], DESCB[4], DESCB[5], DESCB[6], DESCB[7], DESCB[8]);

	SUB_A = (double*)malloc(sizeof(double) * MXLLDA * MXLOCC);//COLUMN-MAJOR matrices due to FORTRAN 
	SUB_B = (double*)malloc(sizeof(double) * MXLLDB * MXRHSC);

	//Copy values from matrix to submatrix and send to appropriate process
	//Algorith is not optimal. It is designed to show cases. Use advanced algorithm from original ScaLAPACK documentation in real world
	//Matrix A
	//Reverse order to skip (0,0) in local memory
	if ((MYROW == 0)&&(MYCOL == 0))
	{
		//printf("ntailc = %i\n", ntailc);
		//printf("ntailr = %i\n", ntailr);
		//printf("nblockc = %i\n", nblockc);
		//printf("nblockr = %i\n", nblockr);
		//printf("nbigc = %i\n", nbigc);
		//printf("nbigr = %i\n", nbigr);
		//printf("nbpc = %i\n", nbpc);
		//printf("nbpr = %i\n", nbpr);
		//printf("MXLLDA = %i\n", MXLLDA);
		//printf("MXLLDB = %i\n", MXLLDB);
		//printf("MXLOCC = %i\n", MXLOCC);
		//printf("MXLOCR = %i\n", MXLOCR);
		//printf("MXRHSC = %i\n", MXRHSC);

		for (i = NPROW - 1; i >= nbigr; i--)
		{
			for (j = NPCOL - 1; j >= nbigc; j--)
			{
				//Small process row, small process column
				for (k = 0; k < nbpr; k++)
				{
					for (l = 0; l < nbpc; l++)
					{
						for (o = 0; o < MB; o++)
						{
							for (p = 0; p < NB; p++)
							{
								SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
							}
						}
					}
				}
				//Send submatrix
				if ((i != 0)||(j != 0))
				{
					dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
				}
			}

			if (nbigc > 0)
			{
				j = nbigc - 1;
				{
					//Small process row, tail process column
					for (k = 0; k < nbpr; k++)
					{
						for (l = 0; l < nbpc; l++)
						{
							for (o = 0; o < MB; o++)
							{
								for (p = 0; p < NB; p++)
								{
									SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];

								}
							}
						}
						l = nbpc;
						for (o = 0; o < MB; o++)
						{
							for (p = 0; p < ntailc; p++)
							{
								SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
							}
						}
					}
					//Send submatrix
					if ((i != 0)||(j != 0))
					{
						dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
					}
				}

				for (j = nbigc - 2; j >= 0; j--)
				{
					//Small process row, big process colummn
					for (k = 0; k < nbpr; k++)
					{
						for (l = 0; l <= nbpc; l++)
						{
							for (o = 0; o < MB; o++)
							{
								for (p = 0; p < NB; p++)
								{
									SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
								}
							}
						}
					}
					//Send submatrix
					if ((i != 0)||(j != 0))
					{
						dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
					}
				}
			}
		}

		if (nbigr > 0)
		{
			i = nbigr - 1;
			{
				for (j = NPCOL - 1; j >= nbigc; j--)
				{
					//Tail process row, small process column
					for (k = 0; k < nbpr; k++)
					{
						for (l = 0; l < nbpc; l++)
						{
							for (o = 0; o < MB; o++)
							{
								for (p = 0; p < NB; p++)
								{
									SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
								}
							}
						}
					}
					k = nbpr;
					for (l = 0; l < nbpc; l++)
					{
						for (o = 0; o < ntailr; o++)
						{
							for (p = 0; p < NB; p++)
							{
								SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
							}
						}
					}
					//Send submatrix
					if ((i != 0)||(j != 0))
					{
						dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
					}
				}

				if (nbigc > 0)
				{
					j = nbigc - 1;
					{
						//Tail process row, tail process column
						for (k = 0; k < nbpr; k++)
						{
							for (l = 0; l < nbpc; l++)
							{
								for (o = 0; o < MB; o++)
								{
									for (p = 0; p < NB; p++)
									{
										SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
									}
								}
							}
							l = nbpc;
							for (o = 0; o < MB; o++)
							{
								for (p = 0; p < ntailc; p++)
								{
									SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
								}
							}
						}
						k = nbpr;
						for (l = 0; l < nbpc; l++)
						{
							for (o = 0; o < ntailr; o++)
							{
								for (p = 0; p < NB; p++)
								{
									SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
								}
							}
						}
						l = nbpc;
						for (o = 0; o < ntailr; o++)
						{
							for (p = 0; p < ntailc; p++)
							{
								SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
							}
						}
						//Send submatrix
						if ((i != 0)||(j != 0))
						{
							dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
						}
					}

					for (j = nbigc - 2; j >= 0; j--)
					{
						//Tail process row, big process column
						for (k = 0; k < nbpr; k++)
						{
							for (l = 0; l <= nbpc; l++)
							{
								for (o = 0; o < MB; o++)
								{
									for (p = 0; p < NB; p++)
									{
										SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
									}
								}
							}
						}
						k = nbpr;
						for (l = 0; l <= nbpc; l++)
						{
							for (o = 0; o < ntailr; o++)
							{
								for (p = 0; p < NB; p++)
								{
									SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
								}
							}
						}
						//Send submatrix
						if ((i != 0)||(j != 0))
						{
							dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
						}
					}
				}
			}

			for (i = nbigr - 2 ; i >= 0; i--)
			{
				for (j = NPCOL - 1; j >= nbpc; j--)
				{
					//Big process row, small process column
					for (k = 0; k <= nbpr; k++)
					{
						for (l = 0; l < nbpc; l++)
						{
							for (o = 0; o < MB; o++)
							{
								for (p = 0; p < NB; p++)
								{
									SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
								}
							}
						}
					}
					//Send submatrix
					if ((i != 0)||(j != 0))
					{
						dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
					}
				}
				if (nbigc > 0)
				{
					j = nbigc - 1;
					{
						//Big process row, tail process column
						for (k = 0; k <= nbpr; k++)
						{
							for (l = 0; l < nbpc; l++)
							{
								for (o = 0; o < MB; o++)
								{
									for (p = 0; p < NB; p++)
									{
										SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
									}
								}
							}
							l = nbpc;
							for (o = 0; o < MB; o++)
							{
								for (p = 0; p < ntailc; p++)
								{
									SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
								}
							}
						}
						//Send submatrix
						if ((i != 0)||(j != 0))
						{
							dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
						}
					}

					for (j = nbigc - 2; j >= 0; j--)
					{
						//Big process row, big process column
						for (k = 0; k <= nbpr; k++)
						{
							for (l = 0; l <= nbpc; l++)
							{
								for (o = 0; o < MB; o++)
								{
									for (p = 0; p < NB; p++)
									{
										SUB_A[(o + k * MB) + (p + l * NB) * MXLLDA] = A[(o + i * MB + k * MB * NPROW) * N + (p + j * NB + l * NB * NPCOL)];
									}
								}
							}
						}
						//Send submatrix
						if ((i != 0)||(j != 0))
						{
							dgesd2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &i, &j);
						}
					}
				}
			}
		}
	}
	else
	{
		//Recv submatrix
		dgerv2d_(&ICTXT, &MXLLDA, &MXLOCC, SUB_A, &LDA, &RSRC, &CSRC);
	}

	//Matrix B
	if (MYCOL == 0)
	{
		if(MYROW == 0)
		{
			for (i = NPROW - 1; i >= nbigr; i--)
			{
				for (k = 0; k < nbpr; k++)
				{
					for (o = 0; o < MB; o++)
					{
							SUB_B[(o + k * MB)] = B[o + i * MB + k * MB * NPROW];
					}
				}
				//Send submatrix
				if (i != 0)
				{
					dgesd2d_(&ICTXT, &MXLLDB, &MXRHSC, SUB_B, &LDA, &i, &MYCOL);
				}
			}

			if (nbigr > 0)
			{

				i = nbigr - 1;
				for (k = 0; k < nbpr; k++)
				{
					for (o = 0; o < MB; o++)
					{
						SUB_B[(o + k * MB)] = B[o + i * MB + k * MB * NPROW];
					}
				}
				k = nbpr;
				for (o = 0; o < ntailr; o++)
				{
					SUB_B[(o + k * MB)] = B[o + i * MB + k * MB * NPROW];
				}
				//Send submatrix
				if (i != 0)
				{
					dgesd2d_(&ICTXT, &MXLLDB, &MXRHSC, SUB_B, &LDA, &i, &MYCOL);
				}

				for (i = nbigr - 2; i >= 0; i--)
				{
					for (k = 0; k <= nbpr; k++)
					{
						for (o = 0; o < MB; o++)
						{
							SUB_B[(o + k * MB)] = B[o + i * MB + k * MB * NPROW];
						}
					}
					//Send submatrix
					if (i != 0)
					{
						dgesd2d_(&ICTXT, &MXLLDB, &MXRHSC, SUB_B, &LDA, &i, &MYCOL);
					}
				}
			}
		}
		else
		{
			//Recv submatrix
			dgerv2d_(&ICTXT, &MXLLDB, &MXRHSC, SUB_B, &LDA, &RSRC, &CSRC);
		}
	}

	//if ((MYCOL == 0)&&(MYROW == 0))
	//{
	//	for (i = 0; i < MXLLDA * MXLOCC; i++)
	//	{
	//		printf("%e ", SUB_A[i]);
	//	}
	//	printf("\n");
	//
	//	for (i = 0; i < MXLLDB * MXRHSC; i++)
	//	{
	//		printf("%e ", SUB_B[i]);
	//	}
	//	printf("\n");
	//}

	//*****
	//Call the SCALAPACK routine
	//*****
	pdgesv_(&N, &NRHS, SUB_A, &IA, &JA, DESCA, IPIV, SUB_B, &IB, &JB, DESCB, &INFO);

	//if ((MYCOL == 0)&&(MYROW == 1))
	//{
	//	for (i = 0; i < MXLLDB * MXRHSC; i++)
	//	{
	//		printf("%e ", SUB_B[i]);
	//	}
	//	printf("\n");
	//}

	//Collect and print the ansver
	//Direct order to take (0,0) from local memory
	if (MYCOL == 0)
	{
		if(MYROW == 0)
		{

			if (nbigr > 0)
			{
				for (i = 0; i < nbigr - 1; i++)
				{
					//Recv submatrix
					if (i != 0)
					{
						dgerv2d_(&ICTXT, &MXLLDB, &MXRHSC, SUB_B, &LDA, &i, &MYCOL);
					}

					for (k = 0; k <= nbpr; k++)
					{
						for (o = 0; o < MB; o++)
						{
							B[o + i * MB + k * MB * NPROW] = SUB_B[(o + k * MB)];
						}
					}
				}

				i = nbigr - 1;
				//Recv submatrix
				if (i != 0)
				{
					dgerv2d_(&ICTXT, &MXLLDB, &MXRHSC, SUB_B, &LDA, &i, &MYCOL);
				}
				for (k = 0; k < nbpr; k++)
				{
					for (o = 0; o < MB; o++)
					{
						B[o + i * MB + k * MB * NPROW] = SUB_B[(o + k * MB)];
					}
				}
				k = nbpr;
				for (o = 0; o < ntailr; o++)
				{
					B[o + i * MB + k * MB * NPROW] = SUB_B[(o + k * MB)];
				}
			}

			for (i = nbigr; i < NPROW; i++)
			{
				//Recv submatrix
				if (i != 0)
				{
					dgerv2d_(&ICTXT, &MXLLDB, &MXRHSC, SUB_B, &LDA, &i, &MYCOL);
				}
				for (k = 0; k < nbpr; k++)
				{
					for (o = 0; o < MB; o++)
					{
						B[o + i * MB + k * MB * NPROW] = SUB_B[(o + k * MB)];
					}
				}
			}
		}
		else
		{
			//Send submatrix
			dgesd2d_(&ICTXT, &MXLLDB, &MXRHSC, SUB_B, &LDA, &RSRC, &CSRC);
		}
	}


	//Distribution cleanup
	free(SUB_A);
	free(SUB_B);
	free(IPIV);

	//Free matrices A and B
	if ((MYROW == 0)&&(MYCOL == 0))
	{
		printf("Matrix X:\n");
		for (i = 0; i < N; i++)
		{
			printf("%e ", B[i]);
		}
		printf("\n");

		free(A);
		free(B);
	}

	//Release the process grid
	blacs_gridexit_(&ICTXT);
	blacs_exit_(&blacs_exitcode);

	return 0;
}
