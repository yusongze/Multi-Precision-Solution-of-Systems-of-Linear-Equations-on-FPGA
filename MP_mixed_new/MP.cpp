/*
 ============================================================================
 Name        : MP.cpp
 Author      : Songze
 Version     :
 Copyright   : Your copyright notice
 Description : Algorithm 1.1 using MPFR under mixed-precision floating point
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "LU.h"

#ifdef fp_count
extern unsigned int fp_counter;
#endif


//AXI is 32bit wide, so ports cannot exceed 32bit
void MP(float A_RX[N][N],
		float b_RX[N],
		unsigned int result[N]) {
#pragma HLS PIPELINE

//#pragma HLS ARRAY_PARTITION variable=result complete dim=1
#pragma HLS INTERFACE s_axilite port=result bundle=yakisoba
#pragma HLS ARRAY_PARTITION variable=A_RX complete dim=0
#pragma HLS ARRAY_PARTITION variable=b_RX complete dim=1
#pragma HLS INTERFACE s_axilite port=return bundle=yakisoba
#pragma HLS INTERFACE s_axilite port=A_RX bundle=yakisoba
#pragma HLS INTERFACE s_axilite port=b_RX bundle=yakisoba


	int i,j,k;
	data_TX_type data_TX;

	half A[N][N];
	float PA[N][N];
	float Pb[N];
	float x[N];

#pragma HLS ARRAY_PARTITION variable=A complete dim=2
#pragma HLS ARRAY_PARTITION variable=PA complete dim=2
#pragma HLS ARRAY_PARTITION variable=Pb complete dim=1
#pragma HLS ARRAY_PARTITION variable=x complete dim=1

#ifdef fp_count
	fp_counter = fp_counter + N*N + 2*N;
#endif


	/****************Declare matrix in floating point*************/
	//A[N]:		coefficient matrix. After LU factorization, it stores L and U
	//PA[N]:	A[N] after pivoting
	//b[N]:		right-hand matrix
	//Pb[N]:	b[N] after pivoting
	//r[Iteration_N][N]:	residue matrix
	//norm_r[Iteration_N]:	norm using Euclidean distance
	//d[N]:		update matrix
	/************************************************************/

	Store_A_RX_Rows:for(i=0;i<N;i++){
		Pb[i] = b_RX[i];
		Store_A_RX_Columns:for(j=0;j<N;j++){
			A[i][j] = A_RX[i][j];
			PA[i][j] = A_RX[i][j];
#ifdef fp_count
			fp_counter = fp_counter + 1;
#endif
		}
	}

/******************************Generate random coefficient matrix A and b******************************/
	//matrix_gen(A, PA, b, Pb);		//generate nonsingular N*N coefficient  matrix, and matrix b;



/*******************************LU factorization *****************************/
	LU_factorization(A, PA, Pb);
//testing
//	printf("Matrix A(LU) is \n");
//			for(i=0;i<N;i++){
//				for(j=0;j<N;j++){
//					printf("%e\t\t", A[i][j]);
//				}
//				putchar ('\n');
//			}
//	putchar ('\n');
//
//	printf("Matrix Pb is \n");
//				for(i=0;i<N;i++){
//					printf("%e\t\t", Pb[i]);
//					putchar ('\n');
//				}
//	putchar ('\n');
/***********************Initial solution*******************************/
	LU_solver(A, Pb, x);

//testing
	double r[N] = {0.0};
	for(i=0;i<N;i++){
		float s_r=0.0;
		for(j=0;j<N;j++){
			s_r += A_RX[i][j]*x[j];
		}
		r[i] = b_RX[i]-s_r;
	}
	float norm_r = 0.0;
	norm_2(r, &norm_r);
		printf("\nSolution ");
		print_out_matrixNby1_single(x); //print x
//		printf("\nResidue ");
//		print_out_matrixNby1(r);
		printf("Residue norm is ");
		printf("%.30f", norm_r);
		putchar('\n');

/***********************Iterative refinement***************/
	Iterative_refine(PA, Pb, x, A);

#ifdef test
//	Refinement: for(i=0;i<Iteration_N;i++){
//		get_residue(A_RX, x, b_RX, r);
//	//		norm_r[i] = norm_2(r[i]);
//	#ifdef test
//		int j=0;
//		float norm_r = 0.0;
//			for (j=0;j<N;j++){
//				norm_r = norm_r + r[j]*r[j];
//			}
//			printf("\nNo. %d Solution ", i);
//			print_out_matrixNby1_single(x); //print x
//			printf("\n\nNo. %d Residue ", i);
//			print_out_matrixNby1(r);
//			printf("No. %d Residue norm is ", i);
//			printf("%e", sqrt(norm_r));
//			putchar('\n');
//	#endif
//
//	//			write_to_file(norm_r, "norm_r.dat");
//
//		LU_solver_refinement(L, U ,r ,d);
//		update_x(x, d);
//	}
#endif
/***********************************************/

#ifdef test
	printf("\nSolution ");
	print_out_matrixNby1_single(x); //print x
#endif
//	printf("Pivoting Matrix is\n");
//	print_A_L_U(Pivoting_Matrix, L, U);	//print A, L, U


	Output: for(i=0;i<N;i++){
		data_TX.f32 = x[i];
		result[i] = data_TX.u32;
#ifdef fp_count
		fp_counter = fp_counter + 2;
#endif
	}

}

