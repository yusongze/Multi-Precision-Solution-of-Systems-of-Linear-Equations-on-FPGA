/*
 * LU.c
 *
 *  Created on: 3 Mar 2018
 *      Author: songze
 */
#include "LU.h"
#include "time.h"
#ifdef fp_count
extern unsigned int fp_counter;
#endif

/***********************************************************************
 * LU factorization
 * Input coefficient matrix A[][], b[] and x[], passing pivoting matrix P
 * Output L, U. A
 *
 ***********************************************************************/
void LU_factorization(	float A[N][N],
						float PA[N][N],
						float Pb[N]){
	int i,j,k;
	int Index = 0;
	float swap[N];
	float swapPA[N];
	float swap1;
#ifdef show_P
	float P[N][N]={0.0};	//testing

	for(i=0;i<N;i++){
		P[i][i] = 1.0;
	}	//testing
#endif
	/********************************************************
	 * Forward elimination without partial pivoting to get L, U*
	 ********************************************************/


	//Forward elimination
	Pivoting_element: for(k=0;k<N-1;k++){
		Index =  maxIndex(A, k);	//get the row index of the max element
		if(Index != k){
//			pivoting_A(A_RX, k, Index);	//swap the pivoting row with the row where max element is
			swap1 = Pb[k];
			Pb[k] = Pb[Index];
			Pb[Index] = swap1;
			LU_factorization_label0:
			for(i=0;i<N;i++){
					swap[i] = A[k][i];
					A[k][i] = A[Index][i];
					A[Index][i] = swap[i];
					swapPA[i] = PA[k][i];
					PA[k][i] = PA[Index][i];
					PA[Index][i] = swapPA[i];
#ifdef show_P
					swap[i] = P[k][i];	testing
					P[k][i] = P[Index][i];
					P[Index][i] = swap[i];
#endif
			#ifdef fp_count
					fp_counter = fp_counter + 6;
			#endif
			}

//			pivoting_A(Pivoting_Matrix, k, Index);	//generate row swapping matrix

		}
		Rows:for(i=k+1;i<N;i++){
			A[i][k] = A[i][k] / A[k][k];
#ifdef fp_count
			fp_counter = fp_counter + 1;
#endif
				//print_out(A, b, x);			//testing
			sub_mul_elimi:for(j=k+1;j<N;j++){
				A[i][j] = A[i][j] - A[i][k]*A[k][j];
#ifdef fp_count
				fp_counter = fp_counter + 2;
#endif
				//print_out(A, b, x);			//testing
			}
		}
	}
#ifdef show_P
	printf("/nP ");		//testing
	print_out_matrixNbyN(P);	//testing
#endif

//	printf("/nPA ");		//testing
//	print_out_matrixNbyN(PA);	//testing
}


void LU_solver(	float A[N][N],
				float b_RX[N],
				float x[N]){
	int i, j;
	float s;	//row partial sum.
	float y[N];


	//Two steps to solve LUx=b
	//Forward solve Ly=b using L and U
	LU_solver_label29:
	for(i=0;i<N;i++){
		s = 0.0;	//intermediate sum initialisation
#ifdef fp_count
		fp_counter = fp_counter + 1;
#endif
		LU_solver_label3:
		for(j=0;j<i;j++){
			s += A[i][j]*y[j];	//L*y
#ifdef fp_count
			fp_counter = fp_counter + 2;
#endif
		}
		y[i] = b_RX[i]-s;
#ifdef fp_count
		fp_counter = fp_counter + 1;
#endif
	}

#ifdef test
	printf("\nNo. 0th y ");
	print_out_matrixNby1_single(y); //print x
#endif

	//Backward solve Ux=y using L and U
	LU_solver_label30:
	for(i=N-1;i>=0;i--){
		s = 0.0;	//intermediate sum initialisation
#ifdef fp_count
		fp_counter = fp_counter + 1;
#endif
		LU_solver_label4:
		for(j=i+1;j<N;j++){
			s += A[i][j]*x[j];
#ifdef fp_count
			fp_counter = fp_counter + 2;
#endif
		}
		x[i] = (y[i]-s)/A[i][i];
#ifdef fp_count
		fp_counter = fp_counter + 2;
#endif
	}

};


void LU_solver_refinement(	float A_RX[N][N],
							float b_RX[N],
							float A[N][N],
							float x[N]){
	int i, j;
	float s;	//row partial sum.
	float s_r;
	float r[N];
	float y[N];
	float d[N];

	//Two steps to solve LUx=b
	//Forward solve Ly=b using L and U
	Rows:for(i=0;i<N;i++){
		s = 0.0;	//intermediate sum initialisation
		s_r = 0.0;
#ifdef fp_count
		fp_counter = fp_counter + 2;
#endif
		Columns:for(j=0;j<N;j++){
			if(j<i){
				s += A[i][j]*y[j];
			}
			s_r += A_RX[i][j]*x[j];
#ifdef fp_count
			fp_counter = fp_counter + 4;
#endif
		}
		r[i] = b_RX[i]-s_r;		//compute residue
		y[i] = r[i]-s;
#ifdef fp_count
		fp_counter = fp_counter + 2;
#endif

	}


	//Backward solve Ux=y using L and U
	Backward_Rows:for(i=N-1;i>=0;i--){
		s = 0.0;	//intermediate sum initialisation
#ifdef fp_count
		fp_counter = fp_counter + 1;
#endif
		Backward_Columns:for(j=i+1;j<N;j++){
			s += A[i][j]*d[j];
#ifdef fp_count
			fp_counter = fp_counter + 2;
#endif
		}
		d[i] = (y[i]-s)/A[i][i];
		x[i] = x[i] + d[i];		//update x
#ifdef fp_count
		fp_counter = fp_counter + 3;
#endif
	}

//testing
	float norm_r = 0.0;
		norm_2(r, &norm_r);
			printf("\nSolution ");
			print_out_matrixNby1_single(x); //print x
			printf("\nResidue ");
			print_out_matrixNby1(r);
			printf("Residue norm is ");
			printf("%.30f", norm_r);
			putchar('\n');
};


void Iterative_refine(	float A_RX[N][N],
						float b_RX[N],
						float x[N],
						float A[N][N]){
	int i;

	Refinement: for(i=0;i<Iteration_N;i++){
			//get_residue(A_RX, x, b_RX, r);


	//			write_to_file(norm_r, "norm_r.dat");

			LU_solver_refinement(A_RX, b_RX, A, x);
		}
}

void get_residue(float A_RX[N][N], float x[N], float b_RX[N], float r[N]){
	float s;

	//compute residue in float precision

	Rows:for(int i=0;i<N;i++){
		s = 0.0;
#ifdef fp_count
		fp_counter = fp_counter + 1;
#endif
		Columns:for(int j=0;j<N;j++){
			s += A_RX[i][j]*x[j];
#ifdef fp_count
			fp_counter = fp_counter + 2;
#endif
		}
		r[i] = b_RX[i]-s;
#ifdef fp_count
		fp_counter = fp_counter + 1;
#endif
	}
};

void update_x(float x[N], float d[N]){
	update_x_label36:
	for(int i=0;i<N;i++){
		x[i] = x[i]+d[i];
#ifdef fp_count
		fp_counter = fp_counter + 1;
#endif
	}
};

void LU_extraction(float A[N][N],
				float L[N][N],
				float U[N][N]){

	int i,j;

#ifdef fp_count
	fp_counter = fp_counter + N;
#endif

	//U is copied from A

	Rows:for(i=0;i<N;i++){
		Columns:for(j=0;j<N;j++){
			/**********Extract L, U************/
			if(j<=i-1){
				L[i][j] = A[i][j];
				U[i][j] = 0.0;
			}
			else{
				U[i][j] = A[i][j];
				L[i][j] = 0.0;
			}
			/**********************************/
#ifdef fp_count
			fp_counter = fp_counter + 2;
#endif
		}
		L[i][i] = 1.0;
#ifdef fp_count
		fp_counter = fp_counter + 1;
#endif
	}
//	print_out_matrixNby1_single(b_RX);	//test
//	print_out_matrixNby1_single(Pb);	//test

}

//
//void get_LU(float A[N][N],
//			float L[N][N],
//			float U[N][N]){
//	//Get L, U which are stored in A_RX, and get the matrix used to pivot
//	LU_factorization(A);
//#ifdef test
//		printf("A after LU factorization ");
//		print_out_matrixNbyN(A);
//#endif
//	//extract L and U from A_RX and change the problem to P*A*x = P*b which is also L*U*x = P*b
//	LU_extraction(A, L, U);
//
//#ifdef test
//		printf("L ");
//		print_out_matrixNbyN(L);
//		printf("U ");
//		print_out_matrixNbyN(U);
//
//		float diff[N][N] = {0.0};
//		float LU[N][N] = {0.0};
//		for(int ii=0;i<N;i++){
//			for(int jj=0;jj<N;jj++){
//				float sum = 0.0;
//				for(int kk=0;kk<N;kk++){
//					sum = sum + L[ii][kk]*U[kk][jj];
//				}
//			LU[ii][jj] = sum;
//			diff[ii][jj] = A_RX[ii][jj] - LU[ii][jj];
//			}
//		}
//		printf("Difference from A_RX and LU is ");
//		print_out_matrixNbyN(diff);
//
//#endif
//}
//


/*void matrix_gen(float A[N][N], float PA[N][N], float b[N] ,float Pb[N]){
	int i,j;

	srand(time(NULL));
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i][j] = (float)rand()/(float)(RAND_MAX/1000);
		}
		b[i] = (float)rand()/(float)(RAND_MAX/1000);
	}
	print_out_matrixNbyN(A);
	print_out_matrixNby1(b);
}

void matrix_gen2(float A[N][N], float PA[N][N], float b[N] ,float Pb[N]){
	int i,j;

	srand(time(NULL));
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i][j] = (float)rand()/(float)(RAND_MAX/1000);
		}
		b[i] = (float)rand()/(float)(RAND_MAX/1000);
	}
	print_out_matrixNbyN(A);
	print_out_matrixNby1(b);
}*/


//Note: use cordic to realize sqrt
void norm_2(float matrix[N], float *result){
	int i;

#ifdef fp_count
	fp_counter = fp_counter + 1;
#endif

	norm_2_label32:
	for (i=0;i<N;i++){
		*result = *result + matrix[i]*matrix[i];
#ifdef fp_count
		fp_counter = fp_counter + 2;
#endif
	}
#ifdef fp_count
	fp_counter = fp_counter + 1;
#endif
	*result = sqrt(*result);
};



int maxIndex(float A[N][N], int columnIndex){
	int i, maxIndex;
	float maxValue;

	maxIndex = columnIndex;
	maxValue = abs(A[columnIndex][columnIndex]);
	maxIndex_label33:
	for(i=0;i<N;i++){
		if(i>=columnIndex){
			if( abs(maxValue) < abs(A[i][columnIndex]) ){
						maxValue = A[i][columnIndex];
						maxIndex = i;
			#ifdef fp_count
					fp_counter = fp_counter + 6;
			#endif
			}
		}
#ifdef fp_count
		fp_counter = fp_counter + 3;
#endif

	}

	return maxIndex;

};


void pivoting_A(float A[N][N], int k, int Index){
	int i;
	float swap[N];


	/**********************
	 * L is stored
	 ************************/
	pivoting_A_label34:
	for(i=0;i<N;i++){
		swap[i] = A[k][i];
		A[k][i] = A[Index][i];
		A[Index][i] = swap[i];
#ifdef fp_count
		fp_counter = fp_counter + 3;
#endif
	}
};

void pivoting_b(float b[N], int k, int Index){
	float swap;

	swap = b[k];
	b[k] = b[Index];
	b[Index] = swap;
#ifdef fp_count
	fp_counter = fp_counter + 3;
#endif

};



void print_out_matrixNbyN(float A[N][N]){
	int i,j;

	printf("Matrix is \n");
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				printf("%e\t\t", float(A[i][j]));
			}
			putchar ('\n');
		}
		putchar ('\n');
}

void print_out_matrixNby1(float A[N]){
	int i;

	printf("Matrix is \n");
		for(i=0;i<N;i++){
			printf("%e\t\t", float(A[i]));
			}
		printf("\n");
}

void print_out_matrixNby1_single(float A[N]){
	int i;

	printf("Matrix is \n");
		for(i=0;i<N;i++){
			printf("%e\t\t", float(A[i]));
			}
		printf("\n");
}

void print_A_L_U(float A[N][N], float L[N][N], float U[N][N]){
	//print A, L and U
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				printf("%e\t\t", float(A[i][j]));
			}
		putchar ('\n');
		}

	//print A, L and U
	printf("L is \n");
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				printf("%e\t\t", float(L[i][j]));
			}
		putchar ('\n');
		}
	printf("U is \n");
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				printf("%e\t\t", float(U[i][j]));
			}
		putchar ('\n');
		}
}

int write_to_file(float norm_r[Iteration_N], char const *fileName)
{
	int i;
	FILE *f = fopen(fileName, "w+");
	if (f == NULL) return -1;
	for(i=0;i<Iteration_N;i++) {
		// you might want to check for out-of-disk-space here, too
		fprintf(f, "%.52f", norm_r[i]);
		fprintf(f, "\n");
	}
	fclose(f);
	return 0;
}
