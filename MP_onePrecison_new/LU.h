/*
 * LU.h
 *
 *  Created on: 28 Jan 2018
 *      Author: songze
 */

#ifndef LU_H_
#define LU_H_
//#define fp_count
//#define test

//#include "/opt/Xilinx/Vivado/2018.2/include/gmp.h"
//#include "/opt/Xilinx/Vivado/2018.2/include/mpfr.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

#include "math.h"
#include "string.h"

//typedef ap_fixed<16,5> float;
//typedef ap_fixed<32,8> float;
//typedef ap_fixed<64,11> float;


#define N 100
#define Iteration_N 10


union data_TX_type
	{
		unsigned int u32;
		float f32;
	};

union data_RX_type
	{
		unsigned int u32;
		float f32;
	};

void MP(float A_RX[N][N],
		float b_RX[N],
		unsigned int result[N]);

void LU_factorization(	float A[N][N],
						float PA[N][N],
						float Pb[N]);		//Gaussian elimination with partial pivoting
void LU_extraction(float A[N][N],
				float L[N][N],
				float U[N][N]);
void get_LU(float A[N][N],
			float L[N][N],
			float U[N][N]);

void LU_solver(	float A[N][N],
				float b_RX[N],
				float x[N]);

void LU_solver_refinement(	float A_RX[N][N],
							float b_RX[N],
							float A[N][N],
							float x[N]);

void norm_2(float matrix[N], float *result);		//Euclidean norm

void get_residue(float A_RX[N][N], float x[N], float b_RX[N], float r[N]);	//compute residue in float precision
void update_x(float x[N], float d[N]);
void Iterative_refine(	float A_RX[N][N],
						float b_RX[N],
						float x[N],
						float A[N][N]);

int maxIndex(float A[N][N], int columnIndex);
void pivoting_A(float A[N][N], int k, int Index);
void pivoting_b(float b[N], int k, int Index);
void matrix_gen(float A[N][N], float PA[N][N], float b[N] ,float Pb[N]);
void print_out_matrixNbyN(float A[N][N]);
void print_out_matrixNby1(float A[N]);
void print_out_matrixNby1_single(float A[N]);
void print_A_L_U(float A[N][N], float L[N][N], float U[N][N]);

int write_to_file(float norm_r[Iteration_N], char const *fileName);

#endif /* LU_H_ */
