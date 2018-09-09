#include <stdio.h>
#include <stdlib.h>
#include "LU.h"
#include "data100by100.h"

#ifdef fp_count
extern unsigned int fp_counter = 0;
#endif
int main(){
	int i,j;
	float A_TX[N][N] = {0.0};
	float b_TX[N] = {0.0};
	float r[N] = {0.0};
	float x[N] = {0.0};
	unsigned int x_u32[N] = {0};
	float norm_RX = 0.0;
	data_RX_type data_RX;

	for(j=0;j<N;j++){
		for(i=0;i<N;i++){
//			A_TX[2*i][j]	= (A_random[i][j] & 0xFFFFFFFF00000000) >> 32;		//first 32bit put in A[i]
//			A_TX[2*i+1][j]	=  A_random[i][j] & 0x00000000FFFFFFFF;		//last 32bit put in A[i+1]
			A_TX[i][j] = A_random[i][j];

		}
	}

	for(i=0;i<N;i++){
//		b_TX[2*i]	= (b_random[i] & 0xFFFFFFFF00000000) >> 32;
//		b_TX[2*i+1]	=  b_random[i] & 0x00000000FFFFFFFF;
		b_TX[i] = b_random[i];
	}


//	printf("Matrix A_TX is \n");
//			for(i=0;i<N;i++){
//				for(j=0;j<N;j++){
//					printf("%e\t\t", A_TX[i][j]);
//				}
//				putchar ('\n');
//			}
//	putchar ('\n');
//
//	printf("Matrix b_TX is \n");
//				for(i=0;i<N;i++){
//					printf("%e\t\t", b_TX[i]);
//					putchar ('\n');
//				}
//	putchar ('\n');

		MP(A_TX, b_TX, x_u32);
		for(i=0;i<N;i++){
			data_RX.u32 = x_u32[i];
			x[i] = data_RX.f32;
		}
		get_residue(A_TX, x, b_TX, r);
		norm_2(r, &norm_RX);

//	printf("Matrix A_TX after computation is \n");
//			for(i=0;i<N;i++){
//				for(j=0;j<N;j++){
//					printf("%e\t\t", A_TX[i][j]);
//				}
//				putchar ('\n');
//			}
//	putchar ('\n');
//
//	printf("Matrix b_TX after computation is \n");
//				for(i=0;i<N;i++){
//					printf("%e\t\t", b_TX[i]);
//					putchar ('\n');
//				}
//	putchar ('\n');
//
	printf("Solution x is \n");
				for(i=0;i<N;i++){
					printf("%e\t\t", x[i]);
					putchar ('\n');
				}
	putchar ('\n');

#ifdef fp_count
	printf("Number of floating point operation is %d\n", fp_counter);
	printf("MFlops is %f\n", float(fp_counter)/	9);		//with optimization
//	printf("MFlops is %f\n", float(fp_counter)/66.33);		//without optimization
//	printf("norm_final is %x\n", data_RX.u32);
#endif
	printf("norm at last iteration stage is %e\n", norm_RX);
	return 0;



}
