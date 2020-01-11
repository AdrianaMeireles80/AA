#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>


#define SIZE 32
#define block_size 16
//#define num_blocks (SIZE/block_size)

double clearcache [30000000];

float A[SIZE][SIZE] __attribute__((aligned(16)));
float B[SIZE][SIZE] __attribute__((aligned(16)));
float C[SIZE][SIZE] __attribute__((aligned(16)));
float aux2[SIZE][SIZE] __attribute__((aligned(16)));

void clearCache () {
	unsigned i;
	for ( i = 0; i < 30000000; i++)
		clearcache[i] = i;
}

void fillMatrices(){

	int i,j;

	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			
			A[i][j] = (float) (rand()%10);
			B[i][j] = 1;
			C[i][j] = 0;
			aux2[i][j] = 0;
		}
	}
}

void matrixMult() {
	int i,j,k;
    //float x;

	for (i = 0; i < SIZE; i++) {
		for ( j = 0; j < SIZE; j++) {
			//x= 0;
			for (k = 0; k < SIZE; k++) {
				//x += A[i][k] * B[k][j];
				aux2[i][j] += A[i][k] * B[k][j];
			}
		   //C[i][j] += x;
		   
		
			//printf("%f\n",C[i][j]);
		}

	}
}

/*
//X
void matrixMult_block(){
	int i,j,jj,kk;

    float x[16]={0};

    for(jj=0; jj<SIZE; jj += block_size)
    {
        for(kk=0; kk<SIZE; kk += block_size)
        {
                for(i=0;i<SIZE;i++)
                {
                    for(j = jj; j< (jj + block_size); j++)
                    {
                            x[0] = A[i][kk] * B[kk][j];
                            x[1] = A[i][kk+1] * B[kk+1][j];
                            x[2] = A[i][kk+2] * B[kk+2][j];
                            x[3] = A[i][kk+3] * B[kk+3][j];
                            x[4] = A[i][kk+4] * B[kk+4][j];
                            x[5] = A[i][kk+5] * B[kk+5][j];
                            x[6] = A[i][kk+6] * B[kk+6][j];
                            x[7] = A[i][kk+7] * B[kk+7][j];
                            x[8] = A[i][kk+8] * B[kk+8][j];
                            x[9] = A[i][kk+9] * B[kk+9][j];
                            x[10] = A[i][kk+10] * B[kk+10][j];
                            x[11] = A[i][kk+11] * B[kk+11][j];
                            x[12] = A[i][kk+12] * B[kk+12][j];
                            x[13] = A[i][kk+13] * B[kk+13][j];
                            x[14] = A[i][kk+14] * B[kk+14][j];
                            x[15] = A[i][kk+15] * B[kk+15][j];

                            C[i][j] += x[0] + x[1] + x[2] + x[3] + x[4]+ x[5]+ x[6]+ x[7] + x[8] + x[9]+ x[10]+ x[11]+ x[12]+ x[13] + x[14]+ x[15];
                    }
                }
        }
    }
}
//Dolly
void mult_Block_Vec(){
	int i = 0, j = 0, jj = 0, kk = 0, z = 0;
	
	float x[16],aux[16];

	for(jj = 0; jj < SIZE; jj += block_size){
		for(kk = 0; kk < SIZE; kk += block_size){
			for(i = 0; i < SIZE; i++){
				for(j = jj, z = 0; j < ((jj + block_size) > SIZE ? SIZE : (jj + block_size)); j++, z++){
							x[0] = A[i][kk] * B[kk][j];
                            x[1] = A[i][kk+1] * B[kk+1][j];
                            x[2] = A[i][kk+2] * B[kk+2][j];
                            x[3] = A[i][kk+3] * B[kk+3][j];
                            x[4] = A[i][kk+4] * B[kk+4][j];
                            x[5] = A[i][kk+5] * B[kk+5][j];
                            x[6] = A[i][kk+6] * B[kk+6][j];
                            x[7] = A[i][kk+7] * B[kk+7][j];
                            x[8] = A[i][kk+8] * B[kk+8][j];
                            x[9] = A[i][kk+9] * B[kk+9][j];
                            x[10] = A[i][kk+10] * B[kk+10][j];
                            x[11] = A[i][kk+11] * B[kk+11][j];
                            x[12] = A[i][kk+12] * B[kk+12][j];
                            x[13] = A[i][kk+13] * B[kk+13][j];
                            x[14] = A[i][kk+14] * B[kk+14][j];
                            x[15] = A[i][kk+15] * B[kk+15][j];

					aux[z] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[10] + x[11] + x[12] + x[13] + x[14] + x[15];
				}

				for(j = jj, z = 0; j < ((jj + block_size) > SIZE ? SIZE : (jj + block_size)); j++, z++){
					C[i][j] = aux[z];
				}
			}
		}
	}
}

*/


//Resende
void matrixMult_block() {
	int i = 0, j = 0, jj = 0, kk = 0, z = 0;
  //int block_size = 16;
  float x[16], aux[16];

  for (jj = 0; jj < SIZE; jj += block_size) {
    for (kk = 0; kk < SIZE; kk += block_size) {
      for (i = 0; i < SIZE; i++) {
        for (j = jj, z = 0; j < ((jj + block_size) > SIZE ? SIZE : (jj + block_size)); j++, z++) {
          x[0] = A[i][kk] * B[kk][j];
          x[1] = A[i][kk+1] * B[kk+1][j];
          x[2] = A[i][kk+2] * B[kk+2][j];
          x[3] = A[i][kk+3] * B[kk+3][j];
          x[4] = A[i][kk+4] * B[kk+4][j];
          x[5] = A[i][kk+5] * B[kk+5][j];
          x[6] = A[i][kk+6] * B[kk+6][j];
          x[7] = A[i][kk+7] * B[kk+7][j];
          x[8] = A[i][kk+8] * B[kk+8][j];
          x[9] = A[i][kk+9] * B[kk+9][j];
          x[10] = A[i][kk+10] * B[kk+10][j];
          x[11] = A[i][kk+11] * B[kk+11][j];
          x[12] = A[i][kk+12] * B[kk+12][j];
          x[13] = A[i][kk+13] * B[kk+13][j];
          x[14] = A[i][kk+14] * B[kk+14][j];
          x[15] = A[i][kk+15] * B[kk+15][j];

          aux[z] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[10] + x[11] + x[12] + x[13] + x[14] + x[15];
      }

        for (j = jj, z= 0; j < ((jj + block_size) > SIZE ? SIZE : (jj + block_size)); j++, z++) {
          C[i][j] += aux[z];
        }
      }
    }
  }
}


void compare() {
	int i,j, v = 1;
	for(i = 0;i < SIZE; i++) 
		for(j=0;j<SIZE;j++){
			printf("RES: %f -> RESEQ: %f\n",aux2[i][j],C[i][j]);
			if(aux2[i] != C[i]) {
				v = 0;
				printf("ERRRRRROoooooooooooooooooooooooooooooooooooooooooooo\n");
			}
		}
	printf("DEU: %d\n",v);
	}





int main()
{
	fillMatrices();
	matrixMult();
	matrixMult_block();
	//mult_Block_Vec();
	compare();
}