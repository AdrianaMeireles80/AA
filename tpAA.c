#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <immintrin.h>
#include <papi.h>
#include <unistd.h>

#define NUM_EVENTS 2
struct timeval begin;
//int Events[NUM_EVENTS] = {PAPI_L2_DCR,PAPI_LD_INS};
//int Events[NUM_EVENTS]={PAPI_L3_DCR,PAPI_L2_DCR};
//int Events[NUM_EVENTS]={PAPI_L3_TCM,PAPI_L3_TCA};
//int Events[NUM_EVENTS] = {PAPI_L3_TCM,PAPI_TOT_INS};

//int EventSet = PAPI_NULL;
//long long values[NUM_EVENTS];


/*
1.1

No nodo 662 podem ser identificados vários contadores para medir a performance sendo os mais relevantes para o caso
PAPI_FP_INS,PAPI_SR/LD_INS,PAPI_BR_INS,PAPI_SP/DP_VEC referentes às instruções ;PAPI_FP_OPS,PAPI_SP/DP_OPS,
PAPI_TOT_INS referentes a FLOPS e intensidade operacional e,por último, PAPI_L1/2/3_TCM que se refere às cache misses.

--->Relativamente às instruções é preciso ter em conta vários aspetos tais como se são floating ou não que é mais leve, se trata de
instruções store ou load, estas últimas como trasferem dados da memória para os registos é mais pesado porque é mais custoso 
ir à memória o que vai afetar a performance do xo de execução. Também nos diz se trata de instruções SIMD(single instruction
multiple data)  single/double precision que possuem apenas um fluxo de instrução, mas possuem múltiplas unidades de cálculo.
 Isso significa que a máquina é capaz de executar uma mesma instrução em um conjunto de dados de maneira simultânea o que permite
 obter xos de execução menores. 

-->Relativamente às operações permite contar operações escalares de single/double precision o que nos permite identificar possíveis
bottlenecks pois caso algo seja executado vários vezes, se for ineficiente pode atrasar o programa todo, aumentando o xo de 
execução o que afeta a performance.

-->Relativamente às cache misses assume-se que a CPU deve aguardar um carregamento da memória principal e o número total de "stall cycles" depende 
das cache misses e da miss penalty(Memory stall cycles = Memory accesses x miss rate x miss penalty) que por sua vez vão afetar o
xo de execução(CPU time = (CPU execution cycles + Memory stall cycles) x Cycle time).

1.4

L1 cache:             32K
L2 cache:             256K
L3 cache:             30720K
RAM:                  4*30720K

Ki=1024 bytes

(i) N*N*3*4 < 32768 -> N < 52--------------------------> N=32      |
(ii) N*N*3*4 <  262144 -> 52 < N < 148-----------------> N = 128   |
(iii) N*N*3*4 <  31457280 -> 148 < N < 1619------------> N=1024    |BLOCKSIZE = 16
(iv) N*N*3*4 < 125829120 -> N < 3238  ----------------->N=2048     |
*/


long long unsigned initial,final;
struct timeval t;

float** A;
float** B;
float** C;
float** trans;
float** aux;
double clearcache[30000000];




void freeMatrices(){
	free(A);
	free(B);
	free(C);
	free(trans);
}

void printMatrix(float** matrix, int size) {
    int i, j;
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++)
            printf("%f ",matrix[i][j]); 
        printf("\n");
    }
    printf("\n");
}

void clearCache (void) {
	int i;
	for (i = 0; i < 30000000; i++)
		clearcache[i] = i;
}

void init_mat(int size){

	int i;

    A = malloc(size* sizeof(float *));
    B = malloc(size* sizeof(float *));
    C = malloc(size* sizeof(float *));
    trans = malloc(size* sizeof(float *));
    aux = malloc(size* sizeof(float *));

	for( i=0;i<size;i++){
		A[i] = malloc(size* sizeof(float));
		B[i] = malloc(size* sizeof(float));
		C[i] = malloc(size* sizeof(float));
		trans[i] = malloc(size* sizeof(float));
		aux[i] = malloc(size* sizeof(float));
	}
	
}

void fillMatrices (int size) {
	int i,j;

   
	srand(time(NULL));
   for ( i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			A[i][j] = ((float) rand()) / ((float) RAND_MAX)  * 100;
			B[i][j] = 1.0;
			C[i][j] = 0.0;
			aux[i][j] = 0.0;

			//	if(j==N-1) printf("\n");
			//	printf("%f ",A[i][j]);
				
		}
		
	}

}

void transpose(float** mat, int n){
	int i, j;
	float temp;

	for(i = 0; i < n; i++){
		for(j = 0; j < i; j++){
			temp = mat[i][j];
			mat[i][j] = mat[j][i];
			mat[j][i] = temp;
		}
	}
}

//1.2
void matrixMult(float **A, float **B,float **C, int SIZE ) {
	int i,j,k;
    //float x;

	for (i = 0; i < SIZE; i++) {
		for ( j = 0; j < SIZE; j++) {
			//x= 0;
			for (k = 0; k < SIZE; k++) {
				//x += A[i][k] * B[k][j];
				C[i][j] += A[i][k] * B[k][j];
			}
		   //C[i][j] += x;
		   
		
			//printf("%f\n",C[i][j]);
		}

	}
}


//1.2
void matrixMult_T(float **A, float **B,float **C, int SIZE ) {
	int i,j,k;
	//float x;

	transpose(B,SIZE);

	for (i = 0; i < SIZE; i++) {
		for ( j = 0; j < SIZE; j++) {
			//x = 0;
			for (k = 0; k < SIZE; k++) {
				//x += A[i][k] * B[j][k];
				C[i][j] += A[i][k] * B[j][k];

			}
			//C[i][j] += x;
		//	printf("%f\n",C[i][j]);
		}

	}
}


//1.3
void matrixMult_ikj(float **A, float **B, float **C, int SIZE) {
    int i, j, k;
    // float x=0;

    for (i = 0; i < SIZE; i++) {
        for (k = 0; k < SIZE; k++) {
        	//x = A[i][k];
            for (j = 0; j < SIZE; j++) {
               C[i][j] +=  A[i][k] * B[k][j];
            }
  
   	    }
	}
}

//1.3
void matrixMult_jki(float **A, float **B, float **C, int SIZE) {
    int i, j, k;
    //float x=0;
  

    for (j = 0; j < SIZE; j++){
        for (k = 0; k < SIZE; k++) {
        	//x = B[k][j];
            for (i = 0; i < SIZE; i++){
               C[i][j]+= A[i][k] * B[k][j];
              
            }
           // C[i][j] = x;
            //printf("%f\n",C[i][j]);
    	}
	}
	
}
	
//1.3
// MATRIZ A ACEDIDA POR COLUNA/B ACEDE POR COLUNA MAS SO DEPOIS DE SIZE ITERACOES,POIS É MULTIPLICADO SEMPRE A MESMA POSICAO
void matrixMult_jki_T(float **A, float **B, float **C, int SIZE) {
  
    int i, j, k;
    //float x;

    transpose(A,SIZE);
    transpose(B,SIZE);

    for (j = 0; j < SIZE; j++){
        for (k = 0; k < SIZE; k++) {
        	//x = B[j][k];
            for (i = 0; i < SIZE; i++){
               C[j][i] += A[i][k] * B[j][k];
            }
           
            // printf("%f\n",C[i][j]);
 	    }
	}
	transpose(C,SIZE);
}

//1.5

void start (void) {
	gettimeofday(&begin, NULL);
}


long long stop () {
	struct timeval end;
	gettimeofday(&end, NULL);
	long long duration = (end.tv_sec-begin.tv_sec)*1000000LL + end.tv_usec-begin.tv_usec;

	return duration;
}


long long sort_median(long long* a){
	int i, j, m;
	

	for(i = 0; i < 5; i++){
		for(j = 0; j < 5; j++){
			if(a[j] < a[i]){
				m = a[i];
				a[i] = a[j];
				a[j]  = m;
			}
		}
	}
	return a[2];
	
}

/*
1.6

i),ii) 
--matrixMult

RAM accesses per instruction
N=32 -> RAM / INS: 0.000729
N=128 ->RAM / INS: 0.000124
N=1024 ->RAM / INS: 0.000031
N=2048 ->RAM / INS: 0.000173

number of bytes transferred to/from the RAM
N=32 -> 12
N=128 ->134
N=1024 -> 17183
N=2048 ->759216

--matrixmult_T

RAM accesses per instruction
N=32 ->RAM / INS: 0.000953
N=128 ->RAM / INS: 0.000146
N=1024 ->RAM / INS: 0.000012
N=2048 ->RAM / INS: 0.000090

number of bytes transferred to/from the RAM
N=32 ->14
N=128 ->138
N=1024 ->6003
N=2048 ->347137

--matrixMult_ikj

RAM accesses per instruction
N=32 ->RAM / INS: 0.000977
N=128 ->RAM / INS: 0.000123
N=1024 ->RAM / INS: 0.000006
N=2048 ->RAM / INS: 0.000051

number of bytes transferred to/from the RAM
N=32 ->14
N=128 ->116
N=1024 ->3110
N=2048 ->195798

--matrixMult_jki

RAM accesses per instruction
N=32 ->RAM / INS: 0.000506
N=128 ->RAM / INS: 0.000115
N=1024 ->RAM / INS: 0.000029
N=2048 ->RAM / INS: 0.000155

number of bytes transferred to/from the RAM
N=32 ->10
N=128 ->155
N=1024 ->19631
N=2048 ->850758

--matrixMult_jki_T

RAM accesses per instruction
N=32 ->RAM / INS: 0.000684
N=128 ->RAM / INS: 0.000161
N=1024 -> RAM / INS: 0.000033
N=2048 ->RAM / INS: 0.000243

number of bytes transferred to/from the RAM
N=32 ->12
N=128 ->175
N=1024 ->17986
N=2048 ->1068746

*/

/*
1.7 
i)we will perform two floating point operations on each iteration (one multiplication and one accumulation) .
There is 3 nested cycles where each one iterates size times resulting in size 3 iterations.

N=32-----------------> FP = 2 * 32^3 =65536
N=128----------------->  FP = 2 * 128^3 =4194304
N=1024---------------->  FP = 2 * 1024^3 =2147483648
N=2048---------------->  FP = 2 * 2048^3 =17179869180

ii)

*/

/*
1.8
matrixMult

Cache L1
mr = 0,00329
mr= 0,02491
 mr = 0,60550
mr=  0,72696

CACHE L2
N= 32 miss rate: 0.452096  miss rate: 0.479495
N=128 miss rate: 0.026115  miss rate: 0.023977
N=1024 miss rate: 0.036637  miss rate: 0.037875
N=2048 miss rate: 0.519205  miss rate: 0.500934

CACHE L3

miss rate: 0.000256      miss rate: 0.989011
miss rate: 0.000416      miss rate: 0.548404
miss rate: 0.000107      miss rate: 0.002308
miss rate: 0.006320      miss rate: 0.001932

matrixMult_T

Cache L1
miss rate: 0.003786
miss rate: 0.042448
miss rate: 0.032730
miss rate: 0.032208

Cache L2
miss rate: 0.750929       miss rate: 0.891667
miss rate: 0.021624       miss rate: 0.038353
miss rate: 0.085100       miss rate: 0.086400
miss rate: 0.058554       miss rate: 0.060965

Cache L3
miss rate: 0.004484      miss rate: 0.970464
miss rate: 0.000426      miss rate: 0.213758
miss rate: 0.000853      miss rate: 0.014656
miss rate: 0.107686      miss rate: 0.089760


*/


//1.10

//i)
//Simple rearrangements of the loops can increase spatial locality, but the time per loop iteration increases with increasing array size
//What is happening is that as the array size increases, the xoral locality decreases, and the cache experiences an increasing number of capacity misses. 
//To fix this, we can use a general technique called blocking
//The general idea of blocking is to organize the data structures in a program into large chunks called blocks
//The program is structured so that it loads a chunk into the L1 cache, does all the reads and writes that it needs to on that chunk, then discards the chunk, loads in the next chunk, and so on.
//Blocking a matrix multiply routine works by partitioning the matrices into submatrices and then exploiting the mathematical fact that these submatrices can be manipulated just like scalars.
//This technique can produce big performance gains

//ii)
void matrixMult_block(float **A, float **B, float **C, int SIZE){
	int i, j, k, kk, jj;
	float x;
	int b_SIZE = 16;
	

	for(kk = 0; kk < SIZE; kk += b_SIZE){
		for(jj = 0; jj < SIZE; jj += b_SIZE){
			for(i = 0; i < SIZE; i++){
				for(j = jj; j < ((jj + b_SIZE) > SIZE ? SIZE : (jj + b_SIZE)); j++) {
				
					x=0;
					for(k = kk; k < ((kk + b_SIZE) > SIZE ? SIZE : (kk + b_SIZE)); k++) {
					
						x += A[i][k]*B[k][j];
						
					}
					C[i][j] += x;
				}
			}
		}
	}
}

//1.11

void matrixMult_block_Vec(float** A, float** B, float** C, int SIZE) {
	int i = 0, j = 0, jj = 0, kk = 0;
  float x[16];
  int b_SIZE = 16;


  for (jj = 0; jj < SIZE; jj += b_SIZE) {
    for (kk = 0; kk < SIZE; kk += b_SIZE) {
      for (i = 0; i < SIZE; i++) {       
        for (j = jj; j < ((jj + b_SIZE) > SIZE ? SIZE : (jj + b_SIZE)); j++) {
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

          
          C[i][j] += x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[10] + x[11] + x[12] + x[13] + x[14] + x[15];
        }
      }
    }
  }
}

//1.12
void matrixMult_blockVecOMP (float **A, float **B, float **C, int SIZE) {
  int i = 0, j = 0, jj = 0, kk = 0;
  float x[16];
  int b_SIZE = 16;

#pragma omp parallel for private(i, j, jj, kk, x)
  for (jj = 0; jj < SIZE; jj += b_SIZE) {
    for (kk = 0; kk < SIZE; kk += b_SIZE) {
      for (i = 0; i < SIZE; i++) {

        #pragma omp simd
        for (j = jj; j < ((jj + b_SIZE) > SIZE ? SIZE : (jj + b_SIZE)); j++) {
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

          
          C[i][j] += x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[10] + x[11] + x[12] + x[13] + x[14] + x[15];
        }
      }
    }
  }
}

/*
void compare() {
	int i, v = 1;
	for(i = 0;i < N; i++) {
		printf("RES: %f -> RESEQ: %f\n",aux[i][0],C[i][0]);
		if(trans[i] != C[i]) {
			v = 0;
			printf("ERRRRRROoooooooooooooooooooooooooooooooooooooooooooo\n");
		}
	}printf("DEU: %d\n",v);
}

*/

int main(int argc, char const *argv[]){

	//int i;
     int size;
	size = 32;
        //  size = atoi(argv[1]);
	//long long mul[5], mulT[5], mul_ikj[5],mul_jki[5],mul_jkiT[5], ret;

        long long mulB[5],mulBV[5],mulBVO[5];
	init_mat(size);
	fillMatrices(size);




	for(i = 0; i < 5; i++){
/*
		clearCache();
		start();
		matrixMult(A,B,C,size);
		mul[i] = stop();

		clearCache();
		start();
		matrixMult_T(A,B,C,size);
		mulT[i] = stop();

		clearCache();
		start();
		matrixMult_ikj(A,B,C,size);
		mul_ikj[i] = stop();


		clearCache();
		start();
		matrixMult_jki(A,B,C,size);
		mul_jki[i] = stop();


		clearCache();
		start();
		matrixMult_jki_T(A,B,C,size);
		mul_jkiT[i] = stop();
              */
             clearCache();
                start();
                matrixMult_block(A,B,C,size);
                mulB[i] = stop();

                clearCache();
                start();
                matrixMult_block_Vec(A,B,C,size);
                mulBV[i] = stop();


                clearCache();
                start();
                matrixMult_blockVecOMP(A,B,C,size);
                mulBVO[i] = stop();
                 
	}

/*
	ret = sort_median(mul);
	printf("Tempo Mult = %lld\n", ret);

	ret = sort_median(mulT);
	printf("Tempo Mult Trans = %lld\n", ret);

	ret = sort_median(mul_ikj);
	printf("Tempo ikj = %lld\n", ret);

	ret = sort_median(mul_jki);
	printf("Tempo jki = %lld\n", ret);

	ret = sort_median(mul_jkiT);
	printf("Tempo jkiT = %lld\n", ret);
*/
        ret = sort_median(mulB);
        printf("Tempo block = %lld\n", ret);

        ret = sort_median(mulBV);
        printf("Tempo vector = %lld\n", ret);

        ret = sort_median(mulBVO);
        printf("Tempo omp = %lld\n", ret);

	//PAPI_library_init(PAPI_VER_CURRENT);
	//PAPI_create_eventset(&EventSet);
	//PAPI_add_events(EventSet,Events,NUM_EVENTS);

	//start();

	//PAPI_start(EventSet);

    // printMatrix(C, size);
	//matrixMult(A,B,C,size);
	//matrixMult_T(A,B,C,size);
	//matrixMult_ikj(A,B,C,size);
	//matrixMult_jki(A,B,C,size);
	//matrixMult_jki_T(A,B,C,size);

	//matrixMult_block(A,B,aux,N,16);
	//compare();

    //PAPI_stop(EventSet,values);
    //printf("time: ");
	//stop();
      
	//printf("Cache L1->L2 DCR:%lld LD INS:%lld\n",values[0],values[1]);
        //double r = (double) values[0]/values[1];
        //printf("miss rate: %f\n",r);
        
     // printf("Cache L2->L3 DCR:%lld L2 DCR:%lld\n",values[0],values[1]);
     // double r = (double) values[0]/values[1];
      //printf("miss rate: %f\n",r);

        // printf("Cache L3-> L3 TCM:%lld L3 TCA:%lld\n",values[0],values[1]);
          //double r = (double) values[0]/values[1];
           //printf("miss rate: %f\n",r);

	//printf("RAM ACCESSES: %lld\n",values[0]);
	//printf("INSTR: %lld\n",values[1]);

	//double r = (double) values[0]/values[1];
	//printf("RAM / INS: %f\n",r);


	
}
	



	

	


