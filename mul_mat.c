#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <immintrin.h>
//#include <papi.h>
#include <unistd.h>

//#define NUM_EVENTS 2
struct timeval begin;

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


float** A;
float** B;
float** C;
float** trans;
float** aux;
double clearcache[30000000];

//int Events[NUM_EVENTS] = {PAPI_L3_TCM,PAPI_TOT_INS};

//int EventSet = PAPI_NULL;
//long long values[NUM_EVENTS];


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

void start (void) {
	gettimeofday(&begin, NULL);
}


void stop () {
	struct timeval end;
	gettimeofday(&end, NULL);
	long long duration = (end.tv_sec-begin.tv_sec)*1000000LL + end.tv_usec-begin.tv_usec;
	printf("%lld\n", duration);
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

void transpose(float** source, int SIZE){

	
	int i, j;

	for (i = 0; i < SIZE; i++) {
		
   	 for (j = 0; j < SIZE; j++) {
    	trans[j][i] = source[i][j];
	}
  }

    for(i = 0; i < SIZE; i++) {
        for(j = 0; j < SIZE; j++) 
            source[i][j] = trans[i][j];          
	}	      
}

//1.2
void matrixMult(float **A, float **B,float **C, int SIZE ) {
	int i,j,k;
    float x;

	for (i = 0; i < SIZE; i++) {
		for ( j = 0; j < SIZE; j++) {
			x= 0;
			for (k = 0; k < SIZE; k++) {
				x += A[i][k] * B[k][j];
				//C[i][j] += A[i][k] * B[k][j];
			}
		   C[i][j] += x;
		   
		
			printf("%f\n",C[i][j]);
		}

	}
}


//1.2
void matrixMult_T(float **A, float **B,float **C, int SIZE ) {
	int i,j,k;
	float x;

	transpose(B,SIZE);

	for (i = 0; i < SIZE; i++) {
		for ( j = 0; j < SIZE; j++) {
			x = 0;
			for (k = 0; k < SIZE; k++) {
				x += A[i][k] * B[k][j];
				//C[i][j] += A[i][k] * B[k][j];
			}
			C[i][j] += x;
		//	printf("%f\n",C[i][j]);
		}

	}
}


//1.3
void matrixMult_ikj(float **A, float **B, float **C, int SIZE) {
    int i, j, k;
     float x=0;

    for (i = 0; i < SIZE; i++) {
        for (k = 0; k < SIZE; k++) {
        	x = A[i][k];
            for (j = 0; j < SIZE; j++) {
               C[i][j] += x * B[k][j];
            }
  
   	    }
	}
}

//1.3
void matrixMult_jki(float **A, float **B, float **C, int SIZE) {
    int i, j, k;
    float x=0;
  

    for (j = 0; j < SIZE; j++){
        for (k = 0; k < SIZE; k++) {
        	x = B[k][j];
            for (i = 0; i < SIZE; i++){
               C[i][j]+= A[i][k] * x;
              
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
    float x;

    transpose(A,SIZE);
    transpose(B,SIZE);

    for (j = 0; j < SIZE; j++){
        for (k = 0; k < SIZE; k++) {
        	x = B[k][j];
            for (i = 0; i < SIZE; i++){
               C[i][j] += A[i][k] * x;
            }
           
            // printf("%f\n",C[i][j]);
 	    }
	}
}

//1.5
double double_time() {
    double time_sec = 0.0;
    struct timeval time;
    gettimeofday(&time,(struct timezone*)0);
    time_sec = (double)(time.tv_sec + time.tv_usec*1.0e-6);
    return( time_sec );
}

static int compare_double (const void * a, const void * b)
{
	 if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;

}

int Best_k(double* array,int size, int k, double tol) {
    int i;
    double sigma, top;
    if(size < 3) return 0;
    qsort(array, size, sizeof(double), compare_double);
    sigma = (double)(array[0] * tol);
    top = (double)(array[0] + sigma);
    for(i=1; i<k; i++) {
        if(array[i] > top) return 0;
    }
    return 1;
}

void function_run(void (*function)(float**, float**, float **, int), int size, double* t) {
    int i, t_size=0;
    double t0, t1;
    for(i = 0; i < 8; i++){    
        clearCache();
        t0 = double_time();
        (*function)(A,B,C,size);
        t1 = (double_time() - t0);
        t[t_size] = (double) (t1 * 1000.0);
       printf("%f\n", t[t_size]);
        t_size++;
      
        Best_k(t,t_size,3,0.05);
    }
    printf("%f %f %f \n",t[0],t[1],t[2]);
}

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

void matrixMult_block_Vec(float** __restrict__ A, float** __restrict__ B, float** __restrict__ C, int SIZE, int b_SIZE) {
	int i, j, k, kk, jj;
	float x;
	

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
	}
	printf("DEU: %d\n",v);
}

*/

int main(){
	double *t = malloc(sizeof(double)*8);

	int size = 3;
    
	init_mat(size);
	fillMatrices(size);

	//PAPI_library_init(PAPI_VER_CURRENT);
	//PAPI_create_eventset(&EventSet);
	//PAPI_add_events(EventSet,Events,NUM_EVENTS);

	//clearCache();

	//start();

	//PAPI_start(EventSet);

	//function_run(matrixMult,size,t);
    //function_run(matrixMult_T,size,t);
    function_run(matrixMult_ikj,size,t);
    
    //function_run(matrixMult_jki,size,t);
    //function_run(matrixMult_jki_T,size,t);

     printMatrix(A, size);
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

	//printf("RAM ACCESSES: %lld\n",values[0]);
	//printf("INSTR: %lld\n",values[1]);

	//double r = (double) values[0]/values[1];
	//printf("RAM / INS: %f\n",r);


	
}
	



	

	


