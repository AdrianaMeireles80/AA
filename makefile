SHELL = /bin/sh

BIN_NAME = mul_mat

CC = gcc

PAPI = -lpapi

FLAGS = -O2 -L/share/apps/papi/5.5.0/lib -I/share/apps/papi/5.5.0/include -Wall -Wextra -fopenmp -Wno-unused-parameter
VEC = -O3 -L/share/apps/papi/5.5.0/lib -I/share/apps/papi/5.5.0/include -Wall -Wextra -fopenmp -Wno-unused-parameter -vec-report3

compile: mul_mat.c
	$(CC) -o $(BIN_NAME) mul_mat.c $(FLAGS) $(PAPI)