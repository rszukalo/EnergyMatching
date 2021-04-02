
#ifndef EMATCH_FUNC_H
#define EMATCH_FUNC_H

#include <stdlib.h>
#include <stdio.h>

#include "cgff_types.h"
#include "gromacs_topology.h"
#include "safe_mem.h"
#include "io_read.h"
#include "LDD.h"
#include "ematch_types.h"


extern void dgelss_(int *M, int *N, int *NRHS, double A[], int *LDA,
                    double B[], int *LDB, double S[], double *RCOND,
                    int *RANK, double WORK[], int *LWORK, int *INFO);

void test_stuff();

void summarize_esystem(tW_esystem *sys);

void setup_ematch_sys(tW_esystem *sys);

void do_ematch_bonds(tW_esystem *sys);

void do_ematch_angles(tW_esystem *sys);

void do_ematch_dihedrals(tW_esystem *sys);

void do_ematch_pairs(tW_esystem *sys);

void do_ematch_volume(tW_esystem *sys);

void save_A_matrix(double **A, int nrows, int ncols, char * fnm);

void read_A_matrix(double **A, int nrows, int ncols);

void trim_A(tW_esystem *sys);

void get_answer(tW_esystem *sys);

void save_answer(tW_esystem *sys);

void check_answer(tW_esystem *sys);

#endif
