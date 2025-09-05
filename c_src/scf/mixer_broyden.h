#ifndef TBLITE_SCF_MIXER_BROYDEN_H
#define TBLITE_SCF_MIXER_BROYDEN_H

#include "mixer_type.h"

typedef struct {
    mixer_type base;
    int ndim;
    int memory;
    int iter;
    int iset;
    int idif;
    int iget;
    double damp;
    double *df;
    double *u;
    double *a;
    double *dq;
    double *dqlast;
    double *qlast_in;
    double *omega;
    double *q_in;
} broyden_mixer;

typedef struct {
    int memory;
    double damp;
} broyden_input;

void new_broyden(broyden_mixer *self, int ndim, broyden_input input);

#endif
