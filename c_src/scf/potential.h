#ifndef TBLITE_SCF_POTENTIAL_H
#define TBLITE_SCF_POTENTIAL_H

#include <stdlib.h>

/* Basic structure representing a molecule */
typedef struct {
    int nat; /* number of atoms */
} structure_type;

/* Basis set description */
typedef struct {
    int nat;  /* number of atoms */
    int nsh;  /* number of shells */
    int nao;  /* number of atomic orbitals */
    int *ish_at;   /* starting shell index for each atom (0-based) */
    int *nsh_at;   /* number of shells per atom */
    int *iao_sh;   /* starting AO index for each shell (0-based) */
    int *nao_sh;   /* number of AOs per shell */
    int *ao2at;    /* mapping from AO to atom */
} basis_type;

/* Integral container */
typedef struct {
    double *hamiltonian; /* size nao*nao */
    double *overlap;     /* size nao*nao */
    double *dipole;      /* size 3*nao*nao */
    double *quadrupole;  /* size 6*nao*nao */
} integral_type;

/* Container for density dependent potentials */
typedef struct {
    int grad;          /* gradient flag */
    int nat, nsh, nao, nspin;
    double *vat;       /* size nat*nspin */
    double *vsh;       /* size nsh*nspin */
    double *vao;       /* size nao*nspin */
    double *vdp;       /* size 3*nat*nspin */
    double *vqp;       /* size 6*nat*nspin */
    double *dvatdr;    /* size 3*nat*nat*nspin */
    double *dvatdL;    /* size 3*3*nat*nspin */
    double *dvshdr;    /* size 3*nat*nsh*nspin */
    double *dvshdL;    /* size 3*3*nsh*nspin */
} potential_type;

void new_potential(potential_type *self, const structure_type *mol,
                   const basis_type *bas, int nspin, int grad);
void reset_potential(potential_type *self);
void free_potential(potential_type *self);

void add_pot_to_h1(const basis_type *bas, const integral_type *ints,
                   potential_type *pot, double *h1);

#endif /* TBLITE_SCF_POTENTIAL_H */
