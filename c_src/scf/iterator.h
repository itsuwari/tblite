#ifndef TBLITE_SCF_ITERATOR_H
#define TBLITE_SCF_ITERATOR_H

#include "info.h"
#include "mixer_type.h"
#include "solver.h"

typedef struct {
    int nat;
} structure_type;

typedef struct {
    int nsh;
    int nao;
    int nat;
    int *ao2at;
    int *sh2at;
} basis_type;

typedef struct {
    double *coeff;
    double *emo;
    double *focc;
    double *density;
    double *n0sh;
    double *qsh;
    double *qat;
    double *dpat;
    double *qpat;
    double kt;
} wavefunction_type;

typedef struct {
    double *overlap;
    double *hamiltonian;
    double *dipole;
    double *quadrupole;
} integral_type;

void get_electronic_energy(const double *h0, const double *density, double *energies,
                           int norb, int nspin);
void reduce(double *reduced, const double *full, const int *map, int nmap);
void get_qat_from_qsh(const basis_type *bas, const double *qsh, double *qat, int nsh, int nspin);
int get_mixer_dimension(const structure_type *mol, const basis_type *bas, scf_info info);

#endif
