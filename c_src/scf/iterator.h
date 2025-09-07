#ifndef TBLITE_SCF_ITERATOR_H
#define TBLITE_SCF_ITERATOR_H

#include "info.h"
#include "mixer_type.h"
#include "diag.h"
#include "potential.h"

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

void get_electronic_energy(const double *h0, const double *density, double *energies,
                           int norb, int nspin);
void reduce(double *reduced, const double *full, const int *map, int nmap);
void get_qat_from_qsh(const basis_type *bas, const double *qsh, double *qat, int nsh, int nspin);
int get_mixer_dimension(const structure_type *mol, const basis_type *bas, scf_info info);

void set_mixer(mixer_type *mixer, const wavefunction_type *wfn, scf_info info,
               const basis_type *bas, int nspin);
void get_mixer(mixer_type *mixer, const basis_type *bas, wavefunction_type *wfn,
               scf_info info, int nspin);
void diff_mixer(mixer_type *mixer, const wavefunction_type *wfn, scf_info info,
                const basis_type *bas, int nspin);
void next_density(wavefunction_type *wfn, diag_solver_type *solver,
                  const integral_type *ints, int norb, int nspin, double *ts);
void next_scf(int *iscf, const structure_type *mol, const basis_type *bas,
              wavefunction_type *wfn, diag_solver_type *solver,
              mixer_type *mixer, scf_info info, const integral_type *ints,
              potential_type *pot, double *energies);

#endif
