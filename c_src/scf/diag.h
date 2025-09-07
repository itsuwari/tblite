#ifndef TBLITE_SCF_DIAG_H
#define TBLITE_SCF_DIAG_H

#include <stddef.h>

typedef struct diag_solver_type diag_solver_type;

typedef void (*solve_fn)(diag_solver_type *self, double *hmat, const double *smat,
                         double *eval, int norb);

struct diag_solver_type {
    double kt;
    double nel[2];
    solve_fn solve;
};

void new_diag_solver(diag_solver_type *s, double kt, const double nel[2], solve_fn solve);
void delete_diag_solver(diag_solver_type *s);

void get_density(diag_solver_type *self, const double *smat, double *hmat,
                 double *eval, double *focc, double *density,
                 int norb, int nspin);

void get_wdensity(diag_solver_type *self, const double *smat, double *hmat,
                  double *eval, double *focc, double *density,
                  int norb, int nspin);

void get_density_matrix(const double *focc, const double *coeff,
                        double *pmat, int n);

void get_fermi_filling(double nel, double kt, const double *emo, int n,
                       int *homo, double *focc, double *efermi);

#endif
