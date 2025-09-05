#ifndef TBLITE_SCF_SOLVER_H
#define TBLITE_SCF_SOLVER_H

typedef struct solver_type solver_type;

typedef void (*get_density_fn)(solver_type *self, const double *smat, double *hmat,
                                double *eval, double *focc, double *density,
                                int norb, int nspin);
typedef void (*get_wdensity_fn)(solver_type *self, const double *smat, double *hmat,
                                double *eval, double *focc, double *density,
                                int norb, int nspin);
typedef void (*delete_solver_fn)(solver_type *self);

struct solver_type {
    double kt;
    double nel[2];
    get_density_fn get_density;
    get_wdensity_fn get_wdensity;
    delete_solver_fn del;
};

#endif
