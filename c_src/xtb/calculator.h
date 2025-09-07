#ifndef TBLITE_XTB_CALCULATOR_H
#define TBLITE_XTB_CALCULATOR_H

#include "coulomb.h"
#include "repulsion.h"
#include "singlepoint.h"

typedef struct {
    const structure_type *mol;
    const basis_type *bas;
    tb_coulomb *coul;
    tb_repulsion *rep;
} xtb_calculator;

void xtb_calculator_new(xtb_calculator *calc, const structure_type *mol,
                        const basis_type *bas, tb_coulomb *coul,
                        tb_repulsion *rep);

double xtb_calculator_energy(xtb_calculator *calc, wavefunction_type *wfn,
                             diag_solver_type *solver, mixer_type *mixer,
                             scf_info info, const integral_type *ints,
                             potential_type *pot, int maxiter);

#endif /* TBLITE_XTB_CALCULATOR_H */
