#include "calculator.h"
#include <stdlib.h>

void xtb_calculator_new(xtb_calculator *calc, const structure_type *mol,
                        const basis_type *bas, tb_coulomb *coul,
                        tb_repulsion *rep){
    calc->mol = mol;
    calc->bas = bas;
    calc->coul = coul;
    calc->rep  = rep;
}

double xtb_calculator_energy(xtb_calculator *calc, wavefunction_type *wfn,
                             diag_solver_type *solver, mixer_type *mixer,
                             scf_info info, const integral_type *ints,
                             potential_type *pot, int maxiter){
    double energy = 0.0;
    xtb_singlepoint(calc->mol, calc->bas, wfn, solver, mixer, info, ints, pot,
                    maxiter, &energy);
    tb_coulomb_energy(calc->coul, wfn->qat, wfn->dpat, wfn->qpat, &energy);
    if (calc->rep) {
        double *eat = malloc(sizeof(double) * calc->rep->nat);
        tb_repulsion_energy(calc->rep, eat, &energy);
        free(eat);
    }
    return energy;
}
