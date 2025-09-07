#include "singlepoint.h"

void xtb_singlepoint(const structure_type *mol, const basis_type *bas,
                     wavefunction_type *wfn, diag_solver_type *solver,
                     mixer_type *mixer, scf_info info, const integral_type *ints,
                     potential_type *pot, int maxiter, double *energy){
    double energies[3]={0.0,0.0,0.0};
    int iscf=0;
    for(int iter=0; iter<maxiter; ++iter){
        next_scf(&iscf, mol, bas, wfn, solver, mixer, info, ints, pot, energies);
    }
    *energy = energies[0];
}
