#ifndef TBLITE_XTB_SINGLEPOINT_H
#define TBLITE_XTB_SINGLEPOINT_H

#include "../scf/iterator.h"

void xtb_singlepoint(const structure_type *mol, const basis_type *bas,
                     wavefunction_type *wfn, diag_solver_type *solver,
                     mixer_type *mixer, scf_info info, const integral_type *ints,
                     potential_type *pot, int maxiter, double *energy);

#endif /* TBLITE_XTB_SINGLEPOINT_H */
