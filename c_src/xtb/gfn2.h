#ifndef TBLITE_XTB_GFN2_H
#define TBLITE_XTB_GFN2_H

#include <stddef.h>
#include "h0.h"

/* Hubbard parameters for each element (Z=1..86) */
extern const double gfn2_hubbard_parameter[86];

/* Shell Hubbard derivative weights (l=0..4) */
extern const double gfn2_shell_hubbard_derivs[5];

/* Element-specific Hubbard derivatives (Z=1..86) */
extern const double gfn2_p_hubbard_derivs[86];

/* Retrieve shell-specific Hubbard scaling factor (l=0..2, Z=1..86) */
double gfn2_shell_hubbard_value(int l, int z);

/* Compute shell-resolved hardness values.
 * - natoms: number of atoms
 * - nums: array of atomic numbers length natoms (1-based Z)
 * - nsh_id: array with number of shells per atom length natoms
 * - cgto_ang: array of size (stride*natoms) giving angular momentum for
 *             each shell and atom; stride must be >= max(nsh_id).
 * - stride: leading dimension for shell-resolved arrays
 * - hardness: output array of size (stride*natoms) filled with values
 */
void gfn2_shell_hardness(size_t natoms, const int *nums, const int *nsh_id,
                         const int *cgto_ang, size_t stride, double *hardness);

/* Compute shell-resolved Hubbard derivatives */
void gfn2_hubbard_derivs(size_t natoms, const int *nums, const int *nsh_id,
                          const int *cgto_ang, size_t stride, double *derivs);

/* number of shells per element */
extern const int gfn2_nshell[86];

/* Accessors for Hamiltonian parameters */
double gfn2_selfenergy_value(int l, int z);
double gfn2_kcn_value(int l, int z);
double gfn2_refocc_value(int l, int z);

/* Fill Hamiltonian parameter container for given element list */
void gfn2_init_h0(tb_hamiltonian *h0, const int *zlist);

#endif /* TBLITE_XTB_GFN2_H */
