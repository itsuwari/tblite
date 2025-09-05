#ifndef TBLITE_XTB_GFN1_H
#define TBLITE_XTB_GFN1_H

#include <stddef.h>

/* Hubbard parameters for each element (Z=1..86) */
extern const double gfn1_hubbard_parameter[86];

/* Retrieve shell-specific Hubbard scaling factor (l=0..2, Z=1..86) */
double gfn1_shell_hubbard_value(int l, int z);

/* Compute shell-resolved hardness values.
 * - natoms: number of atoms
 * - nums: array of atomic numbers length natoms (1-based Z)
 * - nsh_id: array with number of shells per atom length natoms
 * - cgto_ang: array of size (stride*natoms) giving angular momentum for
 *             each shell and atom; stride must be >= max(nsh_id).
 * - stride: leading dimension for shell-resolved arrays
 * - hardness: output array of size (stride*natoms) filled with values
 */
void gfn1_shell_hardness(size_t natoms, const int *nums, const int *nsh_id,
                         const int *cgto_ang, size_t stride, double *hardness);

#endif /* TBLITE_XTB_GFN1_H */
