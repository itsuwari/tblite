#ifndef TBLITE_SCF_INFO_H
#define TBLITE_SCF_INFO_H

typedef struct {
    int charge;
    int dipole;
    int quadrupole;
} scf_info;

enum {
    scf_not_used = 0,
    scf_atom_resolved = 1,
    scf_shell_resolved = 2,
    scf_orbital_resolved = 3
};

scf_info max_info(scf_info lhs, scf_info rhs);

#endif
