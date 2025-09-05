#include "info.h"

scf_info max_info(scf_info lhs, scf_info rhs){
    scf_info out;
    out.charge = lhs.charge > rhs.charge ? lhs.charge : rhs.charge;
    out.dipole = lhs.dipole > rhs.dipole ? lhs.dipole : rhs.dipole;
    out.quadrupole = lhs.quadrupole > rhs.quadrupole ? lhs.quadrupole : rhs.quadrupole;
    return out;
}
