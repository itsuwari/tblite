#include "iterator.h"
#include <stdlib.h>
#include <string.h>

#define IDX2(i,j,n1) ((j)*(n1)+(i))
#define IDX3(i,j,k,n1,n2) ((k)*(n1)*(n2)+(j)*(n1)+(i))

void get_electronic_energy(const double *h0, const double *density, double *energies,
                           int norb, int nspin){
    for(int iao=0; iao<norb; ++iao) energies[iao]=0.0;
    for(int spin=0; spin<nspin; ++spin){
        for(int iao=0; iao<norb; ++iao){
            for(int jao=0; jao<norb; ++jao){
                energies[iao] += h0[IDX2(jao, iao, norb)] * density[IDX3(jao, iao, spin, norb, norb)];
            }
        }
    }
}

void reduce(double *reduced, const double *full, const int *map, int nmap){
    for(int i=0;i<nmap;i++) reduced[map[i]] += full[i];
}

void get_qat_from_qsh(const basis_type *bas, const double *qsh, double *qat, int nsh, int nspin){
    for(int i=0;i<bas->nat*nspin;i++) qat[i]=0.0;
    for(int spin=0; spin<nspin; ++spin){
        for(int ish=0; ish<nsh; ++ish){
            int at = bas->sh2at[ish];
            qat[at + spin*bas->nat] += qsh[ish + spin*nsh];
        }
    }
}

int get_mixer_dimension(const structure_type *mol, const basis_type *bas, scf_info info){
    int ndim=0;
    switch(info.charge){
        case scf_atom_resolved: ndim += mol->nat; break;
        case scf_shell_resolved: ndim += bas->nsh; break;
    }
    if(info.dipole==scf_atom_resolved) ndim += 3*mol->nat;
    if(info.quadrupole==scf_atom_resolved) ndim += 6*mol->nat;
    return ndim;
}
