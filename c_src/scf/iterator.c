#include "iterator.h"
#include <math.h>
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
        int ish = 0;
        for(int iat=0; iat<bas->nat; ++iat){
            for(int s=0; s<bas->nsh_at[iat]; ++s){
                qat[iat + spin*bas->nat] += qsh[ish + spin*nsh];
                ish++;
            }
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

void set_mixer(mixer_type *mixer, const wavefunction_type *wfn, scf_info info,
               const basis_type *bas, int nspin){
    int nat = bas->nat;
    int nsh = bas->nsh;
    switch(info.charge){
        case scf_atom_resolved:
            mixer->set(mixer, wfn->qat, nat);
            break;
        case scf_shell_resolved:
            mixer->set(mixer, wfn->qsh, nsh);
            break;
    }
    if(info.dipole==scf_atom_resolved)
        mixer->set(mixer, wfn->dpat, 3*nat);
    if(info.quadrupole==scf_atom_resolved)
        mixer->set(mixer, wfn->qpat, 6*nat);
}

void get_mixer(mixer_type *mixer, const basis_type *bas, wavefunction_type *wfn,
               scf_info info, int nspin){
    int nat = bas->nat;
    int nsh = bas->nsh;
    switch(info.charge){
        case scf_atom_resolved:
            mixer->get(mixer, wfn->qat, nat);
            break;
        case scf_shell_resolved:
            mixer->get(mixer, wfn->qsh, nsh);
            get_qat_from_qsh(bas, wfn->qsh, wfn->qat, nsh, nspin);
            break;
    }
    if(info.dipole==scf_atom_resolved)
        mixer->get(mixer, wfn->dpat, 3*nat);
    if(info.quadrupole==scf_atom_resolved)
        mixer->get(mixer, wfn->qpat, 6*nat);
}

void diff_mixer(mixer_type *mixer, const wavefunction_type *wfn, scf_info info,
                const basis_type *bas, int nspin){
    int nat = bas->nat;
    int nsh = bas->nsh;
    switch(info.charge){
        case scf_atom_resolved:
            mixer->diff(mixer, wfn->qat, nat);
            break;
        case scf_shell_resolved:
            mixer->diff(mixer, wfn->qsh, nsh);
            break;
    }
    if(info.dipole==scf_atom_resolved)
        mixer->diff(mixer, wfn->dpat, 3*nat);
    if(info.quadrupole==scf_atom_resolved)
        mixer->diff(mixer, wfn->qpat, 6*nat);
}

static double get_entropy(const double *occ, double kt, int n){
    double s=0.0;
    for(int i=0;i<n;i++){
        double o=occ[i];
        if(o>0.0 && o<1.0)
            s += log(pow(o,o)*pow(1.0-o,1.0-o));
    }
    return s*kt;
}

void next_density(wavefunction_type *wfn, diag_solver_type *solver,
                  const integral_type *ints, int norb, int nspin, double *ts){
    get_density(solver, ints->overlap, wfn->coeff, wfn->emo, wfn->focc,
                wfn->density, norb, nspin);
    *ts = 0.0;
    for(int spin=0; spin<nspin; ++spin)
        *ts += get_entropy(&wfn->focc[spin*norb], wfn->kt, norb);
}

void next_scf(int *iscf, const structure_type *mol, const basis_type *bas,
              wavefunction_type *wfn, diag_solver_type *solver,
              mixer_type *mixer, scf_info info, const integral_type *ints,
              potential_type *pot, double *energies){
    int nspin = pot->nspin;
    int norb = pot->nao;
    if(*iscf > 0){
        mixer->next(mixer);
        get_mixer(mixer, bas, wfn, info, nspin);
    }
    (*iscf)++;
    reset_potential(pot);
    add_pot_to_h1(bas, ints, pot, wfn->coeff);
    set_mixer(mixer, wfn, info, bas, nspin);
    double ts=0.0;
    next_density(wfn, solver, ints, norb, nspin, &ts);
    get_qat_from_qsh(bas, wfn->qsh, wfn->qat, bas->nsh, nspin);
    diff_mixer(mixer, wfn, info, bas, nspin);
    double *eao = calloc(norb, sizeof(double));
    get_electronic_energy(ints->hamiltonian, wfn->density, eao, norb, nspin);
    for(int i=0;i<mol->nat;i++) energies[i] = ts / mol->nat;
    reduce(energies, eao, bas->ao2at, norb);
    free(eao);
}
