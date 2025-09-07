#include "../c_src/scf/iterator.h"
#include "../c_src/scf/diag.h"
#include "../c_src/scf/mixer.h"
#include "../c_src/scf/mixer_broyden.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define IDX2(i,j,n1) ((j)*(n1)+(i))
#define IDX3(i,j,k,n1,n2) ((k)*(n1)*(n2)+(j)*(n1)+(i))

double rand_double(void){
    return (double)rand()/RAND_MAX*2.0 - 1.0;
}

void identity_solve(diag_solver_type *self, double *hmat, const double *smat,
                    double *eval, int norb){
    (void)self; (void)hmat; (void)smat; (void)eval; (void)norb;
}

int main(){
    srand(0);
    int nat=2, nsh=3, nao=4, nspin=2;
    structure_type mol = {nat};
    basis_type bas; bas.nat=nat; bas.nsh=nsh; bas.nao=nao;
    int ish_at_arr[2]={0,2};
    int nsh_at_arr[2]={2,1};
    int iao_sh_arr[3]={0,1,3};
    int nao_sh_arr[3]={1,2,1};
    int ao2at_arr[4]={0,0,1,1};
    bas.ish_at=ish_at_arr;
    bas.nsh_at=nsh_at_arr;
    bas.iao_sh=iao_sh_arr;
    bas.nao_sh=nao_sh_arr;
    bas.ao2at=ao2at_arr;

    integral_type ints;
    ints.hamiltonian=malloc(sizeof(double)*nao*nao);
    ints.overlap=malloc(sizeof(double)*nao*nao);
    ints.dipole=calloc(3*nao*nao,sizeof(double));
    ints.quadrupole=calloc(6*nao*nao,sizeof(double));
    for(int i=0;i<nao;i++){
        for(int j=0;j<nao;j++){
            ints.hamiltonian[IDX2(j,i,nao)] = (i==j)?1.0:0.0;
            ints.overlap[IDX2(j,i,nao)] = (i==j)?1.0:0.0;
        }
    }

    wavefunction_type wfn, wfn_ref;
    size_t matsz = (size_t)nao*nao*nspin;
    size_t vecsz = (size_t)nao*nspin;
    wfn.coeff = malloc(sizeof(double)*matsz);
    wfn.emo   = malloc(sizeof(double)*vecsz);
    wfn.focc  = calloc(vecsz, sizeof(double));
    wfn.density = calloc(matsz, sizeof(double));
    wfn.n0sh = calloc((size_t)nsh*nspin, sizeof(double));
    wfn.qsh  = malloc(sizeof(double)*(size_t)nsh*nspin);
    wfn.qat  = malloc(sizeof(double)*(size_t)nat*nspin);
    wfn.dpat = calloc((size_t)3*nat*nspin, sizeof(double));
    wfn.qpat = calloc((size_t)6*nat*nspin, sizeof(double));
    wfn.kt = 0.1;

    for(size_t i=0;i<matsz;i++) wfn.coeff[i]=0.0;
    for(size_t i=0;i<vecsz;i++) wfn.emo[i]=rand_double();
    for(size_t i=0;i<(size_t)nsh*nspin;i++) wfn.qsh[i]=rand_double();
    for(size_t i=0;i<(size_t)nat*nspin;i++) wfn.qat[i]=rand_double();

    wfn_ref = wfn;
    wfn_ref.coeff = malloc(sizeof(double)*matsz);
    memcpy(wfn_ref.coeff, wfn.coeff, sizeof(double)*matsz);
    wfn_ref.emo = malloc(sizeof(double)*vecsz); memcpy(wfn_ref.emo, wfn.emo, sizeof(double)*vecsz);
    wfn_ref.focc = malloc(sizeof(double)*vecsz); memset(wfn_ref.focc,0,sizeof(double)*vecsz);
    wfn_ref.density = malloc(sizeof(double)*matsz); memset(wfn_ref.density,0,sizeof(double)*matsz);
    wfn_ref.n0sh = calloc((size_t)nsh*nspin,sizeof(double));
    wfn_ref.qsh = malloc(sizeof(double)*(size_t)nsh*nspin); memcpy(wfn_ref.qsh,wfn.qsh,sizeof(double)*(size_t)nsh*nspin);
    wfn_ref.qat = malloc(sizeof(double)*(size_t)nat*nspin); memcpy(wfn_ref.qat,wfn.qat,sizeof(double)*(size_t)nat*nspin);
    wfn_ref.dpat = calloc((size_t)3*nat*nspin,sizeof(double));
    wfn_ref.qpat = calloc((size_t)6*nat*nspin,sizeof(double));
    wfn_ref.kt = wfn.kt;

    potential_type pot, pot_ref;
    new_potential(&pot, &mol, &bas, nspin, 0);
    new_potential(&pot_ref, &mol, &bas, nspin, 0);

    double nel[2]={2.0,2.0};
    diag_solver_type solver; new_diag_solver(&solver, wfn.kt, nel, identity_solve);

    scf_info info = {scf_atom_resolved, scf_not_used, scf_not_used};
    int ndim = get_mixer_dimension(&mol,&bas,info);
    mixer_type *mix1, *mix2;
    new_mixer(&mix1,2,ndim,0.5);
    new_mixer(&mix2,2,ndim,0.5);

    double energies[2]={0,0};
    double energies_ref[2]={0,0};

    int iscf1=0;
    iscf1++;
    reset_potential(&pot_ref);
    add_pot_to_h1(&bas, &ints, &pot_ref, wfn_ref.coeff);
    set_mixer(mix1, &wfn_ref, info, &bas, nspin);
    double ts_ref=0.0;
    next_density(&wfn_ref, &solver, &ints, bas.nao, nspin, &ts_ref);
    get_qat_from_qsh(&bas, wfn_ref.qsh, wfn_ref.qat, bas.nsh, nspin);
    diff_mixer(mix1, &wfn_ref, info, &bas, nspin);
    double *eao = calloc(nao, sizeof(double));
    get_electronic_energy(ints.hamiltonian, wfn_ref.density, eao, nao, nspin);
    for(int i=0;i<nat;i++) energies_ref[i]=ts_ref/nat;
    reduce(energies_ref, eao, bas.ao2at, nao);
    free(eao);

    int iscf2=0;
    next_scf(&iscf2, &mol, &bas, &wfn, &solver, mix2, info, &ints, &pot, energies);

    assert(iscf1==iscf2);
    for(size_t i=0;i<matsz;i++) assert(fabs(wfn.coeff[i]-wfn_ref.coeff[i])<1e-12);
    for(size_t i=0;i<vecsz;i++) assert(fabs(wfn.focc[i]-wfn_ref.focc[i])<1e-12);
    for(size_t i=0;i<matsz;i++) assert(fabs(wfn.density[i]-wfn_ref.density[i])<1e-12);
    for(size_t i=0;i<(size_t)nat*nspin;i++) assert(fabs(wfn.qat[i]-wfn_ref.qat[i])<1e-12);
    for(int i=0;i<nat;i++) assert(fabs(energies[i]-energies_ref[i])<1e-12);

    broyden_mixer *b1=(broyden_mixer*)mix1;
    broyden_mixer *b2=(broyden_mixer*)mix2;
    for(int i=0;i<ndim;i++){
        assert(fabs(b1->q_in[i]-b2->q_in[i])<1e-12);
        assert(fabs(b1->dq[i]-b2->dq[i])<1e-12);
    }

    mix1->del(mix1); mix2->del(mix2);
    free(mix1); free(mix2);
    free_potential(&pot); free_potential(&pot_ref);
    delete_diag_solver(&solver);
    free(ints.hamiltonian); free(ints.overlap); free(ints.dipole); free(ints.quadrupole);
    free(wfn.coeff); free(wfn.emo); free(wfn.focc); free(wfn.density);
    free(wfn.n0sh); free(wfn.qsh); free(wfn.qat); free(wfn.dpat); free(wfn.qpat);
    free(wfn_ref.coeff); free(wfn_ref.emo); free(wfn_ref.focc); free(wfn_ref.density);
    free(wfn_ref.n0sh); free(wfn_ref.qsh); free(wfn_ref.qat); free(wfn_ref.dpat); free(wfn_ref.qpat);
    printf("solver tests passed\n");
    return 0;
}
