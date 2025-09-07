#include "../c_src/xtb/singlepoint.h"
#include "../c_src/scf/mixer.h"
#include "../c_src/scf/mixer_broyden.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define IDX2(i,j,n1) ((j)*(n1)+(i))

static double rand_double(void){ return (double)rand()/RAND_MAX*2.0-1.0; }

void identity_solve(diag_solver_type *self, double *hmat, const double *smat,
                    double *eval, int norb){ (void)self;(void)hmat;(void)smat;(void)eval;(void)norb; }

int main(){
    srand(0);
    int nat=2,nsh=3,nao=4,nspin=2;
    structure_type mol={nat};
    basis_type bas; bas.nat=nat; bas.nsh=nsh; bas.nao=nao;
    int ish_at[2]={0,2}; int nsh_at[2]={2,1};
    int iao_sh[3]={0,1,3}; int nao_sh[3]={1,2,1};
    int ao2at[4]={0,0,1,1};
    bas.ish_at=ish_at; bas.nsh_at=nsh_at; bas.iao_sh=iao_sh; bas.nao_sh=nao_sh; bas.ao2at=ao2at;
    integral_type ints; ints.hamiltonian=malloc(sizeof(double)*nao*nao);
    ints.overlap=malloc(sizeof(double)*nao*nao);
    ints.dipole=calloc(3*nao*nao,sizeof(double));
    ints.quadrupole=calloc(6*nao*nao,sizeof(double));
    for(int i=0;i<nao;i++) for(int j=0;j<nao;j++){
        ints.hamiltonian[IDX2(j,i,nao)]=(i==j)?1.0:0.0;
        ints.overlap[IDX2(j,i,nao)]=(i==j)?1.0:0.0;
    }
    size_t matsz=(size_t)nao*nao*nspin;
    size_t vecsz=(size_t)nao*nspin;
    wavefunction_type wfn_ref;
    wfn_ref.coeff=calloc(matsz,sizeof(double));
    wfn_ref.emo=malloc(sizeof(double)*vecsz);
    wfn_ref.focc=calloc(vecsz,sizeof(double));
    wfn_ref.density=calloc(matsz,sizeof(double));
    wfn_ref.n0sh=calloc((size_t)nsh*nspin,sizeof(double));
    wfn_ref.qsh=malloc(sizeof(double)*(size_t)nsh*nspin);
    wfn_ref.qat=malloc(sizeof(double)*(size_t)nat*nspin);
    wfn_ref.dpat=calloc((size_t)3*nat*nspin,sizeof(double));
    wfn_ref.qpat=calloc((size_t)6*nat*nspin,sizeof(double));
    wfn_ref.kt=0.1;
    for(size_t i=0;i<vecsz;i++) wfn_ref.emo[i]=rand_double();
    for(size_t i=0;i<(size_t)nsh*nspin;i++) wfn_ref.qsh[i]=rand_double();
    for(size_t i=0;i<(size_t)nat*nspin;i++) wfn_ref.qat[i]=rand_double();
    wavefunction_type wfn;
    wfn.coeff=malloc(sizeof(double)*matsz); memcpy(wfn.coeff,wfn_ref.coeff,sizeof(double)*matsz);
    wfn.emo=malloc(sizeof(double)*vecsz); memcpy(wfn.emo,wfn_ref.emo,sizeof(double)*vecsz);
    wfn.focc=calloc(vecsz,sizeof(double));
    wfn.density=calloc(matsz,sizeof(double));
    wfn.n0sh=calloc((size_t)nsh*nspin,sizeof(double));
    wfn.qsh=malloc(sizeof(double)*(size_t)nsh*nspin); memcpy(wfn.qsh,wfn_ref.qsh,sizeof(double)*(size_t)nsh*nspin);
    wfn.qat=malloc(sizeof(double)*(size_t)nat*nspin); memcpy(wfn.qat,wfn_ref.qat,sizeof(double)*(size_t)nat*nspin);
    wfn.dpat=calloc((size_t)3*nat*nspin,sizeof(double));
    wfn.qpat=calloc((size_t)6*nat*nspin,sizeof(double));
    wfn.kt = wfn_ref.kt;
    potential_type pot, pot_ref;
    new_potential(&pot,&mol,&bas,nspin,0);
    new_potential(&pot_ref,&mol,&bas,nspin,0);
    double nel[2]={2.0,2.0};
    diag_solver_type solver; new_diag_solver(&solver,wfn.kt,nel,identity_solve);
    scf_info info={scf_atom_resolved,scf_not_used,scf_not_used};
    int ndim=get_mixer_dimension(&mol,&bas,info);
    mixer_type *mix1,*mix2; new_mixer(&mix1,2,ndim,0.5); new_mixer(&mix2,2,ndim,0.5);
    int niter=2;
    double energies[3]={0};
    int iscf=0;
    for(int iter=0; iter<niter; ++iter){
        next_scf(&iscf,&mol,&bas,&wfn_ref,&solver,mix2,info,&ints,&pot_ref,energies);
    }
    double e_ref=energies[0];
    double e=0.0;
    xtb_singlepoint(&mol,&bas,&wfn,&solver,mix1,info,&ints,&pot,niter,&e);
    assert(fabs(e-e_ref)<1e-12);
    for(size_t i=0;i<matsz;i++) assert(fabs(wfn.coeff[i]-wfn_ref.coeff[i])<1e-12);
    for(size_t i=0;i<vecsz;i++) assert(fabs(wfn.focc[i]-wfn_ref.focc[i])<1e-12);
    for(size_t i=0;i<matsz;i++) assert(fabs(wfn.density[i]-wfn_ref.density[i])<1e-12);
    for(size_t i=0;i<(size_t)nat*nspin;i++) assert(fabs(wfn.qat[i]-wfn_ref.qat[i])<1e-12);
    mix1->del(mix1); mix2->del(mix2); free(mix1); free(mix2);
    free_potential(&pot); free_potential(&pot_ref);
    delete_diag_solver(&solver);
    free(ints.hamiltonian); free(ints.overlap); free(ints.dipole); free(ints.quadrupole);
    free(wfn.coeff); free(wfn.emo); free(wfn.focc); free(wfn.density); free(wfn.n0sh); free(wfn.qsh); free(wfn.qat); free(wfn.dpat); free(wfn.qpat);
    free(wfn_ref.coeff); free(wfn_ref.emo); free(wfn_ref.focc); free(wfn_ref.density); free(wfn_ref.n0sh); free(wfn_ref.qsh); free(wfn_ref.qat); free(wfn_ref.dpat); free(wfn_ref.qpat);
    printf("singlepoint tests passed\n");
    return 0;
}
