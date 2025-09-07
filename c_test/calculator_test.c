#include "../c_src/xtb/calculator.h"
#include "../c_src/xtb/repulsion.h"
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
    wavefunction_type w0;
    w0.coeff=calloc(matsz,sizeof(double));
    w0.emo=malloc(sizeof(double)*vecsz);
    w0.focc=calloc(vecsz,sizeof(double));
    w0.density=calloc(matsz,sizeof(double));
    w0.n0sh=calloc((size_t)nsh*nspin,sizeof(double));
    w0.qsh=malloc(sizeof(double)*(size_t)nsh*nspin);
    w0.qat=malloc(sizeof(double)*(size_t)nat*nspin);
    w0.dpat=calloc((size_t)3*nat*nspin,sizeof(double));
    w0.qpat=calloc((size_t)6*nat*nspin,sizeof(double));
    w0.kt=0.1;
    for(size_t i=0;i<vecsz;i++) w0.emo[i]=rand_double();
    for(size_t i=0;i<(size_t)nsh*nspin;i++) w0.qsh[i]=rand_double();
    for(size_t i=0;i<(size_t)nat*nspin;i++) w0.qat[i]=rand_double();
    wavefunction_type w_manual=w0;
    w_manual.coeff=malloc(sizeof(double)*matsz); memcpy(w_manual.coeff,w0.coeff,sizeof(double)*matsz);
    w_manual.emo=malloc(sizeof(double)*vecsz); memcpy(w_manual.emo,w0.emo,sizeof(double)*vecsz);
    w_manual.focc=calloc(vecsz,sizeof(double));
    w_manual.density=calloc(matsz,sizeof(double));
    w_manual.n0sh=calloc((size_t)nsh*nspin,sizeof(double));
    w_manual.qsh=malloc(sizeof(double)*(size_t)nsh*nspin); memcpy(w_manual.qsh,w0.qsh,sizeof(double)*(size_t)nsh*nspin);
    w_manual.qat=malloc(sizeof(double)*(size_t)nat*nspin); memcpy(w_manual.qat,w0.qat,sizeof(double)*(size_t)nat*nspin);
    w_manual.dpat=calloc((size_t)3*nat*nspin,sizeof(double));
    w_manual.qpat=calloc((size_t)6*nat*nspin,sizeof(double));
    wavefunction_type w_calc=w_manual;
    w_calc.coeff=malloc(sizeof(double)*matsz); memcpy(w_calc.coeff,w0.coeff,sizeof(double)*matsz);
    w_calc.emo=malloc(sizeof(double)*vecsz); memcpy(w_calc.emo,w0.emo,sizeof(double)*vecsz);
    w_calc.focc=calloc(vecsz,sizeof(double));
    w_calc.density=calloc(matsz,sizeof(double));
    w_calc.n0sh=calloc((size_t)nsh*nspin,sizeof(double));
    w_calc.qsh=malloc(sizeof(double)*(size_t)nsh*nspin); memcpy(w_calc.qsh,w0.qsh,sizeof(double)*(size_t)nsh*nspin);
    w_calc.qat=malloc(sizeof(double)*(size_t)nat*nspin); memcpy(w_calc.qat,w0.qat,sizeof(double)*(size_t)nat*nspin);
    w_calc.dpat=calloc((size_t)3*nat*nspin,sizeof(double));
    w_calc.qpat=calloc((size_t)6*nat*nspin,sizeof(double));
    potential_type pot1,pot2;
    new_potential(&pot1,&mol,&bas,nspin,0);
    new_potential(&pot2,&mol,&bas,nspin,0);
    double nel[2]={2.0,2.0};
    diag_solver_type solver; new_diag_solver(&solver,w0.kt,nel,identity_solve);
    scf_info info={scf_atom_resolved,scf_not_used,scf_not_used};
    int ndim=get_mixer_dimension(&mol,&bas,info);
    mixer_type *mix1,*mix2; new_mixer(&mix1,2,ndim,0.5); new_mixer(&mix2,2,ndim,0.5);
    tb_coulomb coul; tb_coulomb_new(&coul,nat);
    for(int i=0;i<nat*nat;i++) coul.gamma[i]=rand_double();
    int nid=3; int id[2]={0,1}; double xyz[6]={0,0,0, 1.0,0,0};
    tb_repulsion rep; tb_repulsion_new(&rep,nat,id,xyz,nid);
    for(int i=0;i<nid*nid;i++){
        rep.alpha[i]=fabs(rand_double());
        rep.zeff[i]=rand_double();
        rep.kexp[i]=fabs(rand_double())+1.0;
        rep.rexp[i]=fabs(rand_double())+1.0;
    }
    xtb_calculator calc; xtb_calculator_new(&calc,&mol,&bas,&coul,&rep);
    int niter=2;
    double e_manual=0.0; xtb_singlepoint(&mol,&bas,&w_manual,&solver,mix1,info,&ints,&pot1,niter,&e_manual);
    tb_coulomb_energy(&coul,w_manual.qat,w_manual.dpat,w_manual.qpat,&e_manual);
    double *eat = malloc(sizeof(double)*nat);
    tb_repulsion_energy(&rep,eat,&e_manual);
    free(eat);
    double e_calc = xtb_calculator_energy(&calc,&w_calc,&solver,mix2,info,&ints,&pot2,niter);
    assert(fabs(e_manual-e_calc)<1e-12);
    tb_coulomb_free(&coul);
    tb_repulsion_free(&rep);
    mix1->del(mix1); mix2->del(mix2); free(mix1); free(mix2);
    free_potential(&pot1); free_potential(&pot2);
    delete_diag_solver(&solver);
    free(ints.hamiltonian); free(ints.overlap); free(ints.dipole); free(ints.quadrupole);
    free(w0.coeff); free(w0.emo); free(w0.focc); free(w0.density); free(w0.n0sh); free(w0.qsh); free(w0.qat); free(w0.dpat); free(w0.qpat);
    free(w_manual.coeff); free(w_manual.emo); free(w_manual.focc); free(w_manual.density); free(w_manual.n0sh); free(w_manual.qsh); free(w_manual.qat); free(w_manual.dpat); free(w_manual.qpat);
    free(w_calc.coeff); free(w_calc.emo); free(w_calc.focc); free(w_calc.density); free(w_calc.n0sh); free(w_calc.qsh); free(w_calc.qat); free(w_calc.dpat); free(w_calc.qpat);
    printf("calculator tests passed\n");
    return 0;
}
