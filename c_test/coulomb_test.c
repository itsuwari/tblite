#include "../c_src/xtb/coulomb.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

static double rand_double(void){ return (double)rand()/RAND_MAX*2.0-1.0; }

int main(){
    srand(0);
    int nat=3;
    tb_coulomb c;
    tb_coulomb_new(&c,nat);
    for(int i=0;i<nat*nat;i++) c.gamma[i]=rand_double();
    for(int i=0;i<3*nat*nat;i++) c.amat_sd[i]=rand_double();
    for(int i=0;i<3*nat*3*nat;i++) c.amat_dd[i]=rand_double();
    for(int i=0;i<6*nat*nat;i++) c.amat_sq[i]=rand_double();
    double *qat=malloc(sizeof(double)*nat);
    double *dpat=malloc(sizeof(double)*3*nat);
    double *qpat=malloc(sizeof(double)*6*nat);
    for(int i=0;i<nat;i++) qat[i]=rand_double();
    for(int i=0;i<3*nat;i++) dpat[i]=rand_double();
    for(int i=0;i<6*nat;i++) qpat[i]=rand_double();
    double e=0.0;
    tb_coulomb_energy(&c,qat,dpat,qpat,&e);
    double e_ref=0.0;
    effective_charge_energy(nat,qat,c.gamma,&e_ref);
    multipole_energy(nat,c.amat_sd,c.amat_dd,c.amat_sq,qat,dpat,qpat,&e_ref);
    assert(fabs(e-e_ref)<1e-12);
    double *vat=malloc(sizeof(double)*nat);
    double *vdp=malloc(sizeof(double)*3*nat);
    double *vqp=malloc(sizeof(double)*6*nat);
    tb_coulomb_potential(&c,qat,dpat,qpat,vat,vdp,vqp);
    double *vat_ref=malloc(sizeof(double)*nat);
    double *vdp_ref=malloc(sizeof(double)*3*nat);
    double *vqp_ref=malloc(sizeof(double)*6*nat);
    multipole_potential(nat,c.amat_sd,c.amat_dd,c.amat_sq,qat,dpat,qpat,
                        vdp_ref,vat_ref,vqp_ref);
    double *vat_iso=malloc(sizeof(double)*nat);
    effective_charge_potential(nat,qat,c.gamma,vat_iso);
    for(int i=0;i<nat;i++) vat_ref[i]+=vat_iso[i];
    for(int i=0;i<nat;i++) assert(fabs(vat[i]-vat_ref[i])<1e-12);
    for(int i=0;i<3*nat;i++) assert(fabs(vdp[i]-vdp_ref[i])<1e-12);
    for(int i=0;i<6*nat;i++) assert(fabs(vqp[i]-vqp_ref[i])<1e-12);
    free(vat);free(vdp);free(vqp);
    free(vat_ref);free(vdp_ref);free(vqp_ref);free(vat_iso);
    free(qat);free(dpat);free(qpat);
    tb_coulomb_free(&c);
    printf("coulomb container tests passed\n");
    return 0;
}
