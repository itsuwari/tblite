#include "coulomb.h"
#include <stdlib.h>
#include <string.h>

void tb_coulomb_new(tb_coulomb *c, int nat){
    c->nat = nat;
    c->gamma = calloc((size_t)nat*nat, sizeof(double));
    c->amat_sd = calloc((size_t)3*nat*nat, sizeof(double));
    c->amat_dd = calloc((size_t)3*nat*3*nat, sizeof(double));
    c->amat_sq = calloc((size_t)6*nat*nat, sizeof(double));
}

void tb_coulomb_free(tb_coulomb *c){
    free(c->gamma);
    free(c->amat_sd);
    free(c->amat_dd);
    free(c->amat_sq);
    c->gamma=c->amat_sd=c->amat_dd=c->amat_sq=NULL;
    c->nat=0;
}

void tb_coulomb_energy(const tb_coulomb *c, const double *qat,
                       const double *dpat, const double *qpat,
                       double *energy){
    effective_charge_energy(c->nat, qat, c->gamma, energy);
    multipole_energy(c->nat, c->amat_sd, c->amat_dd, c->amat_sq,
                     qat, dpat, qpat, energy);
}

void tb_coulomb_potential(const tb_coulomb *c, const double *qat,
                          const double *dpat, const double *qpat,
                          double *vat, double *vdp, double *vqp){
    multipole_potential(c->nat, c->amat_sd, c->amat_dd, c->amat_sq,
                        qat, dpat, qpat, vdp, vat, vqp);
    double *vat_iso = malloc(sizeof(double)*c->nat);
    effective_charge_potential(c->nat, qat, c->gamma, vat_iso);
    for(int i=0;i<c->nat;i++) vat[i] += vat_iso[i];
    free(vat_iso);
}
