#include "repulsion.h"
#include <stdlib.h>

void tb_repulsion_new(tb_repulsion *r, int nat, const int *id,
                      const double *xyz, int nid){
    r->nat = nat;
    r->id = id;
    r->xyz = xyz;
    r->nid = nid;
    r->cutoff = 10.0; /* reasonable default */
    size_t sz = (size_t)nid * nid;
    r->alpha = calloc(sz, sizeof(double));
    r->zeff  = calloc(sz, sizeof(double));
    r->kexp  = calloc(sz, sizeof(double));
    r->rexp  = calloc(sz, sizeof(double));
}

void tb_repulsion_free(tb_repulsion *r){
    free(r->alpha); free(r->zeff); free(r->kexp); free(r->rexp);
    r->alpha = r->zeff = r->kexp = r->rexp = NULL;
    r->nat = r->nid = 0;
    r->id = NULL; r->xyz = NULL;
}

void tb_repulsion_energy(const tb_repulsion *r, double *eat, double *energy){
    effective_repulsion_energy(r->nat, r->id, r->xyz, r->nid,
                               r->alpha, r->zeff, r->kexp, r->rexp,
                               r->cutoff, eat);
    for(int i=0;i<r->nat;i++) *energy += eat[i];
}
