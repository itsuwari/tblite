#include "../c_src/coulomb/charge.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

static double rand_double(void){ return (double)rand()/RAND_MAX*2.0-1.0; }

int main(){
    srand(0);
    int nat=4;
    double *q=malloc(sizeof(double)*nat);
    double *gamma=malloc(sizeof(double)*nat*nat);
    double *vat=malloc(sizeof(double)*nat);
    for(int i=0;i<nat;i++){
        q[i]=rand_double();
        for(int j=0;j<nat;j++) gamma[i*nat+j]=rand_double();
    }
    double e=0.0;
    effective_charge_energy(nat,q,gamma,&e);
    effective_charge_potential(nat,q,gamma,vat);
    double e_ref=0.0;
    for(int i=0;i<nat;i++){
        double vi=0.0;
        for(int j=0;j<nat;j++) vi+=gamma[i*nat+j]*q[j];
        e_ref+=0.5*q[i]*vi;
        assert(fabs(vi-vat[i])<1e-12);
    }
    assert(fabs(e-e_ref)<1e-12);
    free(q);free(gamma);free(vat);
    printf("charge tests passed\n");
    return 0;
}
