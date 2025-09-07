#include "charge.h"
#include <stdlib.h>

void effective_charge_energy(int nat, const double *q, const double *gamma,
                             double *energy){
    double e = 0.0;
    for(int i=0;i<nat;i++){
        double vi = 0.0;
        for(int j=0;j<nat;j++) vi += gamma[i*nat + j]*q[j];
        e += 0.5*q[i]*vi;
    }
    *energy += e;
}

void effective_charge_potential(int nat, const double *q, const double *gamma,
                                double *vat){
    for(int i=0;i<nat;i++){
        vat[i] = 0.0;
        for(int j=0;j<nat;j++) vat[i] += gamma[i*nat + j]*q[j];
    }
}
