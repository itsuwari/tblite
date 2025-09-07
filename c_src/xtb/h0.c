#include "h0.h"
#include <stdlib.h>

void new_hamiltonian(tb_hamiltonian *h0, int mshell, int nid) {
    h0->mshell = mshell;
    h0->nid = nid;
    size_t sz = (size_t)mshell * nid;
    h0->selfenergy = (double*)calloc(sz, sizeof(double));
    h0->kcn = (double*)calloc(sz, sizeof(double));
    h0->kq1 = (double*)calloc(sz, sizeof(double));
    h0->kq2 = (double*)calloc(sz, sizeof(double));
    h0->refocc = (double*)calloc(sz, sizeof(double));
}

void free_hamiltonian(tb_hamiltonian *h0) {
    free(h0->selfenergy);
    free(h0->kcn);
    free(h0->kq1);
    free(h0->kq2);
    free(h0->refocc);
}

void get_selfenergy(const tb_hamiltonian *h0, int nat,
                    const int *id, const int *ish_at,
                    const int *nshell, const double *cn,
                    const double *qat, double *selfenergy,
                    double *dsedcn, double *dsedq) {
    for (int i=0; i<nat; ++i) {
        int izp = id[i];
        int ii = ish_at[i];
        for (int ish=0; ish<nshell[izp-1]; ++ish) {
            int idx = ish + h0->mshell*(izp-1);
            selfenergy[ii+ish] = h0->selfenergy[idx];
            if (dsedcn) dsedcn[ii+ish] = 0.0;
            if (dsedq) dsedq[ii+ish] = 0.0;
        }
    }
    if (cn) {
        if (dsedcn) {
            for (int i=0; i<nat; ++i) {
                int izp=id[i]; int ii=ish_at[i];
                for (int ish=0; ish<nshell[izp-1]; ++ish) {
                    int idx = ish + h0->mshell*(izp-1);
                    selfenergy[ii+ish] -= h0->kcn[idx]*cn[i];
                    dsedcn[ii+ish] = -h0->kcn[idx];
                }
            }
        } else {
            for (int i=0; i<nat; ++i) {
                int izp=id[i]; int ii=ish_at[i];
                for (int ish=0; ish<nshell[izp-1]; ++ish) {
                    int idx = ish + h0->mshell*(izp-1);
                    selfenergy[ii+ish] -= h0->kcn[idx]*cn[i];
                }
            }
        }
    }
    if (qat) {
        if (dsedq) {
            for (int i=0; i<nat; ++i) {
                int izp=id[i]; int ii=ish_at[i];
                double q=qat[i];
                for (int ish=0; ish<nshell[izp-1]; ++ish) {
                    int idx = ish + h0->mshell*(izp-1);
                    selfenergy[ii+ish] -= h0->kq1[idx]*q + h0->kq2[idx]*q*q;
                    dsedq[ii+ish] = -h0->kq1[idx] - 2*h0->kq2[idx]*q;
                }
            }
        } else {
            for (int i=0; i<nat; ++i) {
                int izp=id[i]; int ii=ish_at[i];
                double q=qat[i];
                for (int ish=0; ish<nshell[izp-1]; ++ish) {
                    int idx = ish + h0->mshell*(izp-1);
                    selfenergy[ii+ish] -= h0->kq1[idx]*q + h0->kq2[idx]*q*q;
                }
            }
        }
    }
}

void get_occupation(const tb_hamiltonian *h0, int nat, int nsh_tot,
                    const int *id, const int *ish_at,
                    const int *nshell, double charge,
                    double *nocc, double *n0at, double *n0sh) {
    *nocc = -charge;
    for (int i=0;i<nat;i++) n0at[i]=0.0;
    for (int i=0;i<nsh_tot;i++) n0sh[i]=0.0;
    for (int iat=0; iat<nat; ++iat) {
        int izp=id[iat];
        int ii=ish_at[iat];
        for (int ish=0; ish<nshell[izp-1]; ++ish) {
            int idx = ish + h0->mshell*(izp-1);
            double ref = h0->refocc[idx];
            *nocc += ref;
            n0at[iat] += ref;
            n0sh[ii+ish] += ref;
        }
    }
}

void shift_operator(const double vec[3], double s,
                    const double di[3], const double qi[6],
                    double dj[3], double qj[6]) {
    for (int k=0;k<3;k++) dj[k] = di[k] + vec[k]*s;
    qj[0] = 2*vec[0]*di[0] + vec[0]*vec[0]*s;
    qj[2] = 2*vec[1]*di[1] + vec[1]*vec[1]*s;
    qj[5] = 2*vec[2]*di[2] + vec[2]*vec[2]*s;
    qj[1] = vec[0]*di[1] + vec[1]*di[0] + vec[0]*vec[1]*s;
    qj[3] = vec[0]*di[2] + vec[2]*di[0] + vec[0]*vec[2]*s;
    qj[4] = vec[1]*di[2] + vec[2]*di[1] + vec[1]*vec[2]*s;
    double tr = 0.5*(qj[0] + qj[2] + qj[5]);
    qj[0] = qi[0] + 1.5*qj[0] - tr;
    qj[1] = qi[1] + 1.5*qj[1];
    qj[2] = qi[2] + 1.5*qj[2] - tr;
    qj[3] = qi[3] + 1.5*qj[3];
    qj[4] = qi[4] + 1.5*qj[4];
    qj[5] = qi[5] + 1.5*qj[5] - tr;
}
