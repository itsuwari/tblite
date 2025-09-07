#include "disp/d4.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double frand(void) { return rand() / (double)RAND_MAX; }

int main(void) {
    srand(1);
    int nat = 3;
    int zmax = 3;
    int mref = 2;
    int id[3];
    double xyz[9];
    for (int i = 0; i < nat; ++i) {
        id[i] = 1 + rand() % zmax;
        for (int k = 0; k < 3; ++k)
            xyz[3*i+k] = 2.0*frand() - 1.0;
    }
    int ref[3];
    for (int z = 0; z < zmax; ++z) ref[z] = mref;
    double r4r2[3];
    for (int z = 0; z < zmax; ++z) r4r2[z] = frand() + 0.5;
    double c6[mref*mref*zmax*zmax];
    for (size_t i = 0; i < sizeof(c6)/sizeof(c6[0]); ++i)
        c6[i] = frand();
    double trans[3] = {0.0, 0.0, 0.0};
    struct d4_param param = {frand(), frand(), frand(), frand()};
    double cutoff = 10.0;

    size_t sz = (size_t)mref*nat*mref*nat;
    double *dispmat1 = calloc(sz, sizeof(double));
    double *dispmat2 = calloc(sz, sizeof(double));

    d4_get_dispersion_matrix(nat, mref, zmax, ref, id, xyz,
                             &param, r4r2, c6, trans, 1, cutoff, dispmat1);

    double cutoff2 = cutoff*cutoff;
    for (int iat = 0; iat < nat; ++iat) {
        int izp = id[iat]-1;
        for (int jat = 0; jat <= iat; ++jat) {
            int jzp = id[jat]-1;
            double rrij = 3.0 * r4r2[izp] * r4r2[jzp];
            double r0ij = param.a1 * sqrt(rrij) + param.a2;
            double r0_2 = r0ij*r0ij;
            double r0_3 = r0_2*r0ij;
            double r0_6 = r0_3*r0_3;
            double r0_8 = r0_6*r0_2;
            double dE = 0.0;
            double vecx = xyz[3*iat+0] - (xyz[3*jat+0] + trans[0]);
            double vecy = xyz[3*iat+1] - (xyz[3*jat+1] + trans[1]);
            double vecz = xyz[3*iat+2] - (xyz[3*jat+2] + trans[2]);
            double r2 = vecx*vecx + vecy*vecy + vecz*vecz;
            if (r2 <= cutoff2 && r2 > 1e-12) {
                double r2_2 = r2*r2;
                double r2_3 = r2_2*r2;
                double t6 = 1.0/(r2_3 + r0_6);
                double t8 = 1.0/(r2_2*r2_2 + r0_8);
                double edisp = param.s6*t6 + param.s8*rrij*t8;
                dE -= edisp;
            }
            for (int iref = 0; iref < ref[izp]; ++iref) {
                for (int jref = 0; jref < ref[jzp]; ++jref) {
                    size_t cidx1 = (((size_t)iref*mref + jref)*zmax + izp)*zmax + jzp;
                    size_t cidx2 = (((size_t)jref*mref + iref)*zmax + jzp)*zmax + izp;
                    double val1 = dE * c6[cidx1];
                    double val2 = dE * c6[cidx2];
                    size_t idx1 = (((size_t)iref*nat + iat)*mref + jref)*nat + jat;
                    size_t idx2 = (((size_t)jref*nat + jat)*mref + iref)*nat + iat;
                    dispmat2[idx1] = val1;
                    dispmat2[idx2] = val2;
                }
            }
        }
    }

    int ok = 1;
    for (size_t i = 0; i < sz; ++i) {
        if (fabs(dispmat1[i] - dispmat2[i]) > 1e-12) {
            ok = 0;
            break;
        }
    }
    if (!ok) {
        fprintf(stderr, "d4 dispersion matrix mismatch\n");
        free(dispmat1);
        free(dispmat2);
        return 1;
    }

    size_t n = (size_t)mref*nat;
    double *gw = malloc(n * sizeof(double));
    double *dgwdq = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
        gw[i] = frand();
        dgwdq[i] = frand();
    }

    double *eat1 = calloc(nat, sizeof(double));
    double *eat2 = calloc(nat, sizeof(double));
    d4_self_energy(nat, mref, dispmat1, gw, eat1);
    for (int i = 0; i < nat; ++i) eat2[i] = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j)
            sum += dispmat1[i*n + j] * gw[j];
        sum *= 0.5;
        int at = i % nat;
        eat2[at] += sum * gw[i];
    }
    for (int i = 0; i < nat; ++i) {
        if (fabs(eat1[i] - eat2[i]) > 1e-12) {
            fprintf(stderr, "d4 self energy mismatch\n");
            free(dispmat1); free(dispmat2); free(gw); free(dgwdq); free(eat1); free(eat2);
            return 1;
        }
    }

    double *vat1 = calloc(nat, sizeof(double));
    double *vat2 = calloc(nat, sizeof(double));
    d4_self_potential(nat, mref, dispmat1, gw, dgwdq, vat1);
    for (int i = 0; i < nat; ++i) vat2[i] = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j)
            sum += dispmat1[i*n + j] * gw[j];
        int at = i % nat;
        vat2[at] += sum * dgwdq[i];
    }
    for (int i = 0; i < nat; ++i) {
        if (fabs(vat1[i] - vat2[i]) > 1e-12) {
            fprintf(stderr, "d4 self potential mismatch\n");
            free(dispmat1); free(dispmat2); free(gw); free(dgwdq); free(eat1); free(eat2); free(vat1); free(vat2);
            return 1;
        }
    }

    free(dispmat1);
    free(dispmat2);
    free(gw);
    free(dgwdq);
    free(eat1);
    free(eat2);
    free(vat1);
    free(vat2);
    return 0;
}
