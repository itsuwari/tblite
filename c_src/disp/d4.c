#include "disp/d4.h"
#include <math.h>
#include <string.h>

static size_t idx4(int iref, int iat, int jref, int jat, int mref, int nat) {
    return ((size_t)iref*nat + iat)*mref*nat + (size_t)jref*nat + jat;
}

void d4_get_dispersion_matrix(int nat, int mref, int zmax,
                              const int *ref, const int *id,
                              const double *xyz, const struct d4_param *param,
                              const double *r4r2, const double *c6,
                              const double *trans, int ntrans,
                              double cutoff, double *dispmat)
{
    size_t total = (size_t)mref * nat * mref * nat;
    for (size_t i = 0; i < total; ++i) dispmat[i] = 0.0;
    double cutoff2 = cutoff * cutoff;

    for (int iat = 0; iat < nat; ++iat) {
        int izp = id[iat] - 1; /* assume id are 1-based */
        for (int jat = 0; jat <= iat; ++jat) {
            int jzp = id[jat] - 1;
            double rrij = 3.0 * r4r2[izp] * r4r2[jzp];
            double r0ij = param->a1 * sqrt(rrij) + param->a2;
            double r0_2 = r0ij * r0ij;
            double r0_3 = r0_2 * r0ij;
            double r0_6 = r0_3 * r0_3;
            double r0_8 = r0_6 * r0_2;
            double dE = 0.0;
            for (int jtr = 0; jtr < ntrans; ++jtr) {
                double vecx = xyz[3*iat+0] - (xyz[3*jat+0] + trans[3*jtr+0]);
                double vecy = xyz[3*iat+1] - (xyz[3*jat+1] + trans[3*jtr+1]);
                double vecz = xyz[3*iat+2] - (xyz[3*jat+2] + trans[3*jtr+2]);
                double r2 = vecx*vecx + vecy*vecy + vecz*vecz;
                if (r2 > cutoff2 || r2 < 1e-12) continue;
                double r2_2 = r2 * r2;
                double r2_3 = r2_2 * r2;
                double t6 = 1.0 / (r2_3 + r0_6);
                double t8 = 1.0 / (r2_2 * r2_2 + r0_8);
                double edisp = param->s6 * t6 + param->s8 * rrij * t8;
                dE -= edisp;
            }
            for (int iref = 0; iref < ref[izp]; ++iref) {
                for (int jref = 0; jref < ref[jzp]; ++jref) {
                    size_t cidx1 = (((size_t)iref*mref + jref)*zmax + izp)*zmax + jzp;
                    size_t cidx2 = (((size_t)jref*mref + iref)*zmax + jzp)*zmax + izp;
                    double val1 = dE * c6[cidx1];
                    double val2 = dE * c6[cidx2];
                    dispmat[idx4(iref, iat, jref, jat, mref, nat)] = val1;
                    dispmat[idx4(jref, jat, iref, iat, mref, nat)] = val2;
                }
            }
        }
    }
}

void d4_self_energy(int nat, int mref, const double *dispmat,
                    const double *gw, double *eat)
{
    int n = mref * nat;
    for (int i = 0; i < nat; ++i) eat[i] = 0.0;
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        const double *row = dispmat + (size_t)i * n;
        for (int j = 0; j < n; ++j) sum += row[j] * gw[j];
        sum *= 0.5;
        int at = i % nat;
        eat[at] += sum * gw[i];
    }
}

void d4_self_potential(int nat, int mref, const double *dispmat,
                       const double *gw, const double *dgwdq,
                       double *vat)
{
    int n = mref * nat;
    for (int i = 0; i < nat; ++i) vat[i] = 0.0;
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        const double *row = dispmat + (size_t)i * n;
        for (int j = 0; j < n; ++j) sum += row[j] * gw[j];
        int at = i % nat;
        vat[at] += sum * dgwdq[i];
    }
}
