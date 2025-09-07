#include "multipole.h"
#include <math.h>
#include <stdlib.h>

static inline double norm3(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void multipole_mrad(int nat, const int *id, const double *cn,
                    const double *rad, const double *vcn,
                    double shift, double kexp, double rmax,
                    double *mrad, double *dmrdcn) {
    for (int i=0;i<nat;i++) {
        int izp = id[i]-1;
        double arg = cn[i] - vcn[izp] - shift;
        double t1 = exp(-kexp*arg);
        double t2 = (rmax - rad[izp]) / (1.0 + t1);
        mrad[i] = rad[izp] + t2;
        dmrdcn[i] = -t2 * kexp * t1 / (1.0 + t1);
    }
}

#define SD(k,j,i) amat_sd[(k)*nat*nat + (j)*nat + (i)]
#define DD(k,j,m,i) amat_dd[((k)*nat + (j))*3*nat + (m)*nat + (i)]
#define SQ(c,j,i) amat_sq[(c)*nat*nat + (j)*nat + (i)]

void multipole_matrix_0d(int nat, const double *xyz, const double *mrad,
                         double kdmp3, double kdmp5,
                         double *amat_sd, double *amat_dd, double *amat_sq) {
    for (int k=0; k<3*nat*nat; ++k) amat_sd[k]=0.0;
    for (int k=0; k<3*nat*3*nat; ++k) amat_dd[k]=0.0;
    for (int k=0; k<6*nat*nat; ++k) amat_sq[k]=0.0;
    for (int i=0;i<nat;i++) {
        for (int j=0;j<nat;j++) {
            if (i==j) continue;
            double vec[3];
            vec[0]=xyz[3*i+0]-xyz[3*j+0];
            vec[1]=xyz[3*i+1]-xyz[3*j+1];
            vec[2]=xyz[3*i+2]-xyz[3*j+2];
            double r1 = norm3(vec);
            double g1 = 1.0/r1;
            double g3 = g1*g1*g1;
            double g5 = g3*g1*g1;
            double rr = 0.5*(mrad[j]+mrad[i])*g1;
            double fdmp3 = 1.0/(1.0 + 6.0*pow(rr,kdmp3));
            double fdmp5 = 1.0/(1.0 + 6.0*pow(rr,kdmp5));
            for (int k=0;k<3;k++) {
                SD(k,j,i) += vec[k]*g3*fdmp3;
                for (int m=0;m<3;m++) {
                    double unity = (k==m)?1.0:0.0;
                    DD(k,j,m,i) += unity*g3*fdmp5 - vec[k]*vec[m]*3.0*g5*fdmp5;
                }
            }
            double tc[6];
            tc[0] = vec[0]*vec[0]*g5*fdmp5;
            tc[1] = 2*vec[0]*vec[1]*g5*fdmp5;
            tc[2] = vec[1]*vec[1]*g5*fdmp5;
            tc[3] = 2*vec[0]*vec[2]*g5*fdmp5;
            tc[4] = 2*vec[1]*vec[2]*g5*fdmp5;
            tc[5] = vec[2]*vec[2]*g5*fdmp5;
            for (int c=0;c<6;c++) {
                SQ(c,j,i) += tc[c];
            }
        }
    }
}

void multipole_energy(int nat, const double *amat_sd, const double *amat_dd,
                      const double *amat_sq, const double *qat,
                      const double *dpat, const double *qpat, double *energy) {
    double *vd = (double*)calloc(3*nat, sizeof(double));
    double *vq = (double*)calloc(6*nat, sizeof(double));
    for (int i=0;i<nat;i++) {
        for (int k=0;k<3;k++) {
            for (int j=0;j<nat;j++) vd[k*nat+i] += SD(k,i,j) * qat[j];
            for (int m=0;m<3;m++) {
                for (int j=0;j<nat;j++) vd[k*nat+i] += 0.5*DD(k,i,m,j)*dpat[m*nat+j];
            }
        }
        for (int c=0;c<6;c++) {
            for (int j=0;j<nat;j++) vq[c*nat+i] += SQ(c,i,j) * qat[j];
        }
    }
    double e = 0.0;
    for (int i=0;i<nat;i++) {
        for (int k=0;k<3;k++) e += dpat[k*nat+i]*vd[k*nat+i];
        for (int c=0;c<6;c++) e += qpat[c*nat+i]*vq[c*nat+i];
    }
    *energy += e;
    free(vd); free(vq);
}

void multipole_potential(int nat, const double *amat_sd, const double *amat_dd,
                         const double *amat_sq, const double *qat,
                         const double *dpat, const double *qpat,
                         double *vdp, double *vat, double *vqp) {
    for (int k=0;k<3*nat;k++) vdp[k]=0.0;
    for (int i=0;i<nat;i++) vat[i]=0.0;
    for (int k=0;k<6*nat;k++) vqp[k]=0.0;
    for (int i=0;i<nat;i++) {
        for (int k=0;k<3;k++) {
            for (int j=0;j<nat;j++) vdp[k*nat+i] += SD(k,i,j)*qat[j];
            for (int m=0;m<3;m++) {
                for (int j=0;j<nat;j++) vdp[k*nat+i] += DD(k,i,m,j)*dpat[m*nat+j];
            }
        }
        for (int c=0;c<6;c++) {
            for (int j=0;j<nat;j++) vqp[c*nat+i] += SQ(c,i,j)*qat[j];
        }
    }
    for (int i=0;i<nat;i++) {
        for (int k=0;k<3;k++) {
            for (int j=0;j<nat;j++) vat[i] += SD(k,j,i)*dpat[k*nat+j];
        }
        for (int c=0;c<6;c++) {
            for (int j=0;j<nat;j++) vat[i] += SQ(c,j,i)*qpat[c*nat+j];
        }
    }
}

#undef SD
#undef DD
#undef SQ
