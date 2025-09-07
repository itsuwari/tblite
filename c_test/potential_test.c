#include "../c_src/scf/potential.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define IDX2(i,j,n1) ((j)*(n1)+(i))
#define IDX3(i,j,k,n1,n2) ((k)*(n1)*(n2) + (j)*(n1) + (i))

double rand_double(void) {
    return (double)rand() / RAND_MAX * 2.0 - 1.0;
}

int main(void) {
    srand(0);
    int nat = 2, nsh = 3, nao = 4, nspin = 2;
    structure_type mol = { nat };
    basis_type bas;
    bas.nat = nat; bas.nsh = nsh; bas.nao = nao;
    int ish_at_arr[2] = {0,2};
    int nsh_at_arr[2] = {2,1};
    int iao_sh_arr[3] = {0,1,3};
    int nao_sh_arr[3] = {1,2,1};
    int ao2at_arr[4] = {0,0,1,1};
    bas.ish_at = ish_at_arr;
    bas.nsh_at = nsh_at_arr;
    bas.iao_sh = iao_sh_arr;
    bas.nao_sh = nao_sh_arr;
    bas.ao2at = ao2at_arr;

    integral_type ints;
    ints.hamiltonian = malloc(sizeof(double)*nao*nao);
    ints.overlap = malloc(sizeof(double)*nao*nao);
    ints.dipole = malloc(sizeof(double)*3*nao*nao);
    ints.quadrupole = malloc(sizeof(double)*6*nao*nao);
    for (int i=0;i<nao*nao;i++) {
        ints.hamiltonian[i] = rand_double();
        ints.overlap[i] = rand_double();
    }
    for (int i=0;i<3*nao*nao;i++) ints.dipole[i] = rand_double();
    for (int i=0;i<6*nao*nao;i++) ints.quadrupole[i] = rand_double();

    potential_type pot;
    new_potential(&pot, &mol, &bas, nspin, 0);

    for (int i=0;i<nat*nspin;i++) pot.vat[i] = rand_double();
    for (int i=0;i<nsh*nspin;i++) pot.vsh[i] = rand_double();
    for (int i=0;i<nao*nspin;i++) pot.vao[i] = rand_double();
    for (int i=0;i<3*nat*nspin;i++) pot.vdp[i] = rand_double();
    for (int i=0;i<6*nat*nspin;i++) pot.vqp[i] = rand_double();

    /* Copies for reference computation */
    double *ref_vsh = malloc(sizeof(double)*nsh*nspin);
    double *ref_vao = malloc(sizeof(double)*nao*nspin);
    memcpy(ref_vsh, pot.vsh, sizeof(double)*nsh*nspin);
    memcpy(ref_vao, pot.vao, sizeof(double)*nao*nspin);

    /* Reference expansion from vat to vsh */
    for (int spin=0; spin<nspin; ++spin) {
        for (int iat=0; iat<nat; ++iat) {
            int ii = bas.ish_at[iat];
            double val = pot.vat[IDX2(iat, spin, nat)];
            for (int ish=0; ish<bas.nsh_at[iat]; ++ish) {
                ref_vsh[IDX2(ii+ish, spin, nsh)] += val;
            }
        }
    }

    /* Reference expansion from vsh to vao */
    for (int spin=0; spin<nspin; ++spin) {
        for (int ish=0; ish<nsh; ++ish) {
            int ii = bas.iao_sh[ish];
            double val = ref_vsh[IDX2(ish, spin, nsh)];
            for (int iao=0; iao<bas.nao_sh[ish]; ++iao) {
                ref_vao[IDX2(ii+iao, spin, nao)] += val;
            }
        }
    }

    size_t hsz = (size_t)nao*nao*nspin;
    double *h1 = malloc(sizeof(double)*hsz);
    double *ref_h1 = malloc(sizeof(double)*hsz);

    /* Reference h1 */
    memcpy(ref_h1, ints.hamiltonian, sizeof(double)*nao*nao);
    memset(ref_h1 + nao*nao, 0, sizeof(double)*nao*nao*(nspin-1));

    for (int spin=0; spin<nspin; ++spin) {
        for (int iao=0; iao<nao; ++iao) {
            for (int jao=0; jao<nao; ++jao) {
                double s = ints.overlap[IDX2(jao, iao, nao)];
                double v = ref_vao[IDX2(jao, spin, nao)] + ref_vao[IDX2(iao, spin, nao)];
                ref_h1[IDX3(jao, iao, spin, nao, nao)] -= 0.5 * s * v;
            }
        }
    }

    /* Dipole */
    for (int spin=0; spin<nspin; ++spin) {
        for (int iao=0; iao<nao; ++iao) {
            for (int jao=0; jao<nao; ++jao) {
                int iat = bas.ao2at[iao];
                int jat = bas.ao2at[jao];
                double dot1=0.0, dot2=0.0;
                for (int k=0;k<3;k++) {
                    double m1 = ints.dipole[IDX3(k, jao, iao, 3, nao)];
                    double m2 = ints.dipole[IDX3(k, iao, jao, 3, nao)];
                    double v1 = pot.vdp[IDX3(k, iat, spin, 3, nat)];
                    double v2 = pot.vdp[IDX3(k, jat, spin, 3, nat)];
                    dot1 += m1*v1;
                    dot2 += m2*v2;
                }
                ref_h1[IDX3(jao, iao, spin, nao, nao)] -= 0.5*(dot1+dot2);
            }
        }
    }

    /* Quadrupole */
    for (int spin=0; spin<nspin; ++spin) {
        for (int iao=0; iao<nao; ++iao) {
            for (int jao=0; jao<nao; ++jao) {
                int iat = bas.ao2at[iao];
                int jat = bas.ao2at[jao];
                double dot1=0.0, dot2=0.0;
                for (int k=0;k<6;k++) {
                    double m1 = ints.quadrupole[IDX3(k, jao, iao, 6, nao)];
                    double m2 = ints.quadrupole[IDX3(k, iao, jao, 6, nao)];
                    double v1 = pot.vqp[IDX3(k, iat, spin, 6, nat)];
                    double v2 = pot.vqp[IDX3(k, jat, spin, 6, nat)];
                    dot1 += m1*v1;
                    dot2 += m2*v2;
                }
                ref_h1[IDX3(jao, iao, spin, nao, nao)] -= 0.5*(dot1+dot2);
            }
        }
    }

    /* Compute with implementation */
    add_pot_to_h1(&bas, &ints, &pot, h1);

    /* Compare */
    for (size_t i=0;i<nsh*nspin;i++) {
        assert(fabs(pot.vsh[i] - ref_vsh[i]) < 1e-12);
    }
    for (size_t i=0;i<nao*nspin;i++) {
        assert(fabs(pot.vao[i] - ref_vao[i]) < 1e-12);
    }
    for (size_t i=0;i<hsz;i++) {
        assert(fabs(h1[i] - ref_h1[i]) < 1e-12);
    }

    /* Test reset */
    reset_potential(&pot);
    for (size_t i=0;i<nat*nspin;i++) assert(pot.vat[i]==0.0);
    for (size_t i=0;i<nsh*nspin;i++) assert(pot.vsh[i]==0.0);
    for (size_t i=0;i<nao*nspin;i++) assert(pot.vao[i]==0.0);
    for (size_t i=0;i<3*nat*nspin;i++) assert(pot.vdp[i]==0.0);
    for (size_t i=0;i<6*nat*nspin;i++) assert(pot.vqp[i]==0.0);

    free_potential(&pot);
    free(ref_vsh); free(ref_vao); free(ref_h1); free(h1);
    free(ints.hamiltonian); free(ints.overlap); free(ints.dipole); free(ints.quadrupole);
    printf("potential tests passed\n");
    return 0;
}
