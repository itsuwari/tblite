#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xtb/gfn2.h"
#include "xtb/h0.h"

int main(void) {
    srand(0);
    for (int t = 0; t < 10; ++t) {
        int nat = 1 + rand() % 5;
        int nums[5];
        int nsh_id[5];
        int max_nsh = 0;
        for (int i = 0; i < nat; ++i) {
            nums[i] = 1 + rand() % 86;
            nsh_id[i] = 1 + rand() % 3;
            if (nsh_id[i] > max_nsh) max_nsh = nsh_id[i];
        }
        int *ang = malloc(sizeof(int) * nat * max_nsh);
        double *hard = malloc(sizeof(double) * nat * max_nsh);
        double *deriv = malloc(sizeof(double) * nat * max_nsh);
        for (int i = 0; i < nat; ++i) {
            for (int ish = 0; ish < max_nsh; ++ish) {
                ang[i*max_nsh + ish] = (ish < nsh_id[i]) ? rand() % 3 : 0;
            }
        }
        gfn2_shell_hardness(nat, nums, nsh_id, ang, max_nsh, hard);
        gfn2_hubbard_derivs(nat, nums, nsh_id, ang, max_nsh, deriv);
        for (int i = 0; i < nat; ++i) {
            for (int ish = 0; ish < nsh_id[i]; ++ish) {
                int il = ang[i*max_nsh + ish];
                double expect = gfn2_hubbard_parameter[nums[i]-1] *
                                gfn2_shell_hubbard_value(il, nums[i]);
                if (fabs(expect - hard[i*max_nsh + ish]) > 1e-12) {
                    fprintf(stderr, "hard mismatch atom %d shell %d\n", i, ish);
                    return 1;
                }
                double expd = gfn2_p_hubbard_derivs[nums[i]-1] *
                               gfn2_shell_hubbard_derivs[il];
                if (fabs(expd - deriv[i*max_nsh + ish]) > 1e-12) {
                    fprintf(stderr, "deriv mismatch atom %d shell %d\n", i, ish);
                    return 1;
                }
            }
            for (int ish = nsh_id[i]; ish < max_nsh; ++ish) {
                if (fabs(hard[i*max_nsh + ish]) > 1e-12 ||
                    fabs(deriv[i*max_nsh + ish]) > 1e-12) {
                    fprintf(stderr, "tail not zero atom %d shell %d\n", i, ish);
                    return 1;
                }
            }
        }
        free(ang);
        free(hard);
        free(deriv);
    }
    /* Hamiltonian parameter wiring test */
    int zlist[2] = {1, 8};
    tb_hamiltonian h0;
    new_hamiltonian(&h0, 3, 2);
    gfn2_init_h0(&h0, zlist);
    int nshell[2] = {gfn2_nshell[zlist[0]-1], gfn2_nshell[zlist[1]-1]};
    int nat = 2;
    int id[2] = {1,2};
    int ish_at[2] = {0, nshell[0]};
    double cn[2] = {0.3, -0.1};
    double q[2] = {0.05, -0.02};
    int nsh_tot = nshell[0] + nshell[1];
    double *se = calloc(nsh_tot, sizeof(double));
    get_selfenergy(&h0, nat, id, ish_at, nshell, cn, q, se, NULL, NULL);
    for (int i=0;i<nat;i++) {
        int z = zlist[i];
        for (int l=0;l<nshell[i];++l) {
            double expect = gfn2_selfenergy_value(l,z) - gfn2_kcn_value(l,z)*cn[i];
            if (fabs(se[ish_at[i]+l] - expect) > 1e-12) {
                fprintf(stderr, "selfenergy mismatch\n");
                return 1;
            }
        }
    }
    free(se);
    double nocc, n0at[2], n0sh[6];
    get_occupation(&h0, nat, nsh_tot, id, ish_at, nshell, 0.0, &nocc, n0at, n0sh);
    double ref_nocc = 0.0; double ref_n0at[2] = {0,0};
    for (int i=0;i<nat;i++) {
        int z = zlist[i];
        for (int l=0;l<nshell[i];++l) {
            double val = gfn2_refocc_value(l,z);
            ref_nocc += val;
            ref_n0at[i] += val;
            if (fabs(n0sh[ish_at[i]+l] - val) > 1e-12) {
                fprintf(stderr, "refocc shell mismatch\n");
                return 1;
            }
        }
        if (fabs(n0at[i] - ref_n0at[i]) > 1e-12) {
            fprintf(stderr, "refocc atom mismatch\n");
            return 1;
        }
    }
    if (fabs(nocc - ref_nocc) > 1e-12) {
        fprintf(stderr, "refocc total mismatch\n");
        return 1;
    }
    free_hamiltonian(&h0);
    printf("OK\n");
    return 0;
}
