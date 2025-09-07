#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xtb/gfn1.h"

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
        for (int i = 0; i < nat; ++i) {
            for (int ish = 0; ish < max_nsh; ++ish) {
                ang[i*max_nsh + ish] = (ish < nsh_id[i]) ? rand() % 3 : 0;
            }
        }
        gfn1_shell_hardness(nat, nums, nsh_id, ang, max_nsh, hard);
        for (int i = 0; i < nat; ++i) {
            for (int ish = 0; ish < nsh_id[i]; ++ish) {
                int il = ang[i*max_nsh + ish];
                double expect = gfn1_hubbard_parameter[nums[i]-1] *
                                gfn1_shell_hubbard_value(il, nums[i]);
                if (fabs(expect - hard[i*max_nsh + ish]) > 1e-12) {
                    fprintf(stderr, "mismatch at atom %d shell %d\n", i, ish);
                    return 1;
                }
            }
            for (int ish = nsh_id[i]; ish < max_nsh; ++ish) {
                if (fabs(hard[i*max_nsh + ish]) > 1e-12) {
                    fprintf(stderr, "tail not zero at atom %d shell %d\n", i, ish);
                    return 1;
                }
            }
        }
        free(ang);
        free(hard);
    }
    printf("OK\n");
    return 0;
}
