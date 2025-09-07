#include "../c_src/xtb/h0.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double rand_double(void) {
    return (double)rand() / RAND_MAX * 2.0 - 1.0;
}

int main(void) {
    srand(0);
    int nat = 3;
    int nid = 2;
    int mshell = 3;
    int nshell[2] = {2,3}; /* element 1 has 2 shells, element 2 has 3 */
    int id[3] = {1,2,1};
    int ish_at[3] = {0,2,5};
    int nsh_tot = 7;

    tb_hamiltonian h0;
    new_hamiltonian(&h0, mshell, nid);
    for (int iz=0; iz<nid; ++iz) {
        for (int ish=0; ish<mshell; ++ish) {
            int idx = ish + mshell*iz;
            h0.selfenergy[idx] = rand_double();
            h0.kcn[idx] = rand_double();
            h0.kq1[idx] = rand_double();
            h0.kq2[idx] = rand_double();
            h0.refocc[idx] = rand_double();
        }
    }

    double cn[3], qat[3];
    for (int i=0;i<nat;i++) { cn[i]=rand_double(); qat[i]=rand_double(); }

    double *se = calloc(nsh_tot, sizeof(double));
    double *dsedcn = calloc(nsh_tot, sizeof(double));
    double *dsedq = calloc(nsh_tot, sizeof(double));
    double *ref_se = calloc(nsh_tot, sizeof(double));
    double *ref_dcn = calloc(nsh_tot, sizeof(double));
    double *ref_dq = calloc(nsh_tot, sizeof(double));

    /* Reference selfenergy */
    for (int i=0;i<nat;i++) {
        int izp=id[i]; int ii=ish_at[i];
        for (int ish=0; ish<nshell[izp-1]; ++ish) {
            int idx = ish + mshell*(izp-1);
            ref_se[ii+ish] = h0.selfenergy[idx]
                - h0.kcn[idx]*cn[i]
                - h0.kq1[idx]*qat[i] - h0.kq2[idx]*qat[i]*qat[i];
            ref_dcn[ii+ish] = -h0.kcn[idx];
            ref_dq[ii+ish] = -h0.kq1[idx] - 2*h0.kq2[idx]*qat[i];
        }
    }

    get_selfenergy(&h0, nat, id, ish_at, nshell, cn, qat, se, dsedcn, dsedq);
    for (int i=0;i<nsh_tot;i++) {
        assert(fabs(se[i] - ref_se[i]) < 1e-12);
        assert(fabs(dsedcn[i] - ref_dcn[i]) < 1e-12);
        assert(fabs(dsedq[i] - ref_dq[i]) < 1e-12);
    }

    /* Occupation */
    double charge = 0.5;
    double nocc, n0at[3], n0sh[7];
    double ref_nocc = -charge; double ref_n0at[3] = {0}; double ref_n0sh[7] = {0};
    for (int i=0;i<nat;i++) {
        int iz=id[i]; int ii=ish_at[i];
        for (int ish=0; ish<nshell[iz-1]; ++ish) {
            int idx = ish + mshell*(iz-1);
            double ref = h0.refocc[idx];
            ref_nocc += ref;
            ref_n0at[i] += ref;
            ref_n0sh[ii+ish] += ref;
        }
    }
    get_occupation(&h0, nat, nsh_tot, id, ish_at, nshell, charge,
                   &nocc, n0at, n0sh);
    assert(fabs(nocc - ref_nocc) < 1e-12);
    for (int i=0;i<nat;i++) assert(fabs(n0at[i] - ref_n0at[i]) < 1e-12);
    for (int i=0;i<nsh_tot;i++) assert(fabs(n0sh[i] - ref_n0sh[i]) < 1e-12);

    /* shift_operator test */
    double vec[3], di[3], qi[6];
    for (int k=0;k<3;k++) { vec[k]=rand_double(); di[k]=rand_double(); }
    for (int k=0;k<6;k++) qi[k]=rand_double();
    double s = rand_double();
    double dj[3], qj[6];
    shift_operator(vec, s, di, qi, dj, qj);

    double ref_dj[3], ref_qj[6];
    for (int k=0;k<3;k++) ref_dj[k]=di[k]+vec[k]*s;
    ref_qj[0]=2*vec[0]*di[0]+vec[0]*vec[0]*s;
    ref_qj[2]=2*vec[1]*di[1]+vec[1]*vec[1]*s;
    ref_qj[5]=2*vec[2]*di[2]+vec[2]*vec[2]*s;
    ref_qj[1]=vec[0]*di[1]+vec[1]*di[0]+vec[0]*vec[1]*s;
    ref_qj[3]=vec[0]*di[2]+vec[2]*di[0]+vec[0]*vec[2]*s;
    ref_qj[4]=vec[1]*di[2]+vec[2]*di[1]+vec[1]*vec[2]*s;
    double tr=0.5*(ref_qj[0]+ref_qj[2]+ref_qj[5]);
    ref_qj[0]=qi[0]+1.5*ref_qj[0]-tr;
    ref_qj[1]=qi[1]+1.5*ref_qj[1];
    ref_qj[2]=qi[2]+1.5*ref_qj[2]-tr;
    ref_qj[3]=qi[3]+1.5*ref_qj[3];
    ref_qj[4]=qi[4]+1.5*ref_qj[4];
    ref_qj[5]=qi[5]+1.5*ref_qj[5]-tr;
    for (int k=0;k<3;k++) assert(fabs(dj[k]-ref_dj[k])<1e-12);
    for (int k=0;k<6;k++) assert(fabs(qj[k]-ref_qj[k])<1e-12);

    free(se); free(dsedcn); free(dsedq); free(ref_se); free(ref_dcn); free(ref_dq);
    free_hamiltonian(&h0);
    printf("h0 tests passed\n");
    return 0;
}
