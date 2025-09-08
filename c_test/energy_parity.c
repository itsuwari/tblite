#include <stdio.h>
#include <math.h>
#include "tblite.h"

static int check_energy(tblite_context ctx, tblite_error err,
                        const char *name, int nat, const int *nums,
                        const double *xyz, double ref){
    tblite_structure mol = tblite_new_structure(err, nat, nums, xyz, NULL, NULL, NULL, NULL);
    tblite_calculator calc = tblite_new_gfn2_calculator(ctx, mol);
    tblite_result res = tblite_new_result();
    tblite_get_singlepoint(ctx, mol, calc, res);
    double energy = 0.0;
    tblite_get_result_energy(err, res, &energy);
    printf("%s energy: %.12f Eh\n", name, energy);
    int fail = fabs(energy - ref) > 1e-6;
    if (fail){
        fprintf(stderr, "Energy mismatch for %s: %.12f vs %.12f\n", name, energy, ref);
    }
    tblite_delete(res);
    tblite_delete(calc);
    tblite_delete(mol);
    return fail;
}

int main(void){
    tblite_error err = tblite_new_error();
    tblite_context ctx = tblite_new_context();
    int fail = 0;

    /* H2O */
    int nums_h2o[3] = {8,1,1};
    double xyz_h2o[9] = {
        0.0, 0.0, 0.0,
        0.0, -1.430825032522, 1.107870837823,
        0.0,  1.430825032522, 1.107870837823};
    fail |= check_energy(ctx, err, "H2O", 3, nums_h2o, xyz_h2o, -5.070371370817);

    /* NH3 */
    int nums_nh3[4] = {7,1,1,1};
    double xyz_nh3[12] = {
        0.0, 0.0, 0.0,
        0.0, 1.771996187062, -0.721119489157,
       -1.533890695359, -0.885903607225, -0.721119489157,
        1.533890695359, -0.885903607225, -0.721119489157};
    fail |= check_energy(ctx, err, "NH3", 4, nums_nh3, xyz_nh3, -4.426173749101);

    /* CH4 */
    int nums_ch4[5] = {6,1,1,1,1};
    double xyz_ch4[15] = {
        0.0, 0.0, 0.0,
        0.0, 0.0, 2.057911749717,
        1.940217716950, 0.0, -0.685970583239,
       -0.970109803338, 1.680278329603, -0.685970583239,
       -0.970109803338, -1.680278329603, -0.685970583239};
    fail |= check_energy(ctx, err, "CH4", 5, nums_ch4, xyz_ch4, -4.175098770821);

    tblite_delete(ctx);
    tblite_delete(err);
    return fail;
}
