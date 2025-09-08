#include <stdio.h>
#include <math.h>
#include "tblite.h"

int main(void){
    tblite_error err = tblite_new_error();
    tblite_context ctx = tblite_new_context();

    const int nums[3] = {8,1,1};
    const double xyz[9] = {
        0.0, 0.0, 0.0,
        0.0, -1.430825032522, 1.107870837823,
        0.0,  1.430825032522, 1.107870837823};

    tblite_structure mol = tblite_new_structure(err, 3, nums, xyz, NULL, NULL, NULL, NULL);
    tblite_calculator calc = tblite_new_gfn2_calculator(ctx, mol);
    tblite_result res = tblite_new_result();
    tblite_get_singlepoint(ctx, mol, calc, res);

    double energy = 0.0;
    tblite_get_result_energy(err, res, &energy);
    printf("H2O total energy: %.12f Eh\n", energy);
    const double ref = -5.070371370817;
    if (fabs(energy - ref) > 1e-6){
        fprintf(stderr, "Energy mismatch: expected %.12f Eh\n", ref);
    }

    tblite_delete(res);
    tblite_delete(calc);
    tblite_delete(mol);
    tblite_delete(ctx);
    tblite_delete(err);
    return 0;
}
