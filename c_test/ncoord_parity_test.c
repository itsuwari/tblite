#include "../c_src/ncoord/gfn.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Fortran reference routines */
extern double ncoord_count_ref(double rc_i, double rc_j, double r);
extern double ncoord_dcount_ref(double rc_i, double rc_j, double r);

static double rand_pos(void){ return (double)rand()/RAND_MAX + 0.1; }

int main(){
    srand(0);
    gfn_ncoord_type cn;
    double rc[2] = {rand_pos()*2.0, rand_pos()*2.0};
    new_gfn_ncoord(&cn, 2, rc, 0.0);
    for(int t=0;t<100;t++){
        double r = rand_pos()*10.0;
        double c_val = ncoord_count(&cn,0,1,r);
        double f_val = ncoord_count_ref(rc[0],rc[1],r);
        assert(fabs(c_val-f_val) < 1e-12);
        double cd_val = ncoord_dcount(&cn,0,1,r);
        double fd_val = ncoord_dcount_ref(rc[0],rc[1],r);
        assert(fabs(cd_val-fd_val) < 1e-12);
    }
    free(cn.rcov);
    printf("ncoord parity tests passed\n");
    return 0;
}
