#include "../c_src/scf/mixer.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>

int main(){
    mixer_type *mix;
    new_mixer(&mix, 5, 10, 0.5);
    double q0[10], q1[10], out[10];
    for(int i=0;i<10;i++){ q0[i]=((double)rand()/RAND_MAX); q1[i]=((double)rand()/RAND_MAX); }
    mix->set(mix, q0,10);
    mix->diff(mix, q1,10);
    mix->next(mix);
    mix->get(mix, out,10);
    for(int i=0;i<10;i++){
        double expected = q0[i] + 0.5*(q1[i]-q0[i]);
        assert(fabs(out[i]-expected) < 1e-12);
    }
    double err = mix->get_error(mix);
    double ref=0.0; for(int i=0;i<10;i++){ double d=q1[i]-q0[i]; ref += d*d/10.0; }
    ref = sqrt(ref);
    assert(fabs(err-ref) < 1e-12);
    mix->del(mix);
    free(mix);
    return 0;
}
