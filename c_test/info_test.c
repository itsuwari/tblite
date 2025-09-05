#include "../c_src/scf/info.h"
#include <assert.h>
#include <stdlib.h>
#include <time.h>

int main(){
    srand(0);
    for(int i=0;i<100;i++){
        scf_info a={rand()%4, rand()%4, rand()%4};
        scf_info b={rand()%4, rand()%4, rand()%4};
        scf_info c=max_info(a,b);
        assert(c.charge== (a.charge>b.charge?a.charge:b.charge));
        assert(c.dipole== (a.dipole>b.dipole?a.dipole:b.dipole));
        assert(c.quadrupole== (a.quadrupole>b.quadrupole?a.quadrupole:b.quadrupole));
    }
    return 0;
}
