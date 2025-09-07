#include "../c_src/scf/iterator.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>

int main(){
    int norb=4, nspin=2;
    double *h0 = malloc(sizeof(double)*norb*norb);
    double *density = malloc(sizeof(double)*norb*norb*nspin);
    double *energies = calloc(norb, sizeof(double));
    for(int i=0;i<norb*norb;i++) h0[i]=((double)rand()/RAND_MAX);
    for(int i=0;i<norb*norb*nspin;i++) density[i]=((double)rand()/RAND_MAX);
    get_electronic_energy(h0,density,energies,norb,nspin);
    for(int iao=0; iao<norb; ++iao){
        double ref=0.0;
        for(int spin=0; spin<nspin; ++spin){
            for(int jao=0; jao<norb; ++jao){
                ref += h0[jao + iao*norb]*density[jao + iao*norb + spin*norb*norb];
            }
        }
        assert(fabs(ref-energies[iao])<1e-12);
    }
    free(h0); free(density); free(energies);

    double red[3]={0,0,0};
    double full[4]={1,2,3,4};
    int map[4]={0,1,1,2};
    reduce(red,full,map,4);
    assert(red[0]==1 && red[1]==5 && red[2]==4);

    basis_type bas; bas.nat=3; bas.nsh=4; bas.nao=0;
    int ish_at_arr[3]={0,1,3};
    int nsh_at_arr[3]={1,2,1};
    int iao_sh_arr[4]={0,0,0,0};
    int nao_sh_arr[4]={0,0,0,0};
    bas.ish_at=ish_at_arr; bas.nsh_at=nsh_at_arr;
    bas.iao_sh=iao_sh_arr; bas.nao_sh=nao_sh_arr; bas.ao2at=NULL;
    double qsh[8]; for(int i=0;i<8;i++) qsh[i]=i+1;
    double qat[6];
    get_qat_from_qsh(&bas,qsh,qat,4,2);
    assert(qat[0]==1 && qat[1]==2+3 && qat[2]==4);
    assert(qat[3]==5 && qat[4]==6+7 && qat[5]==8);

    structure_type mol; mol.nat=3; bas.nsh=4;
    scf_info info={scf_shell_resolved, scf_atom_resolved, scf_not_used};
    int ndim=get_mixer_dimension(&mol,&bas,info);
    assert(ndim== bas.nsh + 3*mol.nat);
    return 0;
}
