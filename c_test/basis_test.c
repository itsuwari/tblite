#include "../c_src/basis/type.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static double randu(void){return (double)rand()/RAND_MAX;}

static double clip(double x,double mn,double mx){return x<mn?mn:(x>mx?mx:x);} 
static double integral_cutoff_ref(double acc){
    const double min_intcut=5.0,max_intcut=25.0,max_acc=1e-4,min_acc=1e3;
    double val = max_intcut - 10.0*log10(clip(acc,min_acc,max_acc));
    if (val<min_intcut) val=min_intcut; if (val>max_intcut) val=max_intcut; return val;
}

int main(void){
    srand(0);
    int nspecies=3; int nat=4;
    int nshell[3];
    for(int i=0;i<nspecies;i++) nshell[i]=rand()%3+1;
    cgto_type **cgto = malloc(nspecies*sizeof(cgto_type*));
    for(int isp=0; isp<nspecies; ++isp){
        cgto[isp]=malloc(nshell[isp]*sizeof(cgto_type));
        for(int ish=0; ish<nshell[isp]; ++ish){
            cgto[isp][ish].ang = rand()%3;
            cgto[isp][ish].nprim = rand()%3 + 1;
            for(int ip=0; ip<cgto[isp][ish].nprim; ++ip){
                cgto[isp][ish].alpha[ip] = randu()*2+0.1;
                cgto[isp][ish].coeff[ip] = randu();
            }
        }
    }
    int id[4];
    for(int i=0;i<nat;i++) id[i]=rand()%nspecies + 1;
    basis_type bas; new_basis(&bas, nat, nspecies, id, nshell, cgto, 1e-1);

    int *nsh_at_ref = malloc(sizeof(int)*nat);
    for(int i=0;i<nat;i++) nsh_at_ref[i]=nshell[id[i]-1];
    int nsh=0; for(int i=0;i<nat;i++) nsh+=nsh_at_ref[i];
    int *ish_at_ref=malloc(sizeof(int)*nat);
    int *sh2at_ref=malloc(sizeof(int)*nsh);
    int ii=0;
    for(int i=0;i<nat;i++){
        ish_at_ref[i]=ii;
        for(int ish=0; ish<nsh_at_ref[i]; ++ish) sh2at_ref[ii+ish]=i;
        ii+=nsh_at_ref[i];
    }
    int *nao_sh_ref=malloc(sizeof(int)*nsh);
    for(int i=0;i<nat;i++){
        int isp=id[i]-1; ii=ish_at_ref[i];
        for(int ish=0; ish<nsh_at_ref[i]; ++ish)
            nao_sh_ref[ii+ish]=2*cgto[isp][ish].ang+1;
    }
    int nao=0; for(int ish=0; ish<nsh; ++ish) nao+=nao_sh_ref[ish];
    int *iao_sh_ref=malloc(sizeof(int)*nsh);
    int *ao2sh_ref=malloc(sizeof(int)*nao);
    int *ao2at_ref=malloc(sizeof(int)*nao);
    ii=0;
    for(int ish=0; ish<nsh; ++ish){
        iao_sh_ref[ish]=ii;
        for(int iao=0; iao<nao_sh_ref[ish]; ++iao){
            ao2sh_ref[ii+iao]=ish; ao2at_ref[ii+iao]=sh2at_ref[ish];
        }
        ii+=nao_sh_ref[ish];
    }
    int maxl_ref=0; double min_alpha_ref=1e308;
    for(int isp=0; isp<nspecies; ++isp){
        for(int ish=0; ish<nshell[isp]; ++ish){
            if(cgto[isp][ish].ang>maxl_ref) maxl_ref=cgto[isp][ish].ang;
            for(int ip=0; ip<cgto[isp][ish].nprim; ++ip)
                if(cgto[isp][ish].alpha[ip]<min_alpha_ref) min_alpha_ref=cgto[isp][ish].alpha[ip];
        }
    }

    assert(bas.nsh == nsh);
    assert(bas.nao == nao);
    assert(bas.maxl == maxl_ref);
    assert(fabs(bas.min_alpha - min_alpha_ref) < 1e-12);
    for(int i=0;i<nat;i++) assert(bas.nsh_at[i]==nsh_at_ref[i]);
    for(int i=0;i<nat;i++) assert(bas.ish_at[i]==ish_at_ref[i]);
    for(int ish=0; ish<nsh; ++ish){
        assert(bas.sh2at[ish]==sh2at_ref[ish]);
        assert(bas.nao_sh[ish]==nao_sh_ref[ish]);
        assert(bas.iao_sh[ish]==iao_sh_ref[ish]);
    }
    for(int iao=0; iao<nao; ++iao){
        assert(bas.ao2sh[iao]==ao2sh_ref[iao]);
        assert(bas.ao2at[iao]==ao2at_ref[iao]);
    }

    double cutoff = basis_get_cutoff(&bas, 0.0);
    double ref_cutoff = sqrt(2.0*bas.intcut/bas.min_alpha);
    if(ref_cutoff>40.0) ref_cutoff=40.0;
    assert(fabs(cutoff - ref_cutoff) < 1e-12);

    double cutoff2 = basis_get_cutoff(&bas, 1e-2);
    double ref_int = integral_cutoff_ref(1e-2);
    double ref_cutoff2 = sqrt(2.0*ref_int/bas.min_alpha);
    if(ref_cutoff2>40.0) ref_cutoff2=40.0;
    assert(fabs(cutoff2 - ref_cutoff2) < 1e-12);

    free_basis(&bas);
    for(int isp=0; isp<nspecies; ++isp) free(cgto[isp]);
    free(cgto);
    free(nsh_at_ref); free(ish_at_ref); free(sh2at_ref);
    free(nao_sh_ref); free(iao_sh_ref); free(ao2sh_ref); free(ao2at_ref);
    printf("basis test passed\n");
    return 0;
}
