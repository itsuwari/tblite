#include "../c_src/coulomb/multipole.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static double randu(void){return (double)rand()/RAND_MAX;}

int main(void){
    srand(0);
    int nat=3; int nspecies=2;
    int id[3]; for(int i=0;i<nat;i++) id[i]=rand()%nspecies+1;
    double xyz[9]; for(int i=0;i<9;i++) xyz[i]=randu();
    double cn[3]; for(int i=0;i<3;i++) cn[i]=randu();
    double rad[2]; double vcn[2]; for(int i=0;i<2;i++){rad[i]=randu()+1.0; vcn[i]=randu();}
    double shift=0.1, kexp=1.2, rmax=3.0;
    double mrad[3], dmrdcn[3];
    multipole_mrad(nat,id,cn,rad,vcn,shift,kexp,rmax,mrad,dmrdcn);
    for(int i=0;i<nat;i++){
        int iz=id[i]-1; double arg=cn[i]-vcn[iz]-shift; double t1=exp(-kexp*arg);
        double t2=(rmax-rad[iz])/(1.0+t1);
        double ref_m=mrad[i]; double ref_d=dmrdcn[i];
        assert(fabs(ref_m-(rad[iz]+t2))<1e-12);
        assert(fabs(ref_d-(-t2*kexp*t1/(1.0+t1)))<1e-12);
    }
    double kdmp3=2.0,kdmp5=4.0;
    double *sd=calloc(3*nat*nat,sizeof(double));
    double *dd=calloc(3*nat*3*nat,sizeof(double));
    double *sq=calloc(6*nat*nat,sizeof(double));
    multipole_matrix_0d(nat,xyz,mrad,kdmp3,kdmp5,sd,dd,sq);
    for(int i=0;i<nat;i++){
      for(int j=0;j<nat;j++) if(i!=j){
        double vec[3]={xyz[3*i]-xyz[3*j], xyz[3*i+1]-xyz[3*j+1], xyz[3*i+2]-xyz[3*j+2]};
        double r1=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
        double g1=1.0/r1; double g3=g1*g1*g1; double g5=g3*g1*g1;
        double rr=0.5*(mrad[j]+mrad[i])*g1;
        double fd3=1.0/(1.0+6.0*pow(rr,kdmp3));
        double fd5=1.0/(1.0+6.0*pow(rr,kdmp5));
        for(int k=0;k<3;k++){
          double ref = vec[k]*g3*fd3;
          assert(fabs(sd[k*nat*nat + j*nat + i]-ref)<1e-12);
          for(int m=0;m<3;m++){
            double unity=(k==m)?1.0:0.0;
            double refdd = unity*g3*fd5 - vec[k]*vec[m]*3.0*g5*fd5;
            assert(fabs(dd[((k)*nat + j)*3*nat + m*nat + i]-refdd)<1e-12);
          }
        }
        double tc[6]={vec[0]*vec[0]*g5*fd5,2*vec[0]*vec[1]*g5*fd5,vec[1]*vec[1]*g5*fd5,2*vec[0]*vec[2]*g5*fd5,2*vec[1]*vec[2]*g5*fd5,vec[2]*vec[2]*g5*fd5};
        for(int c=0;c<6;c++) assert(fabs(sq[c*nat*nat + j*nat + i]-tc[c])<1e-12);
      }
    }
    double qat[3]; for(int i=0;i<3;i++) qat[i]=randu()-0.5;
    double dpat[9]; for(int i=0;i<9;i++) dpat[i]=randu()-0.5;
    double qpat[18]; for(int i=0;i<18;i++) qpat[i]=randu()-0.5;
    double energy=0.0; multipole_energy(nat,sd,dd,sq,qat,dpat,qpat,&energy);
    double ref_e=0.0; double vd[9]={0}; double vq[18]={0};
    for(int i=0;i<nat;i++){
      for(int k=0;k<3;k++){
        for(int j=0;j<nat;j++) vd[k*nat+i]+=sd[k*nat*nat + i*nat + j]*qat[j];
        for(int m=0;m<3;m++) for(int j=0;j<nat;j++) vd[k*nat+i]+=0.5*dd[((k)*nat + i)*3*nat + m*nat + j]*dpat[m*nat+j];
      }
      for(int c=0;c<6;c++) for(int j=0;j<nat;j++) vq[c*nat+i]+=sq[c*nat*nat + i*nat + j]*qat[j];
      for(int k=0;k<3;k++) ref_e+=dpat[k*nat+i]*vd[k*nat+i];
      for(int c=0;c<6;c++) ref_e+=qpat[c*nat+i]*vq[c*nat+i];
    }
    assert(fabs(energy-ref_e)<1e-12);
    double vdp[9], vat[3], vqp[18];
    multipole_potential(nat,sd,dd,sq,qat,dpat,qpat,vdp,vat,vqp);
    double rvdp[9]={0}, rvat[3]={0}, rvqp[18]={0};
    for(int i=0;i<nat;i++){
      for(int k=0;k<3;k++){
        for(int j=0;j<nat;j++) rvdp[k*nat+i]+=sd[k*nat*nat + i*nat + j]*qat[j];
        for(int m=0;m<3;m++) for(int j=0;j<nat;j++) rvdp[k*nat+i]+=dd[((k)*nat + i)*3*nat + m*nat + j]*dpat[m*nat+j];
      }
      for(int c=0;c<6;c++) for(int j=0;j<nat;j++) rvqp[c*nat+i]+=sq[c*nat*nat + i*nat + j]*qat[j];
    }
    for(int i=0;i<nat;i++){
      for(int k=0;k<3;k++) for(int j=0;j<nat;j++) rvat[i]+=sd[k*nat*nat + j*nat + i]*dpat[k*nat+j];
      for(int c=0;c<6;c++) for(int j=0;j<nat;j++) rvat[i]+=sq[c*nat*nat + j*nat + i]*qpat[c*nat+j];
    }
    for(int i=0;i<9;i++) assert(fabs(vdp[i]-rvdp[i])<1e-12);
    for(int i=0;i<3;i++) assert(fabs(vat[i]-rvat[i])<1e-12);
    for(int i=0;i<18;i++) assert(fabs(vqp[i]-rvqp[i])<1e-12);
    free(sd); free(dd); free(sq);
    printf("multipole test passed\n");
    return 0;
}
