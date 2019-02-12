
////////////////////////////////////////////////////////////////////////
//
//  Field calculation
//
////////////////////////////////////////////////////////////////////////
#define Ncoilmax 100
#define nmaxmag 4000 //500
#define Ncenmax 5000
#define Ncenaxisymmmax 5000



void magfieldtraj(double *x,double *B, string filedir);
void elfieldtraj(double *x,double *E,double *Phi);
void elfield(double z,double r,double *Phi,double *Ez,double *Er);

double common[100],facB,facE;
double earth_Bx, earth_By, earth_Bz;

#include "magfield3.cpp"


////////////////////////////////////////////////////////

void magfieldtraj(double *x,double *B, string filedir){
// Computes in a point x[1],x[2],x[3] (Descartes coordinates)
// the magnetic field B[1],B[2],B[3]
  double z,r,Bz,Br,A,sine,cosine,Bm[4]; // not used
  int k,j;         // not used
//
  magfield(x,B,filedir);        // calculate B-field with magfield
  B[1]*=facB; B[2]*=facB; B[3]*=facB;   // rescale B-field
  // add earth magnetic field
  B[1] += earth_Bx*10e-6;
  B[2] += earth_By*10e-6;
  B[3] += earth_Bz*10e-6;

  return;
}

////////////////////////////////////////////////////////

void elfieldtraj(double *x,double *E,double *Phi)
{
// Computes in a point x[1],x[2],x[3] the electric field
//  E[1],E[2],E[3] and the electric potential Phi.
  double z,r,Ez,Er,sine,cosine;
  int k;
//
  z=x[3];
  r=sqrt(x[1]*x[1]+x[2]*x[2]);
//
  if(facE<1.e-10)
  {
    E[1]=E[2]=E[3]=*Phi=0.;
    return;
  }
  elfield(z,r,Phi,&Ez,&Er);
  E[3]=Ez;
  if(r<1.e-12)
    E[1]=E[2]=0.;
  else
  {
    cosine=x[1]/r; sine=x[2]/r;
    E[1]=cosine*Er;
    E[2]=sine*Er;
  }
  E[1]*=facE; E[2]*=facE; E[3]*=facE;  (*Phi)*=facE;
//
  return;
}

///////////////////////////////////////////////////////////////////

void elfield(double z,double r,double *Phi,double *Ez,double *Er)
// This subroutine computes the electric potential and the axial and radial
//   electric field components Ez and Er, in a point with cylindrical
//   coordinates z and r.
// Method: Legendre polynomial expansion around the source point with index k.
//  Nsp: number of z0 source points;
//  nmax: maximal index of the source constants (maximum of n).
//  nmax is given by global #define command.
//  Important (if elfield is used separated from elsource):
//    the same nmax number used for elsource
//    should also be used for elfield !!!
{
  FILE *fp;
  int kloop,kx,n,nmaxtest,k,kk;
  static double z0[Nspmax+1],rocen[Nspmax+1],Phin[Nspmax+1][nmax+1];
  static double c1[nmax+1],c2[nmax+1],c3[nmax+1],c4[nmax+1];
  static int Nsp,klast;
  static int iff=0;
  double ro,u,delz,s,rcmin;
  double P[nmax+1],Pp[nmax+1],rc,rcn;
  int iPhi,iEz,iEr;
  double Phiplus,Ezplus,Erplus;
  double Phiplus1,Ezplus1,Erplus1;
  double Phi1,Ez1,Er1;
//
// Input from file elsource.dat:
  if(iff==0)
  {
    fp=fopen("elsource.dat","r");
    fscanf(fp,"%i %i",&Nsp,&nmaxtest);
    if(Nsp > Nspmax || nmaxtest != nmax)
    {
      printf("Message from subroutine elfield:\
             Nsp > Nspmax or different nmax values used in elsource\
	     and elfield !!!\
             Computation is  stopped !!! \n");
      exit(0);
    }
    for(kloop=1;kloop<=Nsp;kloop++)
    {
      fscanf(fp,"%i",&kx);
      fscanf(fp,"%le %le",&z0[kloop],&rocen[kloop]);
      for(n=0;n<=nmax;n++)
        fscanf(fp,"%le",&Phin[kloop][n]);
    }
    fclose(fp);
// Initialization of c1,c2,c3,c4 vectors:
    for(n=2;n<=nmax;n++)
    {
      c1[n]=(2.*n-1.)/(1.*n);
      c2[n]=(n-1.)/(1.*n);
      c3[n]=(2.*n-1.)/(1.*(n-1.));
      c4[n]=(1.*n)/(1.*(n-1.));
    }
// The best source point is searched here (with minimal
//    convergence ratio rc):
    rcmin=1.e20;
    for(k=1;k<=Nsp;k++)
    {
      delz=z-z0[k]; ro=sqrt(r*r+delz*delz); rc=ro/rocen[k];
      if(rc<rcmin)
      { rcmin=rc; klast=k; }
    }
// End of source point searching
    iff=1;
  }
//
// The best source point is searched here
//   (starting from the last source point)
  k=klast;
  delz=z-z0[k]; ro=sqrt(r*r+delz*delz); rcmin=ro/rocen[k];
  kk=k+1;
  if(kk<=Nsp)
  { delz=z-z0[kk]; ro=sqrt(r*r+delz*delz); rc=ro/rocen[kk];
    if(rc<rcmin)
      { rcmin=rc; k=kk; }
  }
  kk=klast-1;
  if(kk>=1)
  { delz=z-z0[kk]; ro=sqrt(r*r+delz*delz); rc=ro/rocen[kk];
    if(rc<rcmin)
      k=kk;
  }
  klast=k;
// If rc>0.99: new searching:
  delz=z-z0[k]; ro=sqrt(r*r+delz*delz); rc=ro/rocen[k];
  if(rc>0.99)
  {
    rcmin=1.e20;
    for(k=1;k<=Nsp;k++)
    {
      delz=z-z0[k]; ro=sqrt(r*r+delz*delz); rc=ro/rocen[k];
      if(rc<rcmin)
      { rcmin=rc; klast=k; }
    }
    k=klast;
  }
// End of source point searching
//////////////////////////////////////
// If the field point is very close to the source point:
  if(r<1.e-12 && fabs(z-z0[k])<1.e-12)
  {
    *Phi=Phin[k][0];
    *Ez=-Phin[k][1]/rocen[k];
    *Er=0.;
    return;
  }
// ro,u,s:
  delz=z-z0[k];
  ro=sqrt(r*r+delz*delz);
  u=delz/ro;
  s=r/ro;
// First 2 terms of Legendre polynomial P and its derivative Pp (P-primed)
  P[0]=1.; P[1]=u;
  Pp[0]=0.; Pp[1]=1.;
// Convergence ratio:
  rc=ro/rocen[k];
  common[0]=rc;
// If rc>0.99: computation is stopped
//   (the series is not convergent)
  if(rc>0.99)
  {
   *Phi=1.e31;
   *Ez=0.;
   *Er=0.;
    return;
  }
// First 2 terms of the series:
  rcn=rc;
  *Phi=Phin[k][0]+Phin[k][1]*rc*u;
  *Ez=Phin[k][1]+Phin[k][2]*2.*rc*u;
  *Er=Phin[k][2]*rc;
//
  iPhi=0; iEz=0; iEr=0;
  Phiplus1=1.e30; Ezplus1=1.e30; Erplus1=1.e30;
// We start here the series expansion:
  for(n=2;n<=nmax-1;n++)
  {
    rcn*=rc;
    P[n]=c1[n]*u*P[n-1]-c2[n]*P[n-2];
    Pp[n]=c3[n]*u*Pp[n-1]-c4[n]*Pp[n-2];
    Phiplus=Phin[k][n]*rcn*P[n];
    Ezplus=Phin[k][n+1]*(n+1)*rcn*P[n];
    Erplus=Phin[k][n+1]*rcn*Pp[n];
    *Phi+=Phiplus; *Ez+=Ezplus; *Er+=Erplus;
    Phi1=1.e-15*fabs(*Phi); Ez1=1.e-15*fabs(*Ez); Er1=1.e-15*fabs(*Er);
    if(fabs(Phiplus)<Phi1 && fabs(Phiplus1)<Phi1) iPhi=1;
    if(fabs(Ezplus)<Ez1 && fabs(Ezplus1)<Ez1) iEz=1;
    if(fabs(Erplus)<Er1 && fabs(Erplus1)<Er1) iEr=1;
    if(r<1.e-12) iEr=1;
    if(iPhi*iEz*iEr == 1) break;
    Phiplus1=Phiplus; Ezplus1=Ezplus; Erplus1=Erplus;
  }
  *Ez*=-1./rocen[k];
  *Er*=s/rocen[k];
//
  common[1]=n;
  return;
}

/////////////////////////////////////////////////////////////////


