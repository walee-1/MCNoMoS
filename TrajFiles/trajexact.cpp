#include <stdio.h>
#include <math.h>

void trajexact(double *x,double *v);
void deriv(double *y,double *f);
void rungekutta4(int n,double h,double *y);
void rungekutta8(int n,double h,double *y);

// Mass and charge of the particle, to be defined by the user:
long double MASS,CHARGE;

struct typetrajexact
{
// Variables needed to be defined by the user:
  int trajstart,typeh,typerungekutta,ntimestep;
  double errsmall,errlarge,slimit,huser;
  string filepath;

// Output variables:
  double energy,Ekin,h,s,Phi;
  int index,indexrungekutta8;
}
commontrajexact;

////////////////////////////////////////////////////////

void trajexact(double *x,double *v)
{
// Computes new space coordinates and velocity
// components after a Runge-Kutta step.
// x[1], x[2], x[3] : input and output space coordinates
// v[1], v[2], v[3]: input and output velocity components
//  (input: at the beginning of the Runge-Kutta step,
//   output: at the end of the Runge-Kutta step).
// Input global parameters:
// commontrajexact.trajstart:  in the beginning of each new trajectory
//             calculation the user has to initialize it to zero
// commontrajexact.typeh=1:  the timestep h is defined at each call by the user
//            (h=commontrajexact.huser)
// commontrajexact.typeh=2:  the timestep h is computed automatically by the
//           program: h=cyclotron period / commontrajexact.ntimestep
// commontrajexact.typeh=3:  the timestep h is computed automatically by the
//           program: at one Runge-Kutta step
//           the energy conservation error should be between
//           commontrajexact.errsmall and commontrajexact.errlarge,
//           and the pathlength of the
//           particle should be not larger than commontrajexact.slimit.
// commontrajexact.typerungekutta=4: diff. eq. solution with rungekutta4
// commontrajexact.typerungekutta=8: diff. eq. solution with rungekutta8
// Needs mass and charge of particle:  MASS, CHARGE (SI).
// Needs magnetic field calculation:  magfieldtraj(x,B)
//     (input: space point x[1], x[2], x[3],
//      output: magnetic field B[1], B[2], B[3]).
// Needs electric field calculation:  elfieldtraj(x,E,Phi)
//     (input: space point x[1], x[2], x[3],
//      output: electric field E[1], E[2], E[3], potential Phi).
// Additional output parameters:
//     commontrajexact.energy: relativistic kinetic + potential energy (eV)
//     commontrajexact.Ekin:   relativistic kinetic energy (eV)
//     commontrajexact.h:    time value of RK-step (s)
//     commontrajexact.s:    pathlength of RK-step(m)
//     commontrajexact.Phi:    potential (V)
//  SI units are used (except at energy and kin_energy).
// In the beginning of the RK-step commontrajexact.index=0. If the potential
// becomes unphysical (Phi>1.e30) (f.e.: because the particle hits an electrode),
// the value of commontrajexact.index will be 1, and
// the trajectory evaluation stops.
  double PI = 3.1415926;
    const double c=299792458.;  // velocity of light in SI units
  const double e=1.602177e-19;  // electron charge (in SI, without sign)
  double v2,V,gamma,p[4],p2,y[7],delt,ystart[7],B[4],b,omega,y0[7];
  double energy0,energy,E[4],Phi,chper,Ekin,err,s,h1,h2;
  static double h;
  int j,iloop,k,i;
//
  commontrajexact.index=0;
    // Velocity:
    v2=v[1]*v[1]+v[2]*v[2]+v[3]*v[3];
    V=sqrt(v2);
    if(V<1.e-10) V=1.e-10;
    // Relativistic momentum components:
    gamma=1./sqrt(1.-v2/(c*c));
    for(j=1;j<=3;j++)
        p[j]=gamma*MASS*v[j];  // relativistic momentum (in SI)
    p2=p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
// Timestep initialization for commontrajexact.typeh=3:
  if(commontrajexact.trajstart==0)
  {
    if(commontrajexact.typeh==3)
    {
      h1=commontrajexact.slimit/V;
      magfieldtraj(x,B,commontrajexact.filepath);
      b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);
      if(b<1.e-10) b=1.e-10;
      omega=e*b/(MASS*gamma);
      if(commontrajexact.ntimestep<1)
      {
        printf("Message from subroutine trajexact:\
               commontrajexact.ntimestep  should be positive !!!\
               Computation is  stopped !!! \n");
        exit(0);
      }
      h2=(2.*PI/omega)/commontrajexact.ntimestep;
      if(h1<h2)
        h=h1;
      else
        h=h2;
    }
    commontrajexact.trajstart=1;
  }
//
// Timestep for commontrajexact.typeh=1 and commontrajexact.typeh=2:
  if(commontrajexact.typeh==1)
    h=commontrajexact.huser;
  else if(commontrajexact.typeh==2)
  {
    if(commontrajexact.ntimestep<1)
    {
        printf("Message from subroutine trajexact:\
               commontrajexact.ntimestep  should be positive !!!\
               Computation is  stopped !!! \n");
        exit(0);
    }
    magfieldtraj(x,B,commontrajexact.filepath);
    b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);
    omega=e*b/(MASS*gamma);
    h=(2.*PI/omega)/commontrajexact.ntimestep;
    s=V*h;
    if(s>commontrajexact.slimit)
      h=commontrajexact.slimit/V;
  }
// Starting y vector (needed for the Runge-Kutta program):
  for(j=1;j<=3;j++)
  {
    y[j]=x[j];
    y[j+3]=p[j];
  }
// we save the starting parameters y:
  for(j=1;j<=6;j++)
    ystart[j]=y[j];
// Starting total energy:
  chper=1./e;
  if(commontrajexact.typeh==3)
  {
    Ekin=p2/((gamma+1.)*MASS)*chper;  // eV
    elfieldtraj(x,E,&Phi);
    if(Phi>1.e30)
    { commontrajexact.index=1;  return; }
    energy0=Ekin+CHARGE*Phi*chper;    // eV
  }
// Pathlength-limit criterion:
  s=V*h;
  if(commontrajexact.typeh==3 && s>commontrajexact.slimit)
    h=commontrajexact.slimit/V;
  iloop=0;
label1: ;
// Runge-Kutta step:
  if(commontrajexact.typerungekutta==4)
    rungekutta4(6,h,y);
  if(commontrajexact.typerungekutta==8)
  {
    for(i=1;i<=6;i++)
      y0[i]=y[i];
    for(k=1;k<=20;k++)
    {
      rungekutta8(6,h,y);
      if(commontrajexact.indexrungekutta8==1) break;
      h*=0.5;
      for(i=1;i<=6;i++)
        y[i]=y0[i];
    }
  }
  if(commontrajexact.typerungekutta!=4 && commontrajexact.typerungekutta!=8)
  {
    printf("Message from subroutine trajexact:\
               commontrajexact.typerungekutta  is neither 4 nor 8 !!!\
               Computation is  stopped !!! \n");
    exit(0);
  }
  if(commontrajexact.index==1) return;
// New x[j], p[j], v[j] :
  for(j=1;j<=3;j++)
  {
    x[j]=y[j];
    p[j]=y[j+3];
  }
  p2=p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
  gamma=sqrt(1.+p2/(MASS*c*MASS*c));
  for(j=1;j<=3;j++)
    v[j]=p[j]/(MASS*gamma);
// New energy value:
  Ekin=commontrajexact.Ekin=p2/((gamma+1.)*MASS)*chper;
  elfieldtraj(x,E,&Phi);
  if(Phi>1.e30)
  { commontrajexact.index=1;  return; }
  energy=commontrajexact.energy=Ekin+CHARGE*Phi*chper;
  commontrajexact.Phi=Phi;
//
  if(commontrajexact.typeh==3)
  {
//  energy conservation calc.:
    err=(energy-energy0)/energy0;
// if accuracy of energy conservation is too bad,
//   the RK step is evaluated again, with smaller timestep value;
    if(commontrajexact.errlarge<1.e-13)
      commontrajexact.errlarge=1.e-13;
    if(fabs(err)>commontrajexact.errlarge)
    {
       iloop+=1;
       if(iloop==20)
         printf("Message from subroutine trajexact:\
                 danger of infinite loop !!!\
                 Larger value for 'commontrajexact.errlarge' should be used !!!   iloop= %11i \n",iloop);
       if(iloop==30)
       {
         printf("Message from subroutine trajexact:\
                 danger of infinite loop !!!\
                 Larger value for 'commontrajexact.errlarge' should be used !!!\
		 Program execution is stopped !!!  iloop= %11i \n",iloop);
         exit(0);
       }
       for(j=1;j<=6;j++)  y[j]=ystart[j];
       h=h/2.5;
       goto label1;
    }
  }
//
  commontrajexact.h=h;
  commontrajexact.s=V*h;
//
// if accuracy of energy conservation is too good,
//   the timestep value is increased (for the next RK-step)
  if(commontrajexact.typeh==3 && fabs(err)<commontrajexact.errsmall)
    h=1.5*h;
//
  return;
}


//////////////////////////////////////////////////////////

void deriv(double *y,double *f)
{
// Derivative function for the relativistic particle motion in electric and
// magnetic field.
// Input:
//   y[1], y[2], y[3]:  x, y and z coordinates of the particle;
//   y[4], y[5], y[6]:  relativistic momentum components of the particle.
//                      (px,py,pz).
// Needs mass and charge of particle:  MASS, CHARGE (global variables).
// Needs magnetic field calculation:  magfieldtraj(x,B)
//     (input: space point x[1], x[2], x[3],
//      output: magnetic field B[1], B[2], B[3]).
// Needs electric field calculation:  elfieldtraj(x,E,Phi)
//     (input: space point x[1], x[2], x[3],
//      output: electric field E[1], E[2], E[3], potential Phi).
// Output:
//   f[1],...,f[6]: derivative function components;
//  f[1],f[2],f[3]: velocity components (vx,vy,vz)
//  f[4],f[5],f[6]: Lorentz force components (Fx,Fy,Fz).
//  SI units are used.
  extern long double MASS,CHARGE;
  double fun;
  int i,j;
  double x[4],v[4],E[4],B[4],p2,s,Phi;
  static int m1[4]={0,2,3,1},m2[4]={0,3,1,2};
  const double c=299792458.;  // velocity of light in SI units
//
  p2=y[4]*y[4]+y[5]*y[5]+y[6]*y[6];
  s=1./sqrt(MASS*MASS+p2/(c*c));
  for(j=1;j<=3;j++)
   {x[j]=y[j]; v[j]=y[j+3]*s;}
  for(i=1;i<=3;i++)
    f[i]=v[i];
    // Magnetic and electric field calculation:
    magfieldtraj(x,B,commontrajexact.filepath);
    elfieldtraj(x,E,&Phi);
    if(Phi>1.e30)
  { commontrajexact.index=1;   return; }
    // Lorentz force:
    for(i=4;i<=6;i++){
        j=i-3;
        f[i]=CHARGE*(E[j]+v[m1[j]]*B[m2[j]]-v[m2[j]]*B[m1[j]]);
        }
    return;
}


///////////////////////////////////////////////////////////

void rungekutta4(int n,double h,double *y)
{
/* General 4th order Runge-Kutta program for the calculation of the
  first order differential equation system
  der_t y_i(t)=f_i(y_1,...,y_n)    /i=1,...n/
    with initial conditions
      y_i(t_0)=y_i^0                     /i=1,...n/.
  Input:
    n: number of equations in the system;
       here n cannot be larger than 100.
    y0=y: n-vector, the initial condition (at starting the subr.)
    h: the "time" step
  Output:
    y: n-vector, the y_i values at "time" t0+h (at the end)
  This subr. uses the subr. deriv, which computes
    the "derivative" function of the system: f[1],...,f[n]
  The y vector have indices i from 1 to n
  (in the calling program it should be defined as y[n+1]),
   but the i=0 values are not used in the calculation.
*/
  int i;
  double k1[101],k2[101],k3[101],k4[101],f[101],y0[101];
  for(i=1;i<=n;i++)
    y0[i]=y[i];
// k1[i] calc.:
  deriv(y,f);
  if(commontrajexact.index==1) return;
  for(i=1;i<=n;i++)
    k1[i]=f[i];
// k2[i] calc.:
  for(i=1;i<=n;i++)
    y[i]=y0[i]+h*k1[i]/2.;
  deriv(y,f);
  if(commontrajexact.index==1) return;
  for(i=1;i<=n;i++)
    k2[i]=f[i];
// k3[i] calc.:
  for(i=1;i<=n;i++)
    y[i]=y0[i]+h*k2[i]/2.;
  deriv(y,f);
  if(commontrajexact.index==1) return;
  for(i=1;i<=n;i++)
    k3[i]=f[i];
// k4[i] calc.:
  for(i=1;i<=n;i++)
    y[i]=y0[i]+h*k3[i];
  deriv(y,f);
  if(commontrajexact.index==1) return;
  for(i=1;i<=n;i++)
    k4[i]=f[i];
// New y[i] calculation (i=1,...,n):
  for(i=1;i<=n;i++)
    y[i]=y0[i]+h*(k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.);
  return;
}

///////////////////////////////////////////////////////////////////////


void rungekutta8(int n,double h,double *y)
{
/* General 8th order Runge-Kutta program for the calculation of the
  first order differential equation system
  der_t y_i(t)=f_i(y_1,...,y_n)    /i=1,...n/
    with initial conditions
      y_i(t_0)=y_i^0                     /i=1,...n/.
  Input:
    n: number of equations in the system;
       here n cannot be larger than 100.
    y0=y: n-vector, the initial condition (at starting the subr.)
    h: the "time" step
  Output:
    y: n-vector, the y_i values at "time" t0+h (at the end)
  This subr. uses the subr. deriv, which computes
    the "derivative" function of the system: f[1],...,f[n]
      (f_i, i=1,...,n)
  The y vector have indices i from 1 to n
  (in the calling program it should be defined as y[n+1]),
   but the i=0 values are not used in the calculation.
*/
//////////////////////////////////////////////////////////
// 8th order m, A, a, b value initialization:
  const int m=13;
  const double A[14]={0.,31./720.,0.,0.,0.,0.,16./75.,16807./79200.,
    16807./79200.,243./1760.,0.,0.,243./1760.,31./720.};
  const double a[14]={0.,0.,1./4.,1./12.,1./8.,2./5.,1./2.,6./7.,
    1./7.,2./3.,2./7.,1.,1./3.,1.};
  const double b[14][13]={{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},   // 0

    {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},                          // 1

    {0.,1./4.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},                       // 2

    {0.,5./72.,1./72.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},                  // 3

    {0.,1./32.,0.,3./32.,0.,0.,0.,0.,0.,0.,0.,0.,0.},                  // 4

    {0.,106./125.,0.,-408./125.,352./125.,0.,0.,0.,0.,0.,0.,0.,0.},    // 5

    {0.,1./48.,0.,0.,8./33.,125./528.,0.,0.,0.,0.,0.,0.,0.},           // 6

    {0.,-1263./2401.,0.,0.,39936./26411.,-64125./26411.,5520./2401.,   // 7
     0.,0.,0.,0.,0.,0.},                                               // 7

    {0.,37./392.,0.,0.,0.,1625./9408.,-2./15.,61./6720.,               // 8
     0.,0.,0.,0.,0.},                                                  // 8

    {0.,17176./25515.,0.,0.,-47104./25515.,1325./504.,                 // 9
     -41792./25515.,20237./145800.,4312./6075.,0.,0.,0.,0.},           // 9

    {0.,-23834./180075.,0.,0.,-77824./1980825.,-636635./633864.,       // 10
     254048./300125.,-183./7000.,8./11.,-324./3773.,0.,0.,0.},         // 10

    {0.,12733./7600.,0.,0.,-20032./5225.,456485./80256.,               // 11
     -42599./7125.,339227./912000.,-1029./4180.,1701./1408.,           // 11
     5145./2432.,0.,0.,},                                              // 11

    {0.,-27061./204120.,0.,0.,40448./280665.,-1353775./1197504.,       // 12
     17662./25515.,-71687./1166400.,98./225.,1./16.,3773./11664.,      // 12
     0.,0.},                                                           // 12

    {0.,11203./8680.,0.,0.,-38144./11935.,2354425./458304.,            // 13
     -84046./16275.,673309./1636800.,4704./8525.,                      // 13
     9477./10912.,-1029./992.,0.,729./341.}};                          // 13
// End of m, A, a, b initialization
//  (the A[0], a[0], b[0,j], b[i,0] values are not used !)
///////////////////////////////////////////////////////////
  int i,j,ip,l;
  double sumj,suml;
  double k[101][14],f[101],y0[101],s;
// Starting y values:
  for(i=1;i<=n;i++)
    y0[i]=y[i];
  commontrajexact.indexrungekutta8=1;
// k[i][j] matrix element calculation (i=1,...,n;j=1,...,m):
// j=1:
  deriv(y,f);
  if(commontrajexact.index==1) return;
  for(i=1;i<=n;i++)
    k[i][1]=f[i];
// j>1:
  for(j=2;j<=m;j++)
  {
    for(ip=1;ip<=n;ip++)
    {
      suml=0;
      for(l=1;l<=j-1;l++)
        suml+=b[j][l]*k[ip][l];
      y[ip]=y0[ip]+h*suml;
    }
    s=sqrt((y[1]-y0[1])*(y[1]-y0[1])+(y[2]-y0[2])*(y[2]-y0[2])+(y[3]-y0[3])*(y[3]-y0[3]));
//    printf("s= %12.4f   \t\n",s);
    if(s>1.2*commontrajexact.slimit)
    {
      commontrajexact.indexrungekutta8=0;
      return;
    }
    deriv(y,f);
    if(commontrajexact.index==1) return;
    for(i=1;i<=n;i++)
      k[i][j]=f[i];
  }
// New y[i] calculation (i=1,...,n):
  for(i=1;i<=n;i++)
  { sumj=0;
    for(j=1;j<=m;j++)
      sumj+=A[j]*k[i][j];
    y[i]=y0[i]+h*sumj;
  }
  return;
}

/////////////////////////////////////////////////////////////

