
void trajelectronN(string conclusionfilename,int N,double StartgyraR,ofstream& OutStream);
void trajelectron1(double *x,double *v,int& electronindex);
void velocity(double *n,double kin_energy,double costheta,double phi,double *v);
void electronstart(double z,double r,double costheta,double phidir,double Ekin,double *x,double *v);
void electrongeneration(double *x,double *v, double phidir);
void subrn(double *u,int len);
double randomnumber();
double PhaseAngleFunc(double *v, double *B, int charge);
double ThetaLocal(double theta0, double B0, double Blocal);
double PointDiff(double *x1, double *x2);
void GuidingCenterFunc(double *x, double *v, double *B, double costheta, double gamma, double b_norm, double v_norm, double *GC);

int VectorFindMax(vector<double> &MaxList);
int VectorFindMin(vector<double> &MinList);
int FindTransZone1EndPosMax( double Threshold, vector<double> &MaxList );
int FindTransZone1EndPosMin( double Threshold, vector<double> &MinList );
int FindTransZone2StartPosMax( double Threshold, vector<double> &MaxList );
int FindTransZone2StartPosMin( double Threshold, vector<double> &MinList );
int IJKLRANDOM;


#include "fieldnew.cpp"
#include "trajexact.cpp"

struct typeelectrontraj{
    // Particle properties:
        int  iEkin;             // Energy scaling factor
        double Ekin;            // particle energy [eV]
        double Errenergy;       // Energy error [eV]
        double Time;            // particle flight time [s]
	int hemisphere;
	bool HitApertFlag;
	double StartPhase, ApertPhase, DetPhase, StartPhase2;
	double Apertb, Startb;
	double DetB[4], StartB[4];
	double DeltaTheta;
	double TransitionZone1StartGC[4];
	double TransitionZone1EndGC[4];
	double TransitionZone2StartGC[4];
	double TransitionZone2EndGC[4];
    
	// Parameters of the starting disk = source of the particles representing an appertur or a beam
        double Xstart, ApertGridX;          // x-coordinat of the starting disk center [m]
        double Ystart, ApertGridY;          // y-coordinat of the starting disk center [m]
        double Zstart, ApertGridZ;          // z-coordinat of the starting disk center [m]
        double Rstartmax;       // radius of the starting disk [m]
        double thetaStartmax;   // maximal emission angle of particles in forward direction [rad]
        int theta_fix;          // = 0: theta <= thetaStartmax,  =1: theta = thetaStartmax
	double Startphi;
	double ApertShiftX, ApertShiftY;
    // Detector positions:
        double zdetector;       // axial position of the detector
        double Zend;            // maximal z position
	bool NoMoSOn;
    // NoMoS parameters:
        double R_1;             // radius of the spectrometer
	double alpha;
	bool CompareWBlines;
	bool MonteCarlo;
	bool PointSource;
	bool PercOn;
	bool G4Compare;
	bool Envelope;
	bool ClusterMC;
    // simulation settings:
        int numstepmax;         // maximum number of Runge-Kutta steps
        double Timemax;         // Maximal time [s]
        int typeprint;          // Printing of intermediate quantities  : 1;  no printing:  0
        ofstream outfile;       // path for saving intermediate quantities
        string filepath;        // Path of all data files

  double sizedetector;  // size of the detector in x and y directions
  int step;
  int typeenergy;
  double dstart;
  double pathlength;
  double zstart,rstart,thetastart,excess_energy,Elongstart,Etransstart;
  double phidirstart;
  double Rstart;
  double Deladiabinv; //Delta adiabatic invariant
  int Nstep;
  int pixelring;
}
commonelectrontraj;



/////////////////////////////////////////////////////////

void trajelectronN(string conclusionfilename, int N,double StartgyraR, ofstream& OutStream){     // Trajectory calculation of N number of electrons.

	double PI = 3.1415926;
	int ie;                         // particle counter
	double x[4],v[4], B[4];               // particles spatial and velocity coordinates
	double X,Y,Z;                   // not used
	int i,nel,num0,num1,num2,num3;  // not used
	ofstream outfile;               // path of particle trajectories files
	char numstr1[21],numstr2[21];   // buffers
	string trajfilename;            // file name of particle trajectories files
	int electronindex;              // electronindex: this integer shows the electron exit type
	                                    //  electronindex=0: outside of zonal harmonic convergence
	                                    //  electronindex=1: electron stored
	                                    //  electronindex=2: electron leaves MS
	                                    //  electronindex=3: hit the wall of NoMoS
	                                    //  electronindex=4: reached the detector
	std::stringstream CompareFileName;
	double CrossP[4], CrossPL, GuidX[4], GuidB[4];
	double sign;
	


	// COMPARE WITH BLINES
	if( commonelectrontraj.CompareWBlines ){
		
		outfile.open(conclusionfilename.c_str(),ios::app);

		// Particle loop:
		for(ie=1;ie<=N;ie++){
			
			////////   Writing in CONCLUSION FILE at the start  /////////////////// only for compare with Blines
			if( commonelectrontraj.CompareWBlines ){
			        outfile << commonelectrontraj.ApertGridX << "\t" << commonelectrontraj.ApertGridY;
			        outfile << "\t" << commonelectrontraj.thetaStartmax << "\t" << commonelectrontraj.Ekin;
				outfile << "\t" << commonelectrontraj.Startphi << "\t" << ie;
			}
		
		
			// Detailed tracking filename: definieing the path and the filename, dependent on CompareWBlines
			if(commonelectrontraj.typeprint==1){
			    if( commonelectrontraj.CompareWBlines ) {
			    	CompareFileName << commonelectrontraj.filepath << "trajectories/" << "X" << commonelectrontraj.ApertGridX << "_Y" << commonelectrontraj.ApertGridY;
			    	CompareFileName << "_th" << commonelectrontraj.thetaStartmax << "_Ekin" << int(commonelectrontraj.Ekin/1000.) << "keV_" << "phi"<< commonelectrontraj.Startphi << "_" << ie << ".txt";
			    	CompareFileName >> trajfilename;
			    	commonelectrontraj.outfile.open(trajfilename.c_str(),ios::out);
			    }
			    else{	
			    	sprintf(numstr1, "%d", int(commonelectrontraj.Ekin/1000) );
			    	sprintf(numstr2, "%d", ie);
			    	trajfilename="trajectories/particle";
			    	trajfilename+=numstr2;
			    	trajfilename+="_Ekin";
			    	trajfilename+=numstr1;
			    	trajfilename+="keV.txt";
			    	trajfilename = commonelectrontraj.filepath + trajfilename;
			    	commonelectrontraj.outfile.open(trajfilename.c_str());
			   }
			}
					
			//if you want a random phi, give it here to electrongeneration
			//phidir=2.*PI*randomnumber();          // phi = polar angle

			// Electron generator:
			electrongeneration(x,v,commonelectrontraj.Startphi);            // sets initial condition for x and v, also phi setable
											// uses XStart etc. calcs v from th,phi,Ekin and n from Bfield
		

			/////////////////////////////////////////////////
			// write out first point/starting point to file before going into traj1
        		magfieldtraj(x,B,commonelectrontraj.filepath);                      // Calculate B-field
			// v x B
			if(CHARGE < 0) sign = -1.;
			else sign = 1.;
			CrossP[1] = B[3]*v[2] - B[2]*v[3];
			CrossP[2] = -B[3]*v[1] + B[1]*v[3];
			CrossP[3] = B[2]*v[1] - B[1]*v[2];
			CrossPL = sqrt(CrossP[1]*CrossP[1] + CrossP[2]*CrossP[2] + CrossP[3]*CrossP[3]);
			if(CrossPL == 0.){
				for(int i=1;i<=3;i++){GuidX[i] = x[i];}
			}
			else{
				for(int i=1;i<=3;i++){GuidX[i] = x[i] + sign * CrossP[i]/CrossPL * StartgyraR;};
			}
			
			magfieldtraj(GuidX,GuidB,commonelectrontraj.filepath);

			// write out FIRST POINT traj info
		     	if(commonelectrontraj.typeprint==1){        
        		    // spatial coordinats x[i] of the particle [m]
        		    commonelectrontraj.outfile << x[1] << "\t" << x[2] << "\t" << x[3] << "\t";
        		    // velocity of the particle [m/s]
        		    commonelectrontraj.outfile << v[1] << "\t" << v[2] << "\t" << v[3] << "\t";
			    //magnetic field at particle position
        		    commonelectrontraj.outfile << B[1] << "\t" << B[2] << "\t" << B[3] << "\t";
			    // guiding center
			    commonelectrontraj.outfile << GuidX[1] << "\t" << GuidX[2] << "\t" << GuidX[3] << "\t";
			    // Bfield at guid center
			    commonelectrontraj.outfile << GuidB[1] << "\t" << GuidB[2] << "\t" << GuidB[3] << "\t";
			    // rest of info
			    commonelectrontraj.outfile << commontrajexact.Ekin << setw(16) <<  scientific << 0. << setw(16) << commonelectrontraj.Time << setw(16) << 0. << setw(16) << 0. << endl;
        		}
		
			// Trajectory calculation:
			trajelectron1(x,v,electronindex);
			
			// Writing further in the CONCLUSION FILE at the end: MS index, energy error, storage time
			outfile << "\t" << electronindex << "\t" << commonelectrontraj.Errenergy << "\t" << commonelectrontraj.Time << endl;
			outfile.close();
			
			// Detailed tracking writing: close file
			if(commonelectrontraj.typeprint==1)commonelectrontraj.outfile.close();
		
		}  // end of  electron loop

	} // Compare with blines if end



	//////////////////////////////////////////
	// MONTE CARLO detector data!
	// ///////////////////////////////
	if( commonelectrontraj.MonteCarlo ){

		// Electron generator:
		electrongeneration(x,v,commonelectrontraj.Startphi);            // sets initial condition for x and v, also phi setable
		
		// Trajectory calculation:
		trajelectron1(x,v,electronindex);

		// trajelectron doesnt write out anything -> now write out last point into the right file
		
		OutStream << commonelectrontraj.Xstart << "\t"<< commonelectrontraj.Ystart << "\t" << commonelectrontraj.Zstart << "\t"  << x[1] << "\t" << x[2] << "\t" << x[3] << "\t";
	        OutStream << v[1] << "\t" << v[2] << "\t" << v[3] << "\t" << commontrajexact.Ekin << "\t" << commonelectrontraj.hemisphere;
		OutStream << "\t" << commonelectrontraj.HitApertFlag << endl;

	}




	//////////////////////////////////////////
	// Cluster MONTE CARLO
	// ///////////////////////////////
	if( commonelectrontraj.ClusterMC ){

		// Electron generator:
		electrongeneration(x,v,commonelectrontraj.Startphi);            // sets initial condition for x and v, also phi setable
		commonelectrontraj.Startb = sqrt(commonelectrontraj.StartB[1]*commonelectrontraj.StartB[1] +commonelectrontraj.StartB[2]*commonelectrontraj.StartB[2] +commonelectrontraj.StartB[3]*commonelectrontraj.StartB[3]);
		commonelectrontraj.StartPhase = commonelectrontraj.Startphi+M_PI;// investigating velocity(), we see that phi of v is just a phase offset to phase of starting position	
		// StartB is written in electronstart, to be used as cross check to PhaseAngleFunc with Bfield vector
		commonelectrontraj.StartPhase2 = PhaseAngleFunc(v, commonelectrontraj.StartB, -1);	

		// Trajectory calculation:
		trajelectron1(x,v,electronindex);

		// we write out StartP, Ekin, StartPhase (2 methods), BDV	
		OutStream << commonelectrontraj.Xstart << "\t"<< commonelectrontraj.Ystart << "\t" << commonelectrontraj.Zstart << "\t";
	        OutStream << commontrajexact.Ekin << "\t" << commonelectrontraj.StartPhase << "\t" << commonelectrontraj.StartPhase2 << "\t" << commonelectrontraj.Startb << "\t";

		//reuse Startb for norm at Detector
		commonelectrontraj.Startb = sqrt(commonelectrontraj.DetB[1]*commonelectrontraj.DetB[1] +commonelectrontraj.DetB[2]*commonelectrontraj.DetB[2] +commonelectrontraj.DetB[3]*commonelectrontraj.DetB[3]);
		//
		// we write away ApertPoint, Apert Phase, B_A
		OutStream << commonelectrontraj.ApertGridX << "\t"<< commonelectrontraj.ApertGridY << "\t" << commonelectrontraj.ApertGridZ << "\t";
		OutStream << commonelectrontraj.ApertPhase << "\t" << commonelectrontraj.Apertb << "\t";

		//we write away transition zone 1 data: Start GC, EndGC
		OutStream << commonelectrontraj.TransitionZone1StartGC[1] << "\t" << commonelectrontraj.TransitionZone1StartGC[2] << "\t" << commonelectrontraj.TransitionZone1StartGC[3] << "\t";
		OutStream << commonelectrontraj.TransitionZone1EndGC[1] << "\t" << commonelectrontraj.TransitionZone1EndGC[2] << "\t" << commonelectrontraj.TransitionZone1EndGC[3] << "\t";

		//we write away transition zone 2 data: Start GC, EndGC
		OutStream << commonelectrontraj.TransitionZone2StartGC[1] << "\t" << commonelectrontraj.TransitionZone2StartGC[2] << "\t" << commonelectrontraj.TransitionZone2StartGC[3] << "\t";
		OutStream << commonelectrontraj.TransitionZone2EndGC[1] << "\t" << commonelectrontraj.TransitionZone2EndGC[2] << "\t" << commonelectrontraj.TransitionZone2EndGC[3] << "\t";

		//we write away detector data: DetPoint, DetPhase, B_Det
		OutStream << x[1] << "\t" << x[2] << "\t" << x[3] << "\t";
		OutStream << PhaseAngleFunc(v, commonelectrontraj.DetB, -1) << "\t" << commonelectrontraj.Startb;

		//we don't know theta out here, so we calc it from v and the reference vector (theta relative to detector)
		double V=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
		double vpara = v[2]*-sin(commonelectrontraj.alpha/180.*M_PI) + v[3]*cos(commonelectrontraj.alpha/180.*M_PI);
		double thetaend = acos(vpara/V);
		OutStream << "\t" << thetaend/M_PI*180. << endl;

	}






	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// COMPARE WITH G4
	if( commonelectrontraj.G4Compare ){
		
		outfile.open(conclusionfilename.c_str(),ios::app);

		// Particle loop:
		for(ie=1;ie<=N;ie++){
			
			////////   Writing in CONCLUSION FILE at the start  /////////////////// only for compare with Blines
		        outfile << commonelectrontraj.ApertGridX << "\t" << commonelectrontraj.ApertGridY;
		        outfile << "\t" << commonelectrontraj.thetaStartmax << "\t" << commonelectrontraj.Ekin;
			outfile << "\t" << commonelectrontraj.Startphi;
			
		
		
			// Detailed tracking filename: definieing the path and the filename, dependent on CompareWBlines
			if(commonelectrontraj.typeprint==1){
			    if( commonelectrontraj.G4Compare ) {
			    	CompareFileName << commonelectrontraj.filepath << "trajectories/" << "X" << commonelectrontraj.ApertGridX << "_Y" << commonelectrontraj.ApertGridY;
			    	CompareFileName << "_th" << commonelectrontraj.thetaStartmax << "_Ekin" << int(commonelectrontraj.Ekin/1000.) << "keV_" << "phi"<< commonelectrontraj.Startphi << "_" << ie << ".txt";
			    	CompareFileName >> trajfilename;
			    	commonelectrontraj.outfile.open(trajfilename.c_str(),ios::out);
			    }
			    else{	
			    	sprintf(numstr1, "%d", (int)commonelectrontraj.Ekin/1000);
			    	sprintf(numstr2, "%d", ie);
			    	trajfilename="trajectories/particle";
			    	trajfilename+=numstr2;
			    	trajfilename+="_Ekin";
			    	trajfilename+=numstr1;
			    	trajfilename+="keV.txt";
			    	trajfilename = commonelectrontraj.filepath + trajfilename;
			    	commonelectrontraj.outfile.open(trajfilename.c_str());
			   }
			}
					
			
			// instead of electron generator, where vector is hardcoded to B-field, we define x and v outside here
			// Electron generator:
			//electrongeneration(x,v,commonelectrontraj.Startphi);            // sets initial condition for x and v, also phi setable
			x[1] = 0.;
			x[2] = 0.;
			x[3] = 1.;
			velocity(x,commonelectrontraj.Ekin,cos(commonelectrontraj.thetaStartmax),commonelectrontraj.Startphi,v); // x is used here for normal vector
			
			x[1] = commonelectrontraj.Xstart;
			x[2] = commonelectrontraj.Ystart;
			x[3] = commonelectrontraj.Zstart;


			/////////////////////////////////////////////////
			// write out first point/starting point to file before going into traj1
        		magfieldtraj(x,B,commonelectrontraj.filepath);                      // Calculate B-field
			// v x B
			CrossP[1] = B[3]*v[2] - B[2]*v[3];
			CrossP[2] = -B[3]*v[1] + B[1]*v[3];
			CrossP[3] = B[2]*v[1] - B[1]*v[2];
			CrossPL = sqrt(CrossP[1]*CrossP[1] + CrossP[2]*CrossP[2] + CrossP[3]*CrossP[3]);
			if(CrossPL < 1.e-12 ){
				for(int i=1;i<=3;i++){GuidX[i] = x[i];}
			}
			else{
				for(int i=1;i<=3;i++){GuidX[i] = x[i] - CrossP[i]/CrossPL * StartgyraR;};
			}
			
			magfieldtraj(GuidX,GuidB,commonelectrontraj.filepath);

			// write out FIRST POINT traj info
		     	if(commonelectrontraj.typeprint==1){        
        		    // spatial coordinats x[i] of the particle [m]
        		    commonelectrontraj.outfile << x[1] << "\t" << x[2] << "\t" << x[3] << "\t";
        		    // velocity of the particle [m/s]
        		    commonelectrontraj.outfile << v[1] << "\t" << v[2] << "\t" << v[3] << "\t";
			    //magnetic field at particle position
        		    commonelectrontraj.outfile << B[1] << "\t" << B[2] << "\t" << B[3] << "\t";
			    // guiding center
			    commonelectrontraj.outfile << GuidX[1] << "\t" << GuidX[2] << "\t" << GuidX[3] << "\t";
			    // Bfield at guid center
			    commonelectrontraj.outfile << GuidB[1] << "\t" << GuidB[2] << "\t" << GuidB[3] << "\t";
			    // rest of info
			    commonelectrontraj.outfile << commontrajexact.Ekin << setw(16) <<  scientific << 0. << setw(16) << commonelectrontraj.Time << setw(16) << 0. << setw(16) << 0. << endl;
        		}
		

			// Trajectory calculation:
			trajelectron1(x,v,electronindex);
			
			
			// Writing further in the CONCLUSION FILE at the end: MS index, energy error, storage time
			outfile << "\t" << commonelectrontraj.ApertGridX << "\t" << commonelectrontraj.ApertGridY << "\t" << commonelectrontraj.ApertGridZ;
			outfile << "\t" << x[1] << "\t" << x[2] << "\t" << x[3] << "\t" << electronindex << "\t" << commonelectrontraj.Errenergy << "\t" << commonelectrontraj.Time << endl;
			outfile.close();
			
			// Detailed tracking writing: close file
			if(commonelectrontraj.typeprint==1)commonelectrontraj.outfile.close();
		
		}  // end of  electron loop

	} // Compare with G4 if end
	

	//////////////////////////////
	// Point source
	if( commonelectrontraj.PointSource ){

		// Electron generator:
		electrongeneration(x,v,commonelectrontraj.Startphi);            // sets initial condition for x and v, also phi setable
		
		// Trajectory calculation:
		trajelectron1(x,v,electronindex);

		// trajelectron doesnt write out anything -> now write out last point into the right file
		
		OutStream << commonelectrontraj.Xstart << "\t"<< commonelectrontraj.Ystart << "\t" << commonelectrontraj.Zstart;
	        OutStream << "\t" << commonelectrontraj.ApertGridX << "\t" << commonelectrontraj.ApertGridY << "\t" << commonelectrontraj.ApertGridZ;
		OutStream << "\t" << x[1] << "\t" << x[2] << "\t" << x[3] << "\t" << commontrajexact.Ekin << endl;

	}

	return;
}

//////////////////////////////////////////////////////////



void trajelectron1(double *x, double *v, int& electronindex){
	// Trajectory calculation of 1 electron, with exact electron motion.
	// x[1],x[2],x[3]: starting and final spatial coordinates of the electron (m)
	// v[1],v[2],v[3]: starting and final velocity components of the electron (m/s)
	
	
	// used variables:
	double PI = 3.1415926;
	const double c=299792458.;  // velocity of light in SI units
	int istep, i;               // loop counter (Runge-Kutta, coordinates)
	double Phi,B[4],b,btheta;   // Electro-magnetic parameters
	double gam,p[4],costheta,theta;  // Particle parameters (absolute velocity[m/s],[m²/s²],...)
	double err,energy0;         // Energy error calculation
	char filename[20];
	double z,r,rmax,pathlength,Ef0[4],Phi0;
	double E0,E,kin_energy,phi,nH2,zend,vz,vzprev,d;
	int j,index,numstepmax;
	double rc,energy;
	double Elong,plong2,pbper,pb;
	double pperp2,adiabinv,adiabinv0,deladiabinv,maxdeladiabinv;
	double zmax,Elongmin,zmin,Rmean;
	double gyraR, CrossPL;
	double GuidX[4], CrossP[4], GuidB[4];
	long double V, V_parallel, V_perpend;
	int ZatApert = 0;
	double sign;

	// variables for Transition Zone Investigations
	double ThetaAdiab;
	//bool TransitionZone1Start = false;
	//bool TransitionZone1End = false;
	//bool TransitionZone2Start = false;
	//bool TransitionZone2End = false;
	//double DeltaThetaExtremalOld = 0., DeltaThetaExtremalNew = 0.;
	double DeltaThetaOld = 0.;
	int DeltaThetaSignNew = 1, DeltaThetaSignOld = 1;
	//double ExtremalSlope_old = 0.;
	//double ExtremalSlope_new = 0.0001; 
	double Zone1StartThreshold_rel = 0.8;
	double Zone1EndThreshold_rel = 0.2;
	double Zone2StartThreshold_rel = 0.2;
	double Zone2EndThreshold_rel = 0.8;
	int POIMin, POIMax; // point of interest , reused for all the zone boundaries
	//double xExtremal_old[4];
	//xExtremal_old[1] = 0.;
	//xExtremal_old[2] = 1.;
	//xExtremal_old[3] = -2.;
	
	// vectors for storing the extrema including the corresponding guiding centers for the Transition zones in Cluster MC
	vector<double> MaximaValues;
	vector<double> MinimaValues;
	vector<double> MaxGCX;
	vector<double> MaxGCY;
	vector<double> MaxGCZ;
	vector<double> MinGCX;
	vector<double> MinGCY;
	vector<double> MinGCZ;
	

	/////////
	// initialize calculation:
	commontrajexact.trajstart=0;
	commonelectrontraj.Errenergy=0.;    // energy error [eV]
	commonelectrontraj.Time=0.;         // particle flight time [s]
	
	// calc the end of RxB for arbitrary ALPHA
	double zRxBEnd = commonelectrontraj.R_1*cos((commonelectrontraj.alpha-90.)/180.*M_PI);
	double yRxBEnd = -commonelectrontraj.R_1*sin((commonelectrontraj.alpha-90.)/180.*M_PI);
	double zdetwAlpha = zRxBEnd + commonelectrontraj.zdetector*sin((commonelectrontraj.alpha-90.)/180.*M_PI);
	double ydetwAlpha = yRxBEnd + commonelectrontraj.zdetector*cos((commonelectrontraj.alpha-90.)/180.*M_PI);
	//cout << "TRAJBUG: zdetwAlpha = " << zdetwAlpha << ydetwAlpha << endl;

	commontrajexact.filepath = commonelectrontraj.filepath;
	// Starting of Runge-Kutta loop:
	for(istep=1; istep<=commonelectrontraj.numstepmax; istep++){
		
		// Runge-Kutta step:
		trajexact(x,v);             // calculation of the x and v after the Runge-Kutta step
		if(commontrajexact.index==1){
		    electronindex=0;
		    break;
		}
	
		commonelectrontraj.Time+=commontrajexact.h;     // Advance the flight time clock
		// calculate parameters, that are required most of the time!
		V=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
		Phi=commontrajexact.Phi;                // Electric potential
		magfieldtraj(x,B,commonelectrontraj.filepath);                      // Calculate B-field
		b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // absolute B-field [T]
		if(x[3]>=0.) btheta = acos((B[3]*x[2] - B[2]*x[3])/(b*sqrt(x[2]*x[2]+x[3]*x[3])))*180./PI;
		else btheta = acos(B[3]/b)*180./PI;
		V_parallel = (v[1]*B[1]+v[2]*B[2]+v[3]*B[3])/b;
		costheta=V_parallel/V;
		// case if argument of acos is so close to 1, that error makes acos give NAN				
		if( fabs(costheta - 1.) < 10.e-12 ){
			costheta = 1.;
		}
		ThetaAdiab = ThetaLocal(commonelectrontraj.thetaStartmax, commonelectrontraj.Startb, b)/M_PI*180.; //in degree
		theta=180./PI*acos(costheta);         // [deg]
		gam=1./sqrt(1.-V*V/(c*c));     // gamma factor:
		for(i=1;i<=3;i++)p[i]=MASS*v[i]*gam;  // p: relativistic momentum
		
		// The longitudinal and transversal kinetic energies(in unit eV):
		pb=p[1]*B[1]+p[2]*B[2]+p[3]*B[3];
		pbper=pb/(b*b);
		plong2=0.;
		for(i=1;i<=3;i++)  plong2+=(pbper*B[i])*(pbper*B[i]);
		Elong=plong2/(1.+gam)/MASS/fabs(CHARGE);
		pperp2=0.;
		for(i=1;i<=3;i++) pperp2+=(p[i]-pbper*B[i])*(p[i]-pbper*B[i]);
		adiabinv=pperp2/(2.*MASS*b*fabs(CHARGE));
		if(istep==1) adiabinv0=adiabinv;
		deladiabinv=(adiabinv-adiabinv0)/adiabinv0;
		
		// Calculation of the energy error:
		if(istep==1)energy0=commontrajexact.energy;   // initilization
		err=(commontrajexact.energy-energy0)/energy0; // calculation of energy error
		if(fabs(err)>commonelectrontraj.Errenergy)commonelectrontraj.Errenergy=fabs(err);     // memorize the larges error
	
		// Z, R: z and r coordinates in m:  << setw(10) << Z << (10) << R
		// b: magnetic field absolute value in T
		// Phi: electric potential in V
		// Ekin: kinetic energy in eV
		// err: relative change of total energy
		// commonelectrontraj.Time: time in s



		////////////////////////////////////////////////////////////////////////////////////////
		// BREAK CONDITIONS OF LOOP
		//////////////////////////////////////////////////////////////////////////////////////
		// Time break
		if(istep==commonelectrontraj.numstepmax || commonelectrontraj.Time>commonelectrontraj.Timemax){
		    	cout << "TRAP BREAK" << endl;
			electronindex=1;        // particle is trapped in the experment
		    break;                  // stop trajacking because of steps limit or timelimit
		}
		
		// Z Outside break
		if(x[3] >  commonelectrontraj.Zend ){   // break if particle reaches the seted end or fall behind the start position
		   	cout << "Z Zend BREAK" << endl;
			electronindex=2;        // particle left the
			break;
		} 
		
		// particle was reflected
		if( x[3] <  commonelectrontraj.Zstart - 0.1 && x[2] > 0. ){
		   	cout << "Z Zstart BREAK" << endl;
			electronindex=3;        // particle left the
			break;
		}
		
		//Nomos break, assuming a 300 tube radius!
		//if(x[3]>0. && pow((sqrt(x[2]*x[2]+x[3]*x[3])-commonelectrontraj.R_1),2)+x[1]*x[1] > 0.15*0.15){
		//	cout << "NOMOS WALL BREAK" << endl;
		//	electronindex=4;        // particle hit the RxB tube
		//	break;                  // stop trajacking
		//}
		
		// Detector break (includes PERConly, detector set differently, without y restriction)
		if(commonelectrontraj.NoMoSOn == true){
			// Nomos detector break
			// we distinguish between 180, 90 ° and inbetween
			if(commonelectrontraj.alpha == 180.){
				if(x[3] <= commonelectrontraj.zdetector && x[2] < 0.){
					electronindex=5;        // particle reached the detector
					cout << "NOMOS det reached 180" << endl;
					commonelectrontraj.DetB[1] = B[1];
					commonelectrontraj.DetB[2] = B[2];
					commonelectrontraj.DetB[3] = B[3];
					break;                  // stop trajacking
				}
			}
			else if(commonelectrontraj.alpha == 90.){
				if( x[2] <= commonelectrontraj.zdetector){
					electronindex=5;        // particle reached the detector
					cout << "NOMOS det reached 90" << endl;
					commonelectrontraj.DetB[1] = B[1];
					commonelectrontraj.DetB[2] = B[2];
					commonelectrontraj.DetB[3] = B[3];
					break;                  // stop trajacking
				}
			}
			else{
				if(x[3] <= zdetwAlpha && x[2] <= ydetwAlpha){
					electronindex=5;        // particle reached the detector
					cout << "NOMOS det reached ??" << endl;
					break;                  // stop trajacking
				}
			}
		}
		else{
			// PERC detector break
			if( x[3] >  commonelectrontraj.zdetector ){
				electronindex=5;        // particle reached the detector
				cout << "PERC det reached" << endl;
				break;                  // stop trajacking
			}
		}
		
		// break if hit the wall in Perc
		if(x[3] <  -3. && (fabs(x[2] - 1) >  0.5 || fabs(x[1]) > 0.5)){   
		    	cout << "PERC WALL BREAK" << endl;
			electronindex=6;        // particle left the
		    	break;
		}
		



		////////////////////////////////////////////////////////////////////////////////////////
		// SECTION ABOUT WRITING DATA AWAY
		///////////////////////////////////////////////////////////////////////////////////////
		// write away full track
		if( commonelectrontraj.typeprint==1 ) {

			// calculate guiding center and evaluate guiding center B-field as well
			V_perpend = sin(acos(costheta)) * V;
			if(CHARGE < 0) sign = -1.;
			else sign = 1.;
			gyraR = gam * MASS * V_perpend/ 1.602177e-19 / b;
			CrossP[1] = B[3]*v[2] - B[2]*v[3];
			CrossP[2] = -B[3]*v[1] + B[1]*v[3];
			CrossP[3] = B[2]*v[1] - B[1]*v[2];
			CrossPL = sqrt(CrossP[1]*CrossP[1] + CrossP[2]*CrossP[2] + CrossP[3]*CrossP[3]);
			if(CrossPL == 0.){
				for(int i=1;i<=3;i++){GuidX[i] = x[i];}
			}
			else{
				for(int i=1;i<=3;i++){GuidX[i] = x[i] + sign* CrossP[i]/CrossPL * gyraR;};
			}
			magfieldtraj(GuidX,GuidB,commonelectrontraj.filepath);

			// write out traj info
        		if(commonelectrontraj.typeprint==1){        // save simulation data of the particle track
        		    // spatial coordinats x[i] of the particle [m]
        		    commonelectrontraj.outfile << x[1] << "\t" << x[2] << "\t" << x[3] << "\t";
        		    // velocity of the particle [m/s]
        		    commonelectrontraj.outfile << v[1] << "\t" << v[2] << "\t" << v[3] << "\t";
			    //magnetic field at particle position
        		    commonelectrontraj.outfile << B[1] << "\t" << B[2] << "\t" << B[3] << "\t";
			    // guiding center
			    commonelectrontraj.outfile << GuidX[1] << "\t" << GuidX[2] << "\t" << GuidX[3] << "\t";
			    // Bfield at guid center
			    commonelectrontraj.outfile << GuidB[1] << "\t" << GuidB[2] << "\t" << GuidB[3] << "\t";
			    // rest of info
			    commonelectrontraj.outfile << commontrajexact.Ekin << setw(16) <<  scientific << err << setw(16) << commonelectrontraj.Time << setw(16) << btheta << setw(16) << theta << endl;
        		}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////
		// NORMAL MONTE CARLO
		// in case of normal Monte Carlo, we write away at aperture
		//	first close point to aperture
		if( commonelectrontraj.MonteCarlo ){ //for monte carlo, check if particle went through aperture or not
			if( x[3] >= -0.3 && ZatApert == 0 ){ //for MonteCarlo, we get Aperture size by the global values ApertGridX/Y
				ZatApert = 1;
				
				//we check, if we got through the aperture, if not, we flag and break
				// ApertGridX and Y are reused for MC as the size of the aperture
				if	( 
					commonelectrontraj.ApertShiftX - commonelectrontraj.ApertGridX/2. < x[1] < commonelectrontraj.ApertShiftX + commonelectrontraj.ApertGridX/2. &&
					commonelectrontraj.R_1 + commonelectrontraj.ApertShiftY - commonelectrontraj.ApertGridY/2. < x[2] < 
					commonelectrontraj.R_1 + commonelectrontraj.ApertShiftY + commonelectrontraj.ApertGridY/2. 
					){
				commonelectrontraj.HitApertFlag = false;
				}
				else {
					commonelectrontraj.HitApertFlag = true;
					break;
				}
			
			}
		}


		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// for CLUSTER MC 
		if( commonelectrontraj.ClusterMC ){ 
			
			/////////////////////////////////////////////
			//writes away the pos at apert for manual cut later
			if( x[3] >= -0.3 && ZatApert == 0 ){ 
				ZatApert = 1;
				
				// ApertGridX and Y are reused for MC as the size of the aperture
				commonelectrontraj.ApertGridX = x[1];
				commonelectrontraj.ApertGridY = x[2];
				commonelectrontraj.ApertGridZ = x[3];
				commonelectrontraj.Apertb = b;
				commonelectrontraj.ApertPhase = PhaseAngleFunc(v, B, -1);
			}
			
			///////////////////////////////////////////////////////////////
			//check for transition zones. start search after aperture
			if( ZatApert == 1 ){

				// for every check, we calc the new Delta THETA
				commonelectrontraj.DeltaTheta = ThetaAdiab - theta; // in deg

				
				
				
				
				///////////////////////////////////////////////////////
				//ALTERNATIVE APPROACH just with storing all extrema and then compare to known plateau

				//we search for extrema
				// check the change of Delta Theta, if positive or negative
				if( DeltaThetaOld <= commonelectrontraj.DeltaTheta ) DeltaThetaSignNew = 1;
				else DeltaThetaSignNew = -1;

				// if the old sign is positive and new sign is negative -> MaximumPoint
				if( DeltaThetaSignNew == -1 && DeltaThetaSignOld == 1 ){
					// save maximum value
					MaximaValues.push_back(commonelectrontraj.DeltaTheta);
					// calc guiding center and store
					GuidingCenterFunc(x, v, B, costheta, gam, b, V, GuidX);
					MaxGCX.push_back(GuidX[1]);
					MaxGCY.push_back(GuidX[2]);
					MaxGCZ.push_back(GuidX[3]);
				}
				// if the old sign is neg. and new sign is pos -> MinPoint
				else if( DeltaThetaSignNew == 1 && DeltaThetaSignOld == -1 ){
					MinimaValues.push_back(commonelectrontraj.DeltaTheta);
					GuidingCenterFunc(x, v, B, costheta, gam, b, V, GuidX);
					MinGCX.push_back(GuidX[1]);
					MinGCY.push_back(GuidX[2]);
					MinGCZ.push_back(GuidX[3]);
				}
				// after storing all those vectors, after the runge kutta loop, we extract the correct guiding centers with functions defined on the very bottom

				// update the old values
				DeltaThetaOld = commonelectrontraj.DeltaTheta;
				DeltaThetaSignOld = DeltaThetaSignNew;



				////////////////////////////////////////////////////////////////////
				//old Variant

				/*

				// before TransitionZone 1
				if( !TransitionZone1Start ){
					
					if( abs(commonelectrontraj.DeltaTheta) > Zone1LowThreshold_absolute ){
						TransitionZone1Start = true;
						// now we want the Guiding Center of the starting point
						V_perpend = sin(acos(costheta)) * V;
						if(CHARGE < 0) sign = -1.;
						else sign = 1.;
						gyraR = gam * MASS * V_perpend/ 1.602177e-19 / b;
						// v x B
						CrossP[1] = B[3]*v[2] - B[2]*v[3];
						CrossP[2] = -B[3]*v[1] + B[1]*v[3];
						CrossP[3] = B[2]*v[1] - B[1]*v[2];
						CrossPL = sqrt(CrossP[1]*CrossP[1] + CrossP[2]*CrossP[2] + CrossP[3]*CrossP[3]);
						if( CrossPL < 1e-12 ){
							for(int i=1;i<=3;i++){commonelectrontraj.TransitionZone1StartGC[i] = x[i];}
						}
						else{
							for(int i=1;i<=3;i++){commonelectrontraj.TransitionZone1StartGC[i] = x[i] + sign* CrossP[i]/CrossPL * gyraR;};
						}

					}
				}




				// now the case after TransitionZone 1 has started -> we need to save the extremal points of the oscillation
				// then compare extremal points -> if slope between two extremal points is below a threshold -> Transition Zone end.
				// FOR NOW, we check only maxima! (because minima dont followsame envelope necessarily!)
				// start search after Transition Zone 1 has started, but has not ended yet!
				// we can hardcode that the search for the transitionzone1End should only start after z > 0, we may safe some search time!
				if( TransitionZone1Start && !TransitionZone1End && x[3] > 0. ){
					
					// check the change of Delta Theta, if positive or negative
					if( DeltaThetaOld <= commonelectrontraj.DeltaTheta ) DeltaThetaSignNew = 1;
					else DeltaThetaSignNew = -1;

					// if the old sign is positive and new sign is negative -> MaximumPoint
					if( DeltaThetaSignNew == -1 && DeltaThetaSignOld == 1 ){ //there is a new extremal point (the last point is the extremal point, not the new one
						
						//first we have to check, if we have already found an extremal point before. If not, the x_old distance to the first extremal point
						//strongly suppresses the slope value and therefore would immediately give a small slope and the Zone End!
						//so at first extremal point, set x_old and do nothing else
						if(DeltaThetaExtremalNew == 0.){ // if this is zero, there was no extremal point so far
							//so we just set the old extremal point as this one for a start
							xExtremal_old[1] = x[1];	
							xExtremal_old[2] = x[2];	
							xExtremal_old[3] = x[3];	

							// update the extremal values
							DeltaThetaExtremalNew = DeltaThetaOld;
							//still, search must continue, so we make the old values current
							DeltaThetaOld = commonelectrontraj.DeltaTheta;	
							DeltaThetaSignOld = DeltaThetaSignNew;
							DeltaThetaExtremalOld = DeltaThetaExtremalNew;

						}
						else{

							DeltaThetaExtremalNew = DeltaThetaOld;
							//still, search must continue, so we make the old values current
							DeltaThetaOld = commonelectrontraj.DeltaTheta;	
							DeltaThetaSignOld = DeltaThetaSignNew;

							//now we have to check the slope. To be accurate, we need to calc the positional difference for division
							Extremal_slope =( DeltaThetaExtremalNew - DeltaThetaExtremalOld ) / PointDiff( xExtremal_old, x );
							// now we can also update the old extremal point
							DeltaThetaExtremalOld = DeltaThetaExtremalNew;
							xExtremal_old[1] = x[1];	
							xExtremal_old[2] = x[2];	
							xExtremal_old[3] = x[3];	

							if(Extremal_slope < Zone1HighThreshold_slope){ // END OF ZONE!
								
								TransitionZone1End = true;
								
								// and calculate the guiding center of that endpoint
								V_perpend = sin(acos(costheta)) * V;
								if(CHARGE < 0) sign = -1.;
								else sign = 1.;
								gyraR = gam * MASS * V_perpend/ 1.602177e-19 / b;
								// v x B
								CrossP[1] = B[3]*v[2] - B[2]*v[3];
								CrossP[2] = -B[3]*v[1] + B[1]*v[3];
								CrossP[3] = B[2]*v[1] - B[1]*v[2];
								CrossPL = sqrt(CrossP[1]*CrossP[1] + CrossP[2]*CrossP[2] + CrossP[3]*CrossP[3]);
								if( CrossPL < 1e-12 ){
									for(int i=1;i<=3;i++){commonelectrontraj.TransitionZone1EndGC[i] = x[i];}
								}
								else{
									for(int i=1;i<=3;i++){commonelectrontraj.TransitionZone1EndGC[i] = x[i] + sign* CrossP[i]/CrossPL * gyraR;};
								}
						
							} //endif for END OF ZONE1
						} // end of else, when the extremal point was not the first anymore
					} // end of if for extremal point
					// Otherwise not
					else{ // no new extremal point
						DeltaThetaOld = commonelectrontraj.DeltaTheta;
						DeltaThetaSignOld = DeltaThetaSignNew;
					}

				}



				/////////////////
				// after Zone 1, along the envelope plateau, search start of 
				// Zone 2 hasnt started yet
				if( TransitionZone1End && !TransitionZone2Start ){
					
					//analogeous to end of Zone1, we have to check change of extrema
					// check the change of Delta Theta, if positive or negative
					if( DeltaThetaOld <= commonelectrontraj.DeltaTheta ) DeltaThetaSignNew = 1;
					else DeltaThetaSignNew = -1;


					// if the old sign is positive and new sign is negative -> MaximumPoint
					if( DeltaThetaSignNew == -1 && DeltaThetaSignOld == 1 ){ //there is a new extremal point (the last point is the extremal point, not the new one
						
						DeltaThetaExtremalNew = DeltaThetaOld;
						//still, search must continue, so we make the old values current
						DeltaThetaOld = commonelectrontraj.DeltaTheta;	
						DeltaThetaSignOld = DeltaThetaSignNew;

						//now we have to check the slope. To be accurate, we need to calc the positional difference for division
						Extremal_slope =( DeltaThetaExtremalNew - DeltaThetaExtremalOld ) / PointDiff( xExtremal_old, x );
						// now we can also update the old extremal point
						DeltaThetaExtremalOld = DeltaThetaExtremalNew;
						xExtremal_old[1] = x[1];	
						xExtremal_old[2] = x[2];	
						xExtremal_old[3] = x[3];	

						if(Extremal_slope < Zone2HighThreshold_slope){ // END OF ZONE!
							
							TransitionZone2Start = true;
							
							// and calculate the guiding center of that endpoint
							V_perpend = sin(acos(costheta)) * V;
							if(CHARGE < 0) sign = -1.;
							else sign = 1.;
							gyraR = gam * MASS * V_perpend/ 1.602177e-19 / b;
							// v x B
							CrossP[1] = B[3]*v[2] - B[2]*v[3];
							CrossP[2] = -B[3]*v[1] + B[1]*v[3];
							CrossP[3] = B[2]*v[1] - B[1]*v[2];
							CrossPL = sqrt(CrossP[1]*CrossP[1] + CrossP[2]*CrossP[2] + CrossP[3]*CrossP[3]);
							if( CrossPL < 1e-12 ){
								for(int i=1;i<=3;i++){commonelectrontraj.TransitionZone2StartGC[i] = x[i];}
							}
							else{
								for(int i=1;i<=3;i++){commonelectrontraj.TransitionZone2StartGC[i] = x[i] + sign* CrossP[i]/CrossPL * gyraR;};
							}
						
						} //endif for END OF ZONE1
					} // end of if for extremal point
					// Otherwise not
					else{ // no new extremal point
						DeltaThetaOld = commonelectrontraj.DeltaTheta;
						DeltaThetaSignOld = DeltaThetaSignNew;
					}
				} // end of if for start of transition Zone 2



				/////////////////////////////
				// Zone 2 has started, but not ended yet
				// we cant just simply look at delta theta, because delta theta oscillates. we still have to check the extrema and go below a threshold!
				if( TransitionZone2Start && !TransitionZone2End ){

					//analogeous to end of Zone1, we have to check change of extrema
					// check the change of Delta Theta, if positive or negative
					if( DeltaThetaOld <= commonelectrontraj.DeltaTheta ) DeltaThetaSignNew = 1;
					else DeltaThetaSignNew = -1;


					// if the old sign is positive and new sign is negative -> MaximumPoint
					if( DeltaThetaSignNew == -1 && DeltaThetaSignOld == 1 ){ //there is a new extremal point (the last point is the extremal point, not the new one
						
						DeltaThetaExtremalNew = DeltaThetaOld;
						//still, search must continue, so we make the old values current
						DeltaThetaOld = commonelectrontraj.DeltaTheta;	
						DeltaThetaSignOld = DeltaThetaSignNew;

						if( DeltaThetaExtremalNew < Zone2LowThreshold_absolute ){ // END OF ZONE2!
							
							TransitionZone2End = true;
							
							// and calculate the guiding center of that endpoint
							V_perpend = sin(acos(costheta)) * V;
							if(CHARGE < 0) sign = -1.;
							else sign = 1.;
							gyraR = gam * MASS * V_perpend/ 1.602177e-19 / b;
							// v x B
							CrossP[1] = B[3]*v[2] - B[2]*v[3];
							CrossP[2] = -B[3]*v[1] + B[1]*v[3];
							CrossP[3] = B[2]*v[1] - B[1]*v[2];
							CrossPL = sqrt(CrossP[1]*CrossP[1] + CrossP[2]*CrossP[2] + CrossP[3]*CrossP[3]);
							if( CrossPL < 1e-12 ){
								for(int i=1;i<=3;i++){commonelectrontraj.TransitionZone2EndGC[i] = x[i];}
							}
							else{
								for(int i=1;i<=3;i++){commonelectrontraj.TransitionZone2EndGC[i] = x[i] + sign* CrossP[i]/CrossPL * gyraR;};
							}
						
						} //endif for END OF ZONE1
					} // end of if for extremal point
					// Otherwise not
					else{ // no new extremal point
						DeltaThetaOld = commonelectrontraj.DeltaTheta;
						DeltaThetaSignOld = DeltaThetaSignNew;
					}

				} // end of if search for Zone 2 End


				*/
			

			} // end of if for after aperture

		} // end of if for CLUSTER MC




		//	first close point to aperture for G4 compare and Pitch effects
		if( commonelectrontraj.G4Compare ){ //for monte carlo, check if particle went through aperture or not
			if( x[3] >= -0.29 && ZatApert == 0 ){ //for MonteCarlo, we get Aperture size by the global values ApertGridX/Y
				ZatApert = 1;
				
				// we write the local coordinate to a free variable that we can write to the conclusion file
				commonelectrontraj.ApertGridX = x[1];	
				commonelectrontraj.ApertGridY = x[2];	
				commonelectrontraj.ApertGridZ = x[3];	
			}
		}


		// Point Source Aperture write
		//	first close point to aperture for G4 compare and Pitch effects
		if( commonelectrontraj.PointSource ){ //for monte carlo, check if particle went through aperture or not
			if( x[3] >= -0.29 && ZatApert == 0 ){ //for MonteCarlo, we get Aperture size by the global values ApertGridX/Y
				ZatApert = 1;
				
				// we write the local coordinate to a free variable that we can write to the conclusion file
				commonelectrontraj.ApertGridX = x[1];	
				commonelectrontraj.ApertGridY = x[2];	
				commonelectrontraj.ApertGridZ = x[3];	
			}
		}



	} // end of runge kutta loop



	// after the runge kutta loop, for Cluster MC, we analyse the extrema vectors for the transition zones
	if(commonelectrontraj.ClusterMC){

		if( !MaximaValues.empty() && !MinimaValues.empty() ){ // in specific cases, there might not be extrema??? then don't go into calcs with the vectors!
			// Trans1Start can use the same function for search as 1End, only with another threshold - due to definition of threshold, it is the inverted value, so bigger then the upper ones
			POIMin = FindTransZone1EndPosMin( Zone1StartThreshold_rel, MinimaValues ); // so this threshold should be something like 0.8, meaning when extremum is higher then 0.2
	       		//	relative to MaxCenter
			POIMax = FindTransZone1EndPosMax( Zone1StartThreshold_rel, MaximaValues );
			commonelectrontraj.TransitionZone1StartGC[1] = (MaxGCX[POIMax] + MinGCX[POIMin])/2;
			commonelectrontraj.TransitionZone1StartGC[2] = (MaxGCY[POIMax] + MinGCY[POIMin])/2;
			commonelectrontraj.TransitionZone1StartGC[3] = (MaxGCZ[POIMax] + MinGCZ[POIMin])/2;

			// Trans1End 
			POIMin = FindTransZone1EndPosMin( Zone1EndThreshold_rel, MinimaValues );
			POIMax = FindTransZone1EndPosMax( Zone1EndThreshold_rel, MaximaValues );
			commonelectrontraj.TransitionZone1EndGC[1] = (MaxGCX[POIMax] + MinGCX[POIMin])/2;
			commonelectrontraj.TransitionZone1EndGC[2] = (MaxGCY[POIMax] + MinGCY[POIMin])/2;
			commonelectrontraj.TransitionZone1EndGC[3] = (MaxGCZ[POIMax] + MinGCZ[POIMin])/2;

			// Trans2Start 
			POIMin = FindTransZone2StartPosMin( Zone2StartThreshold_rel, MinimaValues ); 
			POIMax = FindTransZone2StartPosMax( Zone2StartThreshold_rel, MaximaValues );
			commonelectrontraj.TransitionZone2StartGC[1] = (MaxGCX[POIMax] + MinGCX[POIMin])/2;
			commonelectrontraj.TransitionZone2StartGC[2] = (MaxGCY[POIMax] + MinGCY[POIMin])/2;
			commonelectrontraj.TransitionZone2StartGC[3] = (MaxGCZ[POIMax] + MinGCZ[POIMin])/2;

			// Trans2End
			// here we can again reuse the Zone2Start function just with another threshold - again invert value
			POIMin = FindTransZone2StartPosMin( Zone2EndThreshold_rel, MinimaValues ); 
			POIMax = FindTransZone2StartPosMax( Zone2EndThreshold_rel, MaximaValues );
			commonelectrontraj.TransitionZone2EndGC[1] = (MaxGCX[POIMax] + MinGCX[POIMin])/2;
			commonelectrontraj.TransitionZone2EndGC[2] = (MaxGCY[POIMax] + MinGCY[POIMin])/2;
			commonelectrontraj.TransitionZone2EndGC[3] = (MaxGCZ[POIMax] + MinGCZ[POIMin])/2;
		}


	} // end of transition calcs for Cluster MC


		
	//cout << "end z" << x[3] << endl;
	
	
	return;
}





//////////////////////////////////////////////////////////

void velocity(double *n,double kin_energy,double costheta,double phi,double *v){    // calculate particle velocity for given energy, B-field, and direction
/* This subr. computes the particle velocity vector (vx, vy and vz components),
    as a function of its kinetic energy, polar and azimuthal angles theta and
    phi (costheta=cos(theta)), relative to the unit vector n pointing
    from the electrode surface inside to the spectrometer.
Units:
    v[i] (i=1,2,3): meter/s; kin_energy: eV; theta and phi: radian  */
    const double e=1.602177e-19;// electron charge (in SI, without sign) [C]
    const double c=299792458.;  // velocity of light in SI units [m/s]
    double V;                   // magnitude of velocity [m/s]
    double gam,beta;            // usual relativistic parameters
    double gam1,xz,na[4],nb[4],sintheta; // variables used for the calulations
    int i;                      // loop counter
  double Phi,E[4],B[4],b,phirad; // unused variables
    // Magnitude of velocity calculation:
    gam1=kin_energy*e/(MASS*c*c);   // ratio between kinetic and rest energy
    gam=gam1+1.;                    // gamma factor of the charged particle (electron/proton)
    beta=sqrt(gam1*(gam+1.)/(gam*gam)); // calculate beta (v/c) of the particle
    V=beta*c;                       // calculate the magnitute of the particl velocity [m/s]
    // Calc. of unit vectors na, nb (orthogonal to n):
    xz=sqrt(n[1]*n[1]+n[3]*n[3]);   // magnitude of the projection of the B-field's unity vector on the xz plane
    if(xz>1.e-12){na[1]=-n[3]/xz; na[2]=0.; na[3]=n[1]/xz;} // creat a unit vector perpenticular to the project y simply rotate the vector in the xy plane and renormalize
    else{na[1]=1.; na[2]=0.; na[3]=0.;}   // if xy is to small than the B-field points approximatly in the y-direction
    // Calculating nb = n x na (cross/vector product)
    nb[1]=n[2]*na[3]-n[3]*na[2];
    nb[2]=n[3]*na[1]-n[1]*na[3];
    nb[3]=n[1]*na[2]-n[2]*na[1];
    // sintheta = sin(theta):
    //
    //debugging vectors
    //cout << "TRAJBUG: normal vectors:" << endl << n[1] << "\t" << n[2] << "\t" << n[3] << endl;
    //cout << na[1] << "\t" << na[2] << "\t" << na[3] << endl;
    //cout << nb[1] << "\t" << nb[2] << "\t" << nb[3] << endl;
    //
    //
    sintheta=sqrt(fabs(1.-costheta*costheta));
    // Velocity vector computation:
    for(i=1;i<=3;i++)
        v[i]=V*(n[i]*costheta+sintheta*(na[i]*cos(phi)+nb[i]*sin(phi)));
    /*std::cout << setprecision(4) << setw(10) << costheta << setw(16) << v[1]/V << setw(16) << v[2]/V << setw(16) << v[3]/V << endl;
    std::cout << setprecision(4) << setw(12) << n[1]  << setw(12) << na[1]  << setw(12) << nb[1] << endl;
    std::cout << setprecision(4) << setw(12) << n[2]  << setw(12) << na[2]  << setw(12) << nb[2] << endl;
    std::cout << setprecision(4) << setw(12) << n[3]  << setw(12) << na[3]  << setw(12) << nb[3] << endl;*/
    return;         // the particle coordniats x and velocity v are returned
    }

///////////////////////////////////////

void electronstart(double Start_cord[4],double r,double costheta,double phidir,double Ekin,double *x,double *v){
  double PI = 3.1415926;
    double mz,mr,L,nr,vlong,vtrans;
  double phispec,cosphi,sinphi,kin_energy;
    // used variables
    double B[4],b, n[4];
    int k;
    // Spatial coordinates of the Starting point, random phi angle on the DISK!
    //phi=2.*PI*randomnumber();
    x[1]=Start_cord[1];
    x[2]=Start_cord[2];
    x[3]=Start_cord[3];
    // Calculate the B-field in the Starting point
    magfieldtraj(x,B,commonelectrontraj.filepath);
    b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);
	commonelectrontraj.StartB[1] = B[1]; 
	commonelectrontraj.StartB[2] = B[2]; 
	commonelectrontraj.StartB[3] = B[3]; 
    // Calculate the unity vector n of the B-field:
    if(B[3]>0.)for(k=1;k<=3;k++)n[k]=B[k]/b;
    else{for(k=1;k<=3;k++)n[k]=-B[k]/b;}
    // calculate the velocity components:
    velocity(n,Ekin,costheta,phidir,v);
    return;         // the particle coordniats x and velocity v are returned
    }

///////////////////////////////////

void electrongeneration(double *x,double *v, double phidir){       // generating particles at the source
    double PI = 3.1415926;
    double Start_cord[4],r,Ekin,costheta;
    // Position of the starting point:
    Start_cord[1]=commonelectrontraj.Xstart;
    Start_cord[2]=commonelectrontraj.Ystart;
    Start_cord[3]=commonelectrontraj.Zstart;
    // distance on a disk around the starting point (in order to simulat an appertur or beam size)
    r=commonelectrontraj.Rstartmax*sqrt(randomnumber());
    // Direction of the emitted particle:
    if(commonelectrontraj.theta_fix == 0)costheta=1.-(1.-cos(commonelectrontraj.thetaStartmax))*randomnumber();   // theta = angle between z-axis and particle direction
    else costheta= cos(commonelectrontraj.thetaStartmax);
    // define particle energy
    Ekin=commonelectrontraj.Ekin;
    // calculate the actual start position and velocity vector
    electronstart(Start_cord,r,costheta,phidir,Ekin,x,v);
    return;         // the particle coordniats x and velocity v are returned
    }


/////////////////////////////////////////////

void subrn(double *u,int len)       // subroutinefor calculating random numbers
{
// This subroutine computes random numbers u[1],...,u[len]
// in the (0,1) interval. It uses the 0<IJKLRANDOM<900000000
// integer as initialization seed.
//  In the calling program the dimension of the u[] vector
// should be larger than len (the u[0] value is not used).
// For each IJKLRANDOM
// numbers the program computes completely independent random number
// sequences (see: F. James, Comp. Phys. Comm. 60 (1990) 329, sec. 3.3).
  static int iff=0;
  static long ijkl,ij,kl,i,j,k,l,ii,jj,m,i97,j97,ivec;
  static float s,t,uu[98],c,cd,cm,uni;
  if(iff==0)
  {
    if(IJKLRANDOM==0)
    {
      printf("Message from subroutine subrn:\
                 the global integer IJKLRANDOM should be larger than 0 !!!\
                 Computation is  stopped !!! \n");
      exit(0);
    }
    ijkl=IJKLRANDOM;
    if(ijkl<1 || ijkl>=900000000) ijkl=1;
    ij=ijkl/30082;
    kl=ijkl-30082*ij;
    i=((ij/177)%177)+2;
    j=(ij%177)+2;
    k=((kl/169)%178)+1;
    l=kl%169;
    for(ii=1;ii<=97;ii++)
    { s=0; t=0.5;
      for(jj=1;jj<=24;jj++)
      { m=(((i*j)%179)*k)%179;
        i=j; j=k; k=m;
        l=(53*l+1)%169;
        if((l*m)%64 >= 32) s=s+t;
        t=0.5*t;
      }
      uu[ii]=s;
    }
    c=362436./16777216.;
    cd=7654321./16777216.;
    cm=16777213./16777216.;
    i97=97;
    j97=33;
    iff=1;
  }
  for(ivec=1;ivec<=len;ivec++)
  { uni=uu[i97]-uu[j97];
    if(uni<0.) uni=uni+1.;
    uu[i97]=uni;
    i97=i97-1;
    if(i97==0) i97=97;
    j97=j97-1;
    if(j97==0) j97=97;
    c=c-cd;
    if(c<0.) c=c+cm;
    uni=uni-c;
    if(uni<0.) uni=uni+1.;
    if(uni==0.)
    { uni=uu[j97]*0.59604644775391e-07;
      if(uni==0.) uni=0.35527136788005e-14;
    }
    u[ivec]=uni;
  }
  return;
}

////////////////////////////////////////////////////////////////

double randomnumber(){      // random number generater (0,1)
    // This function computes 1 random number in the (0,1) interval,
    // using the subrn subroutine.
    double u[2];
    subrn(u,1);         // calculat random number
    return u[1];        // return random number
    }

///////////////////////////////////////////////////////////
//
//
//
//
//function for calculating the local phase of a particle dependent on: v, B, q,
//is currently under the assumption that x-unit vector stays perpendicular to B (roughly)
//

double PhaseAngleFunc(double *v, double *B, int charge){

	double ForceVector[4], LocalYVector[4];
	double LocalYNorm;
	double px, py;

	// v x B
	ForceVector[1] = (v[2]*B[3] - v[3]*B[2])*(double)charge;
	ForceVector[2] = (v[3]*B[1] - v[1]*B[3])*(double)charge;
	ForceVector[3] = (v[1]*B[2] - v[2]*B[1])*(double)charge;
	
	// calc the local y axis vector and its norm. as its with cross product with x, it is easy to calc
	LocalYNorm = sqrt(pow(B[2],2) + pow(B[3],2));
	LocalYVector[2] = B[3]/LocalYNorm;
	LocalYVector[3] = -B[2]/LocalYNorm;

	//now calc the projections
	px = ForceVector[1];
	py = ForceVector[2] * LocalYVector[2] + ForceVector[3] * LocalYVector[3];
	
	// now return phi
	return atan2(py,px);

}

/////////////////////////////////////////////////////////////////
//
//
//
//function for adiabatic theta local
double ThetaLocal(double theta0, double B0, double Blocal){

	return asin(sin(theta0)*sqrt(Blocal/B0));

}



//////////////////////////////////////////////////////////
//
//
//function to calculate the length between 2 3D-points
double PointDiff(double *x1, double *x2){
	return sqrt( pow( x1[1] - x2[1] , 2) + pow( x1[2] - x2[2] , 2) + pow( x1[3] - x2[3] , 2) );
}






////////////////////////////////////////////////////////////
//
//
//
// GC calculation

//void GuidingCenterFunc(double *x, double *v, double *B, double theta, double gamma, double b_norm, double v_norm, int charge){
//
//	double ForceVector[4];
//	double ForceNorm;
//	double rG;
//	double MASS=9.1093836e-31;      // electron mass (in SI) [kg]
//	
//	// v x B
//	ForceVector[1] = v[2]*B[3] - v[3]*B[2];
//	ForceVector[2] = v[3]*B[1] - v[1]*B[3];
//	ForceVector[3] = v[1]*B[2] - v[2]*B[1];
//	//Norm of force vector (without charge still)
//	ForceNorm = sqrt(ForceVector[1]*ForceVector[1] + ForceVector[2]*ForceVector[2] + ForceVector[3]*ForceVector[3]); 
//	v_norm = sqrt(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
//	b_norm = sqrt(B[1]*B[1] + B[2]*B[2] + B[3]*B[3]);
//	v_para = (v[1]*B[1] + v[2]*B[2] + v[3]*B[3])
//        gam=1./sqrt(1.-pow(v_norm/c,2.));     // gamma factor:
//	rG = gam * MASS * V_perpend/ 1.602177e-19 / b;
//
//}


void GuidingCenterFunc(double *x, double *v, double *B, double costheta, double gamma, double b_norm, double v_norm, double *GC){
	
	double sign;
	double V_perpend;
	double CrossP[4];
	double gyraR;
	double CrossPL;

	V_perpend = sin(acos(costheta)) * v_norm;
	if(CHARGE < 0) sign = -1.;
	else sign = 1.;
	gyraR = gamma * MASS * V_perpend/ 1.602177e-19 / b_norm;
	// v x B
	CrossP[1] = B[3]*v[2] - B[2]*v[3];
	CrossP[2] = -B[3]*v[1] + B[1]*v[3];
	CrossP[3] = B[2]*v[1] - B[1]*v[2];
	CrossPL = sqrt(CrossP[1]*CrossP[1] + CrossP[2]*CrossP[2] + CrossP[3]*CrossP[3]);
	if( CrossPL < 1e-12 ){
		for(int i=1;i<=3;i++){GC[i] = x[i];}
	}
	else{
		for(int i=1;i<=3;i++){GC[i] = x[i] + sign* CrossP[i]/CrossPL * gyraR;};
	}

}



////////////////////////////////
// search max of a vector and return position

int VectorFindMax(vector<double> &MaxList){

	double MaxVal = 0.;
	int MaxIt = 0;
	for( int i = 0; i < MaxList.size(); i++){
		if(MaxList[i] > MaxVal){
		       MaxIt = i;
		       MaxVal = MaxList[i];
		}
	}

	return MaxIt;
}


////////////////////////////////
// search max of a vector and return position
int VectorFindMin(vector<double> &MinList){

	double MinVal = 0.;
	int MinIt = 0;
	for( int i = 0; i < MinList.size(); i++){
		if(MinList[i] < MinVal){
		       MinIt = i;
		       MinVal = MinList[i];
		}
	}

	return MinIt;
}


///////////////////////////////////////
// give position of TransitionZone1End by giving relative threshold relative to the plateau given my the max value of Extrema
// we search both in the minima list as well as maxima list and combine the results by averaging the Guiding center positions of the found points
// we separate in two functions for min and max so that we can just return one integer each
int FindTransZone1EndPosMax( double Threshold, vector<double> &MaxList ){

	double MaxCenter = MaxList[ MaxList.size()/2 ]; // gives center Maximum - rounded down always (int)
	int i = 0;
	for( i = 0; i < MaxList.size(); i++ ){
		if( (MaxCenter - MaxList[i])/MaxCenter < Threshold ) break;
	}

	return i;
}

int FindTransZone1EndPosMin( double Threshold, vector<double> &MinList ){

	double MinCenter = MinList[ MinList.size()/2 ]; // gives center Maximum - rounded down always (int)
	int i = 0;
	for( i = 0; i < MinList.size(); i++ ){
		if( (MinCenter - MinList[i])/MinCenter < Threshold ) break;
	}

	return i;
}


/////////////////////////////////////////
//for start of Zone 2, we consider relative to the maximum value of maxima etc.
int FindTransZone2StartPosMax( double Threshold, vector<double> &MaxList ){

	int MaxMaxPos = VectorFindMax( MaxList );
	double MaxMax = MaxList[MaxMaxPos]; // gives Maximum of Maxima 
	int i = 0;
	for( i = MaxMaxPos; i < MaxList.size(); i++ ){
		if( (MaxMax - MaxList[i])/MaxMax > Threshold ) break;
	}

	return i;
}

int FindTransZone2StartPosMin( double Threshold, vector<double> &MinList ){

	int MinMinPos = VectorFindMin( MinList );
	double MinMin = MinList[MinMinPos]; // gives Maximum of Maxima 
	int i = 0;
	for( i = MinMinPos; i < MinList.size(); i++ ){
		if( (MinMin - MinList[i])/MinMin > Threshold ) break;
	}

	return i;
}












