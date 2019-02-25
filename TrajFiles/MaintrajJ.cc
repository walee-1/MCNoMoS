#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>                      // for stringstream
#include <random>

using namespace std;

#include "../ConfigFiles/config.h"

#define Ngroup 1
#define Nmax 3500
#define nmax 1000
#define Nspmax 2000
#define nmaxmag 4000 // 500
#define Nspmagmax 3000
#define _USE_MATH_DEFINES


#include "Traj.cc"


////////////////////////////////////////////////////////

int main(int argc, char ** argv, char* envp[])
{
    
    double EfromP(double);
    
	int N,iEkin, iRandom, psteps;
	string conclusionfilename;
	char numstr1[21],numstr2[21];
	double pmin, pmax, pstep;
	int horilines, vertilines;
	double ApertX, ApertY;
	ofstream ConclusionHeader;
	double StartZ;

////////////////////// Input tests //////////////////////

	if(argc<2){cout << "TRAJ: No PATH given!" << endl; return 1;}
	string filedir=argv[1];
	cout << endl <<  "TRAJ: Given PATH is: " << filedir << endl;

	
	/////////////Config parser//////////////
	
	if(argc <3) {cout << "TRAJ: no config file given!" << endl; return 1;}
	string configname = argv[2];
	cout << "TRAJ: config given: " << configname << endl;
	Config myconfig(configname,envp);// dir, where executable will be located, so without "../"

	//  Simulation parameters:
	commonelectrontraj.filepath = filedir; 
	commonelectrontraj.PercOn =myconfig.pBool("PercOn");
	commonelectrontraj.NoMoSOn = myconfig.pBool("NoMoSOn");
	commonelectrontraj.CompareWBlines = myconfig.pBool("CompareWBlines");
	commonelectrontraj.MonteCarlo = myconfig.pBool("MonteCarloData");
	commonelectrontraj.G4Compare = myconfig.pBool("G4Compare");
	string particle = myconfig.pString("particle");
	N = myconfig.pInt("ParticleN");
	double apertYshift = myconfig.pDouble("apertYshift");
	double apertXshift = myconfig.pDouble("apertXshift");
	bool onlyCornerCenter = myconfig.pBool("onlyCornerCenter");
	horilines = myconfig.pInt("horilines");
	vertilines = myconfig.pInt("vertilines");
	ApertX = myconfig.pDouble("ApertX");
	ApertY = myconfig.pDouble("ApertY");
	double detpos = myconfig.pDouble("detpos");
	
	int lines;
	if(onlyCornerCenter == true) lines =5;
	else lines = horilines*vertilines;

	double startP[2][lines];
	ofstream OutStream[lines];
    	stringstream OutStringStream;
	int counter = 0;

	// here we set starting points. differentiation if onlyCornerCenter == True
	if( onlyCornerCenter == true){
		//first the center point
		startP[0][counter] = 0.; startP[1][counter] = 0.;
		counter++;
		
		for(int x=-1; x <2; x=x+2){
			for(int y=-1; y<2; y=y+2){
				startP[0][counter]=x*ApertX/2.;
				startP[1][counter]=y*ApertY/2.;
				counter ++;
			}
		}
	}
	else{
    		double xstep;
		if( horilines == 1 ){
			xstep = 0.;
		}
		else{
			xstep = ApertX/(double)(horilines-1);
		}
    		double ystep;
        	if( vertilines == 1 ){
			ystep = 0.;
		}
		else{ ystep = ApertY/(double)(vertilines-1);}

    		// define start points dependent on apertsize and hori/verti-lines. NOT REAL COORDS. because we use these for the file names!
    		for(int x=1;x<=horilines; x++){
    		    for(int y=1;y<=vertilines; y++){
    		        startP[0][counter]= -ApertX/2. + (x-1)*xstep;
    		        startP[1][counter]= -ApertY/2. + (y-1)*ystep;
    		        counter++;
    		    }
    		}
	}


	// read in info file from sim so R1 is read in!
	string infoname = filedir + "info.txt";
	ifstream infofile;
	infofile.open(infoname.c_str(),ios::in);
	double R_1;
	infofile >> infoname >> R_1;
	cout << "TRAJ: Read-in R_1 = " << R_1 << endl;
	infofile >> infoname >> StartZ;
	infofile.close();
	commonelectrontraj.R_1=R_1;
	
	// seed of random number generator :
	iRandom = 10;//atoi(argv[2]);
	IJKLRANDOM=iRandom;
	
	// Particle properties:
	// Electron charge and mass:
	if ( particle == "e" ){
		CHARGE=-1.6021766e-19;   // electron charge (in SI, with sign) [C]
		MASS=9.1093836e-31;      // electron mass (in SI) [kg]
	}
	else{
		CHARGE=1.6021766e-19;   // electron charge (in SI, with sign) [C]
		MASS=1.672621e-27;      // electron mass (in SI) [kg]
	}


	// Start conditions:
	// Starting position (x,y,z) in [m]:
	commonelectrontraj.Zstart= StartZ;
	cout << "TRAJ: Zstart = " << commonelectrontraj.Zstart << endl;
	
	// Detector conditions
    	commonelectrontraj.Zend=R_1 + 0.5;        // Final z position
    	commonelectrontraj.zdetector = detpos;    // detector position

	// Simulation settings
	commontrajexact.typeh=2;            // automatic timestep computation
	commontrajexact.typerungekutta=8;   // typerungekutta=4: rungekutta4;
	commontrajexact.slimit=5.e-4;       // pathlength limit for 1 RK-step (m)
	commontrajexact.ntimestep=40;       // RK-step control by cyclotron period
	commonelectrontraj.numstepmax=100000000;  // maximal number of Runge-Kutta steps
	commonelectrontraj.Timemax=30.e-3;  // Maximal time [s]
	commonelectrontraj.typeprint=1;     // Printing of intermediate quantities  : 1;  no printing:  0
	facB=1.;        // scaling factor for magnetic field (1 = turned on)
	facE=0.;        // scaling factor for electric field (0 = turned off)
	
	//load earth magnetic field vector from config file, variables are initialized in fieldnew.cpp, and also used there in the function magfieldtraj()
	// numbers in config in micro tesla!!!
	earth_Bx = myconfig.pDouble("earth_BN");
	earth_By = myconfig.pDouble("earth_BE");
	earth_Bz = myconfig.pDouble("earth_Bvert");
	
	// we turn around the coordinates of the earth magnetic field dependent on ILL OR PERC
	if( !commonelectrontraj.PercOn ){ //ILL
		earth_Bx = -earth_Bx/sqrt(2.) - earth_By/sqrt(2.);
		earth_By = -earth_Bx/sqrt(2.) + earth_By/sqrt(2.);
	}
	else{
		earth_Bx = -earth_Bx/sqrt(2.) + earth_By/sqrt(2.);
		earth_By = -earth_Bz;
		earth_Bz = earth_Bx/sqrt(2.) + earth_By/sqrt(2.);
	}
    

////////////// trajectory calcs	///////////////////////////////////
	
    
//////////////////////////////////////////
	////// COMPARE traj with Blines /////////
    //////////////////////////////////////
    
	if( commonelectrontraj.CompareWBlines ){
		cout << "TRAJ: Compare with BLines" << endl;
       
       
       		// conclusionfilename for simulation conclusion
       		//sprintf(numstr1, "%d", iEkin);
       		sprintf(numstr2, "%d", iRandom);
       		conclusionfilename="Conclusion_";
       		//conclusionfilename+=numstr1;
       		conclusionfilename+="iRand";
       		conclusionfilename+=numstr2 ;
       		conclusionfilename+=".txt";
       		conclusionfilename = commonelectrontraj.filepath +"trajectories/"+ conclusionfilename;
       		cout << "TRAJ: Conclusion filename: " << conclusionfilename << endl;
       		ConclusionHeader.open(conclusionfilename.c_str(),ios::out);
       		ConclusionHeader << "X[m]" << "\t" << "Y[m]" << "\t" << "theta[rad]" <<"\t"<< "Ekin[eV]" <<"\t"<< "phi[rad]" << "\t" << "particle" << "\t" << "part-index" << "\t" << "E-err" << "\t" << "time" << endl;
       		ConclusionHeader.close();
       
       
       
		commonelectrontraj.theta_fix = 1; //we set the angle in a loop
		// Radius of starting disk in [m]:
		commonelectrontraj.Rstartmax=0.000001;
	
		pmin = myconfig.pDouble("PMIN");
		pmax = myconfig.pDouble("PMAX");
		pstep = myconfig.pDouble("PSTEP");
		psteps = (pmax - pmin)/pstep;

		double thetastep = myconfig.pDouble("thetastep")/180.*M_PI;
		double thetamax = myconfig.pDouble("thetamax")/180.*M_PI;
		int thetasteps = thetamax / thetastep;
		
		double phistep = myconfig.pDouble("phistep")/180.*M_PI;
		double phimax = myconfig.pDouble("phimax")/180.*M_PI;
		int phisteps = phimax / phistep;


		double velo[4], B[4], n[4], GuidStart[4];
		double b, gam;
		double StartgyraR, CrossPL;
		long double v, V_perpend;
    		double c=299792458.;  // velocity of light in SI units
		double n_length;

		
		//loop for different start points
		for( int line = 0; line < lines; line++){
			cout << endl << "TRAJ: line: " << line+1 << endl;
	
			commonelectrontraj.ApertGridX=startP[0][line];
			commonelectrontraj.ApertGridY=startP[1][line];
			commonelectrontraj.ApertGridZ=commonelectrontraj.Zstart;

			//angle variation
			for(int angle_i=0; angle_i<= thetasteps; angle_i++){
				commonelectrontraj.thetaStartmax = angle_i*thetastep;
				cout << endl << "TRAJ: theta= " << commonelectrontraj.thetaStartmax/M_PI*180. << endl;
				// energy variation
				for(int P_i=0; P_i<=psteps; P_i++){
					// again set X and Y start, so that in E loop, Y is again correct after applying gyro in previous step
					commonelectrontraj.Ekin = EfromP( pmin*1000. + P_i*pstep*1000. );
					cout << endl << "TRAJ: P= " << pmin*1000. + P_i*pstep*1000. << " (E=" << commonelectrontraj.Ekin << ")" << endl;

					for(int phi_i=0; phi_i<=phisteps; phi_i++){

						commonelectrontraj.Ystart= startP[1][line] + commonelectrontraj.R_1 + apertYshift; // we add R_1 in Traj.cc aswellas gyraR
						commonelectrontraj.Xstart= startP[0][line] + apertXshift; 

						commonelectrontraj.Startphi = phi_i*phistep;
						cout << endl << "TRAJ: Startphi= " << commonelectrontraj.Startphi/M_PI*180. << endl;

						//for gyrationR, we need perpend. velocity. B at Guiding Center:
						GuidStart[1] = commonelectrontraj.Xstart;
						GuidStart[2] = commonelectrontraj.Ystart;
						GuidStart[3] = commonelectrontraj.Zstart;
						magfieldtraj(GuidStart,B,filedir);
 				                b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);
						//calc unit vector to define direction for the velocity calc, in B direction
						n[1] = B[1]/b;
						n[2] = B[2]/b;
						n[3] = B[3]/b;
						velocity(n,commonelectrontraj.Ekin,cos(commonelectrontraj.thetaStartmax),commonelectrontraj.Startphi,velo); //phi varied in 90° steps
						
						// calculations for gyrationR and finally point, where to start for Guiding center at Bline-position
 				                v=sqrt(velo[1]*velo[1]+velo[2]*velo[2]+velo[3]*velo[3]);
						V_perpend =sin(commonelectrontraj.thetaStartmax)*v;
                  				//cout << "beta = v/c = " << v/c << endl;
                  				gam=1./sqrt(1.-pow(v/c,2.));     // gamma factor:
						StartgyraR = gam * MASS * V_perpend/ 1.602177e-19 / b;
						cout << "StartgyraR= " << StartgyraR << endl;
						
						if( V_perpend > 1.e-12 ){
							//now calc cross product q * v x B for direction of offset of starting point
							n[1] = CHARGE*(velo[2]*B[3] - velo[3]*B[2]);
							n[2] = CHARGE*(velo[3]*B[1] - velo[1]*B[3]);
							n[3] = CHARGE*(velo[1]*B[2] - velo[2]*B[1]);
							n_length = sqrt(n[1]*n[1]+n[2]*n[2]+n[3]*n[3]);
							n[1] /= n_length;
							n[2] /= n_length;
							n[3] /= n_length;

							commonelectrontraj.Xstart -= n[1]*StartgyraR;
							commonelectrontraj.Ystart -= n[2]*StartgyraR;
							commonelectrontraj.Zstart -= n[3]*StartgyraR;
						}


						//last argument is offset of Ystart, so that in filename, gyraR can again be subtracted!	
						trajelectronN(conclusionfilename,N,StartgyraR,OutStream[line]); 
					}
				}
			}
		}
	}
    
    
    
    
    //////////////////////////////////////////////////
    // Monte Carlo Trajectories with only end positions
    //////////////////////////////////////////////
    
    if( commonelectrontraj.MonteCarlo == true ){
    
        cout << endl << "TRAJ: MonteCarlo Data" << endl << endl;
        counter = 0;
	commonelectrontraj.hemisphere = 1;
	OutStringStream.str("");
	ofstream MonteCarloOut;
        double buffer1, buffer2, buffer3;
	stringstream filenr;
	//we save the aperture size to one existing global value to have it in Traj.cc
	commonelectrontraj.ApertGridX = ApertX;
	commonelectrontraj.ApertGridY = ApertY;

	double ApertExpander = 0.016*4; // how much the aperture is assumed larger for starting positions
	
	//random number initialization
	random_device r;
	seed_seq seed{r(),r(),r(),r(),r(),r(),r(),r()};
	mt19937 eng{seed};
	uniform_real_distribution<double> dist(-1.,1.);

	//
	//

	commonelectrontraj.theta_fix = 1;
        commonelectrontraj.Rstartmax=0.000001;
        commonelectrontraj.typeprint=0;
	

	//here we loop over the single partial files of the generated spectrum with increasing index
        // read-in from MC
        ifstream MonteCarloData;
        string SpectraDir = "/home/dmoser/FerencSource+Spectra/Spectra1e7/";

        string SpectraNr;
      
	// now we loop over the single files starting with 0
	// or 1 in the mean while (first 10^5 is already done. now in single files!
	for(int nr=2;nr<100;nr++){

		filenr.str(""); // clearing the stringstream before each iteration
		filenr << nr;
		SpectraNr= "Spectra1e5_" + filenr.str() + "_filtered.txt";
	        MonteCarloData.open( (SpectraDir+SpectraNr).c_str() ,ios::in);
		cout << "MonteCarloBUG: Opened File = " << (SpectraDir+SpectraNr).c_str() << endl;
	
		// now make different output file for each 10^5 package
		OutStringStream.str(""); // clearing
		OutStringStream << commonelectrontraj.filepath << "MonteCarlo/"  << "DetectorData"+ filenr.str() +".txt";
		conclusionfilename = OutStringStream.str(); // conclusionfilename is just reused here for opening the files
		MonteCarloOut.open(conclusionfilename.c_str(),ios::app);
		// Header
		MonteCarloOut << "XStart" << "\t" << "YStart" << "\t" << "ZStart" << "\t" << "x" << "\t" << "y" << "\t" << "z" << "\t" << "vx" << "\t" << "vy" << "\t" << "vz" << "\t" << "Ekin" << "\t" << "hemi" <<"\t" << "apertX" << "\t" << "apertY" <<  endl;
	

	        while ( MonteCarloData.good() ){
	    	    cout << endl << "TRAJ: Decay#: " << counter << endl;	    
	    	    counter ++;
	
		    if ( particle == "e" ){ // if electron, take second part of data, else the first part for proton
	            	MonteCarloData >> buffer1 >> buffer2 >> buffer3 >> commonelectrontraj.thetaStartmax >> commonelectrontraj.Startphi >> commonelectrontraj.Ekin;
		    }
		    else{
	            	MonteCarloData >> commonelectrontraj.thetaStartmax >> commonelectrontraj.Startphi >> commonelectrontraj.Ekin >> buffer1 >> buffer2 >> buffer3 ;
		    }
	
		    commonelectrontraj.hemisphere = 1;
		    if( commonelectrontraj.thetaStartmax >= 90. ){
			    commonelectrontraj.thetaStartmax = 180. - commonelectrontraj.thetaStartmax;
			    commonelectrontraj.hemisphere = -1;
		    }
	
		    // if theta is greater than 45°, filter will reflect -> cutoff
		    if( commonelectrontraj.thetaStartmax < 45. ){
		
			commonelectrontraj.Startphi = commonelectrontraj.Startphi/180.*M_PI;
			commonelectrontraj.thetaStartmax = commonelectrontraj.thetaStartmax /180. * M_PI;
			commonelectrontraj.Ekin *=1000.; // from keV to eV!
		   	

			// now we still have to random variate the starting position, we choose a bigger window than the aperture to see the full edge effect
				

			commonelectrontraj.Ystart= commonelectrontraj.R_1 + apertYshift + (ApertY+ApertExpander)/2. * dist(eng); //dist(eng) gives numb between -1 and 1 for selection of aperture position
			commonelectrontraj.Xstart= apertXshift+ (ApertX+ApertExpander)/2. * dist(eng); 
			//Zstart left constant as defined in beginning of main, for now

	           	trajelectronN(conclusionfilename,1,0., MonteCarloOut ); //filenam, N, startgyraR, apertYshift,apertXshift
	
		    }
		    else cout << "MonteCarloBug: theta > 45" << endl;
	            

	        } //closing the while loop
	
		MonteCarloOut.close();
		MonteCarloData.close();

	}// for loop closing of single monte carlo files

		//dont forget to close the file streams
		//
		//

    }//closing monte carlo simulations




/////////////////////////////////
//we comment out the following code, which is the old first version of Monte Carlo tracking implemented (starting points were still on a grid)

/*
	// define Outstream for all the lines
        for( int line = 0; line < lines; line++){
		
		OutStringStream.str("");
		conclusionfilename = "m";
		OutStringStream << commonelectrontraj.filepath << "MonteCarlo/" << "X" << startP[0][line] << "_Y" << startP[1][line]  << ".txt";
		conclusionfilename = OutStringStream.str(); // conclusionfilename is just reused here for opening the files
		OutStream[line].open(conclusionfilename.c_str(),ios::out);
		// Header
		OutStream[line] << "x" << "\t" << "y" << "\t" << "z" << "\t" << "vx" << "\t" << "vy" << "\t" << "vz" << "\t" << "Ekin" << "\t" << "hemi" << endl;
		
        }	

		
	commonelectrontraj.theta_fix = 1;
        commonelectrontraj.Rstartmax=0.000001;
        commonelectrontraj.typeprint=0;
        
        // read-in from MC
        ifstream MonteCarloData;
        string SpectraDir = "/home/dmoser/Dropbox/PhD/Magfield/FerencSource/Spectra1e7/";
        string SpectraNr = "Spectra1e5_0.txt";
        MonteCarloData.open( (SpectraDir+SpectraNr).c_str() ,ios::in);
 
        double buffer1, buffer2, buffer3;
        
        while ( MonteCarloData.good() ){
    	    cout << endl << "TRAJ: Decay#: " << counter << endl;	    
    	    counter ++;

            // for each monte carlo line, read in data for decay
	    if ( particle == "e" ){ // if electron, take second part of data, else the first part for proton
            	MonteCarloData >> buffer1 >> buffer2 >> buffer3 >> commonelectrontraj.thetaStartmax >> commonelectrontraj.Startphi >> commonelectrontraj.Ekin;
	    }
	    else{
            	MonteCarloData >> commonelectrontraj.thetaStartmax >> commonelectrontraj.Startphi >> commonelectrontraj.Ekin >> buffer1 >> buffer2 >> buffer3 ;
	    }

	    hemisphere = 1;
	    if( commonelectrontraj.thetaStartmax > 90. ){
		    commonelectrontraj.thetaStartmax = 180. - commonelectrontraj.thetaStartmax;
		    hemisphere = -1;
	    }

	    cout << "TRAJ: theta = " << commonelectrontraj.thetaStartmax << endl;

	    // if theta is greater than 45°, filter will reflect -> cutoff
	    if( commonelectrontraj.thetaStartmax < 40. ){
	
		 commonelectrontraj.Startphi = commonelectrontraj.Startphi/180.*M_PI;
		 commonelectrontraj.thetaStartmax = commonelectrontraj.thetaStartmax /180. * M_PI;
	   	 
		 // line loop, so for each line, each decay will be tracked
           	 for( int line = 0; line < lines; line++){
           	     	cout << endl << "TRAJ: line: " << line+1 << endl;
			
			commonelectrontraj.Ystart= startP[1][line] + commonelectrontraj.R_1 + apertYshift; // we add R_1 in Traj.cc aswellas gyraR
			commonelectrontraj.Xstart= startP[0][line] + apertXshift; 

           		trajelectronN(conclusionfilename,1,0., OutStream[line]); //filenam, N, startgyraR, apertYshift,apertXshift

	   		OutStream[line] << "\t" << hemisphere << endl;
           	                 
           	 }
	   	 // break condition for tests, so that only 10 lines of MC are running
	   	 //if( counter == 10 ) break;
	    }
	    else{ cout << "TRAJ: MC_theta > 40 deg" << endl;}
            
        }
        
        
        
        
        MonteCarloData.close();
        
    }
*/
    


    
	//////////////////////////////////////////
	////// COMPARE with G4 /////////
    	//////////////////////////////////////
    
	if( commonelectrontraj.G4Compare ){
		cout << "TRAJ: Compare with G4" << endl;
       
       
       		// conclusionfilename for simulation conclusion
       		//sprintf(numstr1, "%d", iEkin);
       		sprintf(numstr2, "%d", iRandom);
       		conclusionfilename="ConclusionG4_";
       		//conclusionfilename+=numstr1;
       		conclusionfilename+="iRand";
       		conclusionfilename+=numstr2 ;
       		conclusionfilename+=".txt";
       		conclusionfilename = commonelectrontraj.filepath +"trajectories/"+ conclusionfilename;
       		cout << "TRAJ: Conclusion filename: " << conclusionfilename << endl;
       		ConclusionHeader.open(conclusionfilename.c_str(),ios::out);
       		ConclusionHeader << "X[m]" << "\t" << "Y[m]" << "\t" << "theta[rad]" <<"\t"<< "Ekin[eV]" <<"\t"<< "phi[rad]" << "\t" << "particle" << "\t" << "part-index" << "\t" << "E-err" << "\t" << "time" << endl;
       		ConclusionHeader.close();
       
       
       
		commonelectrontraj.theta_fix = 1; //we set the angle in a loop
		// Radius of starting disk in [m]:
		commonelectrontraj.Rstartmax=0.000001;
	
		pmin = myconfig.pDouble("PMIN");
		pmax = myconfig.pDouble("PMAX");
		pstep = myconfig.pDouble("PSTEP");
		psteps = (pmax - pmin)*pstep;

		double thetastep = myconfig.pDouble("thetastep")/180.*M_PI;
		double thetamax = myconfig.pDouble("thetamax")/180.*M_PI;
		int thetasteps = thetamax / thetastep;
		
		double phistep = myconfig.pDouble("phistep")/180.*M_PI;
		double phimax = myconfig.pDouble("phimax")/180.*M_PI;
		int phisteps = phimax / phistep;

		
		//loop for different start points
		for( int line = 0; line < lines; line++){
			cout << endl << "TRAJ: line: " << line+1 << endl;
			
			commonelectrontraj.ApertGridX=startP[0][line];
			commonelectrontraj.ApertGridY=startP[1][line];
			commonelectrontraj.ApertGridZ=commonelectrontraj.Zstart;
			
			//angle variation
			for(int angle_i=0; angle_i<= thetasteps; angle_i++){
				commonelectrontraj.thetaStartmax = angle_i*thetastep;
				cout << endl << "TRAJ: theta= " << commonelectrontraj.thetaStartmax/M_PI*180. << endl;
				// energy variation
				for(int P_i=0; P_i<=psteps; P_i++){
					// again set X and Y start, so that in E loop, Y is again correct after applying gyro in previous step
					commonelectrontraj.Ekin = EfromP( pmin*1000. + P_i*pstep*1000. );
					cout << endl << "TRAJ: P= " << pmin*1000. + P_i*pstep*1000. << " (E=" << commonelectrontraj.Ekin << ")" << endl;

					for(int phi_i=0; phi_i<=phisteps; phi_i++){
			
						commonelectrontraj.Ystart= startP[1][line] + commonelectrontraj.R_1 + apertYshift; // we add R_1 in Traj.cc aswellas gyraR
						commonelectrontraj.Xstart= startP[0][line] + apertXshift; 

						commonelectrontraj.Startphi = phi_i*phistep;
						cout << endl << "TRAJ: Startphi= " << commonelectrontraj.Startphi/M_PI*180. << endl;

						
						//last argument is offset of Ystart, so that in filename, gyraR can again be subtracted!	
						trajelectronN(conclusionfilename,N,0.,OutStream[line]); 
					}
				}
			}
		}
	}

//
//	if( !(commonelectrontraj.CompareWBlines) ){
//
//		// maximal emission angle of particles [rad]
//		commonelectrontraj.theta_fix = myconfig.pBool("thetafix");
//		commonelectrontraj.thetaStartmax=myconfig.pDouble("thetamax")/180.*M_PI; 
//
//		commonelectrontraj.Xstart=0;
//		commonelectrontraj.Ystart=commonelectrontraj.R_1;//*(2.-sin(M_PI/4.));
//		// center start, different energies, random angle
//		for (int i = 1; i <= esteps; i++) {                    
//			iEkin = i;
//			commonelectrontraj.iEkin=iEkin;
//			commonelectrontraj.Ekin=iEkin*1000*estep;
//			cout << "TRAJ: E-Loop Ekin = " << commonelectrontraj.Ekin/1000. << "keV" << endl;
//			trajelectronN(filename,N,0.,0.);
//		}
//	}
//



	/*
	// medial upper line
	 commonelectrontraj.Xstart=0.005;
	 commonelectrontraj.Ystart=commonelectrontraj.R_1 - 0.005;
	for (int i = 10; i < 19; i++) {
	  iEkin = i+1;
	  commonelectrontraj.iEkin=iEkin;
	  commonelectrontraj.Ekin=(iEkin-10.)*100000.;
	  if(1==1)trajelectronN(filename,N);
	  }
	
	// lateral upper line
	 commonelectrontraj.Ystart=commonelectrontraj.R_1 + 0.005;
	for (int i = 20; i < 29; i++) {
	  iEkin = i+1;
	  commonelectrontraj.iEkin=iEkin;
	  commonelectrontraj.Ekin=(iEkin-20.)*100000.;
	  if(1==1)trajelectronN(filename,N);
	  }
	
	// central upper line
	 commonelectrontraj.Ystart=commonelectrontraj.R_1;
	for (int i = 30; i < 39; i++) {
	  iEkin = i+1;
	  commonelectrontraj.iEkin=iEkin;
	  commonelectrontraj.Ekin=(iEkin-30.)*100000.;
	  if(1==1)trajelectronN(filename,N);
	  }
	
	// fix theta at different energies
	  commonelectrontraj.thetaStartmax=15./180.*M_PI; //2.563
	  commonelectrontraj.theta_fix = 1;
	   commonelectrontraj.Xstart=0.;
	 commonelectrontraj.Ystart=commonelectrontraj.R_1;
	for (int i = 40; i < 49; i++) {
	  iEkin = i+1;
	  commonelectrontraj.iEkin=iEkin;
	  commonelectrontraj.Ekin=(iEkin-40.)*100000.;
	  if(1==1)trajelectronN(filename,N);
	  }
	
	// fix theta at different energies
	  commonelectrontraj.thetaStartmax=6./180.*M_PI; //2.563
	  commonelectrontraj.theta_fix = 1;
	
	// streight shot
	   commonelectrontraj.theta_fix = 0;
	    commonelectrontraj.thetaStartmax=0.0001/180.*M_PI; //2.563
	for(int i = 74; i < 75; i++){
	    iEkin = i+1;
	    commonelectrontraj.iEkin=iEkin;
	    commonelectrontraj.Ekin=iEkin*10000.;
	    if(1==0)trajelectronN(filename,N);
	}
	*/
	
	  return 0;

}

double EfromP(double p){
    
    return sqrt(pow(p,2) + pow(510998.94,2)) - 510998.94;
}





