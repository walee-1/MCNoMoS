void geoline_anadrift(double R_1,double alpha, int horilines, int vertilines, double apertsizeX, double apertsizeY, double apertYshift, string filedir, string inputcoil, int prep, double momentum, bool driftOn){

	// first write out file for analysis input
	ofstream AnaInput;
	string filename1, filename2;
	if(driftOn){filename1 = "geo_DOn";}
	else {filename1 = "geo_DOff";}
	filename2 ="_info.txt"; 
	AnaInput.open( (filedir+filename1+filename2).c_str() ,ios::out);
	AnaInput << "R_1" << "\t" << "horilines" << "\t" << "vertilines" << "\t" << "apertsizeX" << "\t" << "apertsizeY" << endl;
	AnaInput << R_1 << "\t" << horilines << "\t" << vertilines << "\t" << apertsizeX << "\t" << apertsizeY << endl;
	AnaInput.close();

	//parameters
	//double momentum = 1187.29; //keV/c
	double c= 299792458; // light speed
	double thetaEm = 5./180.*M_PI; // emission angle for f()
	double B[4];          // B-field components [T]
	double P[4];          // Descartes coordinates [m]
	double b;             // Absolute B-field [T]
	double Btang, Brad;
	int counter =0;

	// Prepare simulation with magfield 3:
	if (prep==1) {
		input_coils(inputcoil);
		test_coils(filedir);
		magsource(filedir);
	}

	cout << endl<< "GEO+DRIFT: f()= " << (cos(thetaEm)+1./cos(thetaEm))/2. << endl;
	
	double startP[2][horilines*vertilines];
	double xstep = apertsizeX/(double)(horilines-1);
	double ystep = apertsizeY/(double)(vertilines-1);

	// define start points dependent on apertsize and hori/verti-lines	
	for(int x=1;x<=horilines; x++){
		for(int y=1;y<=vertilines; y++){
			
			startP[0][counter]= -apertsizeX/2. + (x-1)*xstep;
			startP[1][counter]= apertYshift -apertsizeY/2. + (y-1)*ystep;
		 	counter++;
		}
	}

	double radius_bend;
	int steps = 1000;
	double anglestep = alpha/(double)(steps-1);	

	cout  <<  "GEO+DRIFT: Drift = " << driftOn << endl;

	// Save the B-field lines along the beam axis
	for (int  line = 0; line < horilines*vertilines ; line++) { // different lines shifted radial
	    	
		// open outputfile
	    	std::string fileoutname;
	    	std::stringstream fieldstream;
	    	fieldstream << filename1 <<"_X"<< startP[0][line] << "_Y" << startP[1][line]  << ".txt";
	    	fieldstream >> fileoutname;
	    	fileoutname = filedir + fileoutname;
	    	ofstream fileout;
	    	fileout.open (fileoutname.c_str(), ios::out);
	    	// file header positions:
	    	fileout << "pos_x [m]" << "\t" << "pos_y[m]" << "\t" << "pos_z[m]";
	    	// file header B-field values
	    	fileout <<"\t" << "B-Field" <<"\t"<< "B_x" <<"\t"<< "B_y" <<"\t"<< "B_z" <<"\t"<< "B_t" <<"\t"<< "B_r"<< "\n";
	
	    	// radius for specific cut lines
	    	radius_bend = R_1 + startP[1][line];
	    	P[1] = startP[0][line];
		// status report
	    	std::cout << "GEO+DRIFT: line: " << line +1 << std::endl;

		counter =0;

	    	// calculate B-field values at every step on the line:
	    	for(int  angle= 0; angle < steps; angle++) {
		  	P[2] = radius_bend*cos(angle*anglestep);
	        	P[3] = radius_bend*sin(angle*anglestep);
	      		// calculating the B-field with Magfield 3
	      		magfield_elliptic(P,B,filedir);
	      		b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // claculating abs. B-field
			
			Btang = -sin(angle*anglestep)*B[2] + cos(angle*anglestep)*B[3];
			Brad = cos(angle*anglestep)*B[2] + sin(angle*anglestep)*B[3];
	      		// store the datapoint's position
	      		fileout << P[1] <<"\t"<< P[2]<<"\t"<< P[3];
	      		// store the B-field values in the datapoint
	      		fileout <<"\t"<< b <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< B[3] <<"\t" << Btang <<"\t" << Brad << "\n" ;

			// projection of 2D field components y,z onto normed tangential vector
	   		if(driftOn){ P[1] = P[1] - momentum*1000./c/Btang*anglestep*(cos(thetaEm)+1./cos(thetaEm))/2. ;}
			else {}
			counter ++;
	      	}
	fileout.close();
	}
	  
}
