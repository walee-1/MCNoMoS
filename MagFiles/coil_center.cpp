void coil_center(double R_1, double angle, double start_pos, double bend_start, double length_detect,int horilines,int vertilines,double apertsizeX,double apertsizeY, double apertYshift,double apertXshift, string inputcoil, string filedir, int prep)
{
	//parameters
	double B[4];          // B-field components [T]
	double P[4];          // Descartes coordinates [m]
	double b;             // Absolute B-field [T]
	
	// Adjustable:
	double lin_step_size = 10e-3;     // point distance in [m]
	int angle_steps = 500;
	double angle_step_size = angle/(double)(angle_steps);

	// calc startpoints
	double startP[2][horilines*vertilines];
	int counter=0;
	double xstep = apertsizeX/(double)(horilines-1);
	double ystep = apertsizeY/(double)(vertilines-1);
	double radius_bend;

	// define start points dependent on apertsize and hori/verti-lines	
	for(int x=1;x<=horilines; x++){
		for(int y=1;y<=vertilines; y++){
			
			startP[0][counter]= -apertsizeX/2. + (x-1)*xstep;
			startP[1][counter]= -apertsizeY/2. + (y-1)*ystep;
		 	counter++;
		}
	}
	
	// Prepare simulation with magfield 3:
	if (prep==1) {
		input_coils(inputcoil);
		test_coils(filedir);
		magsource(filedir);
	}

	// first write out file for analysis input
	ofstream AnaInput;
	AnaInput.open( (filedir+"center_info.txt").c_str() ,ios::out);
	AnaInput << "R_1" << "\t" << "horilines" << "\t" << "vertilines" << "\t" << "apertsizeX" << "\t" << "apertsizeY" << endl;
	AnaInput << R_1 << "\t" << horilines << "\t" << vertilines << "\t" << apertsizeX << "\t" << apertsizeY << endl;
	AnaInput.close();
	

	// Save the B-field lines along the beam axis
	for ( int line =0; line <horilines*vertilines; line++) {
		// open outputfile
		std::string fileoutname;
		std::stringstream fieldstream;
		fieldstream << "center_" << "X" << startP[0][line] << "_Y"<< startP[1][line] << ".txt";
		fieldstream >> fileoutname;
		fileoutname = filedir + fileoutname;
		ofstream fileout;
		fileout.open (fileoutname.c_str(), ios::out);
	    	
		// file header positions:
	    	fileout << "pos_x [m]" << "\t" << "pos_y[m]" << "\t" << "pos_z[m]";
	    	// file header B-field values
	    	fileout <<"\t" << "B-Field" <<"\t"<< "B_x" <<"\t"<< "B_y" <<"\t"<< "B_z" << "\n";
		
		// status report
		std::cout << "Center: line = " << line +1 << std::endl;
		// starting point
		P[1] = apertXshift + startP[0][line];
		P[2] = apertYshift + R_1 + startP[1][line];
		P[3] = start_pos;

		// Save the B-field lines along the beam axis	
		//Inlet Loop
		while(P[3] < bend_start){
			magfield_elliptic(P,B,filedir);
			b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // claculating abs. B-field
	      		// store the datapoint's position
	      		fileout << P[1] <<"\t"<< P[2]<<"\t"<< P[3];
	      		// store the B-field values in the datapoint
	      		fileout <<"\t"<< b <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< B[3] <<"\n" ;
			//calc next step		
			P[3] = P[3] + lin_step_size;
		}

		radius_bend=P[2];
		P[3]=0.;
		//RxB Loop
		int anglestep;
		for(anglestep=0;anglestep < angle_steps; anglestep++){

			P[2] = radius_bend*cos(anglestep*angle_step_size);
			P[3] = radius_bend*sin(anglestep*angle_step_size);
			
			magfield_elliptic(P,B,filedir);
			b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // claculating abs. B-field
	      		// store the datapoint's position
	      		fileout << P[1] <<"\t"<< P[2]<<"\t"<< P[3];
	      		// store the B-field values in the datapoint
	      		fileout <<"\t"<< b <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< B[3] <<"\n" ;
		}

		//OUtlet Loop
		P[2] = radius_bend*cos(angle);
		P[3] = radius_bend*sin(angle);
		double RxBEndYZ[3];
		RxBEndYZ[2] = P[2];
		RxBEndYZ[3] = P[3];
		int counter =0;
		int outlet_steps = length_detect / lin_step_size + 1.;
		for(int i=0; i < outlet_steps; i++){
			magfield_elliptic(P,B,filedir);
			b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // claculating abs. B-field
	      		// store the datapoint's position
	      		fileout << P[1] <<"\t"<< P[2]<<"\t"<< P[3];
	      		// store the B-field values in the datapoint
	      		fileout <<"\t"<< b <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< B[3] <<"\n" ;
			
			//calc next step
			P[2] = RxBEndYZ[2] - i*lin_step_size*sin(M_PI - angle);
			P[3] = RxBEndYZ[3] - i*lin_step_size*cos(M_PI - angle);
			counter ++;
		}
		//cout << "CENTER: counter = " << counter << endl;
		fileout.close();
	}
}
