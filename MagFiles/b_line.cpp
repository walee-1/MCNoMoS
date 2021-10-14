void b_line(double R_1,double alpha,  double start_posZ, int horilines, int vertilines, double apertsizeX, double apertsizeY, double apertYshift,double apertXshift, string inputcoil, string filedir, int prep,int driftdir,double momentum, bool driftOn, bool BxDrift, double detposZ)
{
	//parameters: R1, startZ, .., .. , .. , .., inputcoil, filedir, magfield prep, maxdrift, driftOn: X Drift nach winkel dazu addieren, projection ON/OFF: X = xdrift+BX oder ohne BX, detposZ
	
	if (driftOn == false && BxDrift == false) {cout << "B-LINE: meh" << endl;}

	//parameters
	double c= 299792458; // light speed
	double thetaEm = 5./180.*M_PI; // emission angle for f()
	double B[4];          // B-field components [T]
	double P[4];          // Descartes coordinates [m]
	double RxBEnd[4];
	double P_buffer[4];   // Descartes coordinates [m]
	double b, Btang, Brad;             // Absolute B-field [T]
	double angle, angleold;

	// Adjustable:
	double step_size = 5e-4;  // point distance in [m]
	int n_loops = 10;         // number of loops, when values are written
	int i_loop;	

	// calc startpoints
	double startP[2][horilines*vertilines];
	int counter=0;
	double xstep = apertsizeX/(double)(horilines-1);
	double ystep = apertsizeY/(double)(vertilines-1);

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

	
	// change file prefix for setting of driftOn and BxDrift
	string prefix;
	if(driftOn) {
		if(BxDrift){prefix = "bline_DOn_BxDOn_";}
		else	{prefix = "bline_DOn_BxDOff_";}
	}
	else {prefix = "bline_DOff_BxDOn_";}

	// first write out file for analysis input
	ofstream AnaInput;
	AnaInput.open( (filedir+prefix+"info.txt").c_str() ,ios::out);
	AnaInput << "R_1" << "\t" << "horilines" << "\t" << "vertilines" << "\t" << "apertsizeX" << "\t" << "apertsizeY" << "\t" << "driftOn" << "\t" << "BxDrift" << endl;
	AnaInput << R_1 << "\t" << horilines << "\t" << vertilines << "\t" << apertsizeX << "\t" << apertsizeY << "\t" << driftOn<<"\t" << BxDrift << endl;
	AnaInput.close();
	
	cout << endl << "B-LINE: Drift = " << driftOn << " BxDrift = " << BxDrift << endl;

	// Save the B-field lines along the beam axis
	for ( int line =0; line <horilines*vertilines; line++) {
		// open outputfile
		std::string fileoutname;
		std::stringstream fieldstream;
		fieldstream << prefix << "X" << startP[0][line] << "_Y"<< startP[1][line] << ".txt";
		fieldstream >> fileoutname;
		fileoutname = filedir + fileoutname;
		ofstream fileout;
		fileout.open (fileoutname.c_str(), ios::out);
	    	
		// file header positions:
	    	fileout << "pos_x [m]" << "\t" << "pos_y[m]" << "\t" << "pos_z[m]";
	    	// file header B-field values
	    	fileout <<"\t" << "B-Field" <<"\t"<< "B_x" <<"\t"<< "B_y" <<"\t"<< "B_z" <<"\t"<< "B_t" <<"\t"<< "B_r"<< "\n";
		
		// status report
		std::cout << "B-LINE: line = " << line +1 << std::endl;
		
		// starting point
		P[1] = startP[0][line] + apertXshift;
		P[2] = R_1 + startP[1][line] + apertYshift;
		P[3] = start_posZ;
		counter = 0;
			
		angleold = 0.;
		angle =0.;
		i_loop = n_loops;   // loop counter, final value at beginning so that first point is definitely written 
		// calculate B-field values at every step on the cut line:
		while (!(P[2] < 0. && P[3] < detposZ)) {
			
			// calculating the B-field with Magfield 3
			magfield_elliptic(P,B,filedir);
			b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // claculating abs. B-field
			
			// bline starts not just at RxB but before. Apply RxB drift only in RxB!
			//calc angle only in RxB, Btang and Brad are different besides the RxB
			if(P[3] > 0.){
				angle = atan2(P[3],P[2]);
				Btang = -sin(angle)*B[2] + cos(angle)*B[3];
				Brad = cos(angle)*B[2] + sin(angle)*B[3];
			}
			else{
				Btang = B[3];
				Brad = sqrt(B[2]*B[2]+B[1]*B[1]);
			}
	
			if (i_loop == n_loops) {
	      			// store the datapoint's position
		      		fileout << P[1] <<"\t"<< P[2]<<"\t"<< P[3];
		      		// store the B-field values in the datapoint
	      			fileout <<"\t"<< b <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< B[3] <<"\t" << Btang <<"\t" << Brad << "\n" ;
				i_loop = 0;
			} else {i_loop++;}
			
			// store old datapoint
			P_buffer[1] = P[1];
			P_buffer[2] = P[2];
			P_buffer[3] = P[3];
			
			// calculate new position
			P[2] += B[2]/b*step_size;
			P[3] += B[3]/b*step_size;



			// Detector break! //////////////////
			if(P[2] < 0. && P[3] < detposZ ) {cout << "BLINE: finished at Detector" << endl;}

			//if P[3] gets over det positions, while loop is finished and last point will not be written away -> manual write away last point before that
			// but if last point of arc is already written away because of i_loop, then do not write away, because the following point is already out of arc
			if(P[2] < 0. && P[3] < detposZ && i_loop != 0) {
	      			// store the datapoint's position
		      		fileout << P_buffer[1] <<"\t"<< P_buffer[2]<<"\t"<< P_buffer[3];
		      		// store the B-field values in the datapoint
	      			fileout <<"\t"<< b <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< B[3] <<"\t" << Btang <<"\t" << Brad << "\n" ;
				break;
			}
			/////////////////////////////////////
			

			counter = counter +1;
			
			// for x, calc angle to find out, what drift should be added (estimation)
			// x component
			
			// again, make case for RxB or else, due to drift only applying in RxB
			if(P[3] > 0.){
	
				if (driftOn == true && BxDrift == false){ // drift added to Bx "drift"
					P[1] += B[1]/b*step_size + (double)driftdir*momentum*1000./c/Btang * (angleold - angle) * (cos(thetaEm)+1./cos(thetaEm))/2.;
				}
				else if (driftOn == false){ // only "drift" effect through Bx component
					P[1] += B[1]/b*step_size;
				}
				else if (driftOn == true && BxDrift == true){ // only drift without effect from Bx
					P[1] += (double)driftdir*momentum*1000./c/Btang * (angleold - angle) * (cos(thetaEm)+1./cos(thetaEm))/2.;	
				}
			}
			else if(BxDrift == true){P[1] += B[1]/b*step_size;}
			else{}
	
			angleold = angle;
			// break conditions for outside area
			if (P[3] < -13. ) {cout << "B-line: -z outside" << endl; break;}
			if (P[3] > R_1+0.5 ) {cout << "B-line: z outside" << endl; break;}
			if (P[2] > R_1+0.5 ) {cout << "B-line: y outside" << endl; break;}
			if (P[2] < -R_1-0.5 ) {cout << "B-line: -y outside" << endl; break;}
			if (P[1] < -0.5 ) {cout << "B-line: -x outside" << endl; break;}
			if (P[1] > 0.5 ) {cout << "B-line: x outside" << endl; break;}
		}
		cout << "B-LINE: Iterations = " <<  counter << endl;
		fileout.close();
	}

}
