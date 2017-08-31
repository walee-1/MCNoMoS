void b_line_real(double R_1,double alpha,  double start_posZ, int horilines, int vertilines, double apertsizeX, double apertsizeY, double apertYshift,double apertXshift, string inputcoil, string filedir, int prep,int driftdir,double momentum,int driftOn, double detposZ)
{
	//parameters: R1, startZ, .., .. , .. , .., inputcoil, filedir, magfield prep, maxdrift, driftOn: X Drift nach winkel dazu addieren, projection ON/OFF: X = xdrift+BX oder ohne BX, detposZ
	
	//parameters
	double c= 299792458; // light speed
	double thetaEm = 10./180.*M_PI; // emission angle for f()
	double B[4];          // B-field components [T]
	double P[4];          // Descartes coordinates [m]
	double RxBEnd[4];
	double P_buffer[4];   // Descartes coordinates [m]
	double B_old[4];
	double planevec[4];
	double planevec_length, n_length, n_old_length;
	double n_old[4];
	double n[4];
	double b, Btang, Brad, b_old, R_curve, R_min;
	double angle;
	double v_drift, denomi, n_diff;
	int max_coord;
	int button = 0;
	
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
	if(driftOn == 1) {prefix = "bline_D1stOn_";}
	else if(driftOn == 2){prefix = "bline_D2ndOn_";}
	else{prefix = "bline_Doff_";}

	// first write out file for analysis input
	ofstream AnaInput;
	AnaInput.open( (filedir+prefix+"info.txt").c_str() ,ios::out);
	AnaInput << "R_1" << "\t" << "horilines" << "\t" << "vertilines" << "\t" << "apertsizeX" << "\t" << "apertsizeY" << "\t" << "driftOn" << endl;
	AnaInput << R_1 << "\t" << horilines << "\t" << vertilines << "\t" << apertsizeX << "\t" << apertsizeY << "\t" << driftOn << endl;
	AnaInput.close();
	
	cout << endl << "B-LINE: Drift = " << driftOn << endl;

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
			
		angle =0.;
		R_min = 100.;
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
			B_old[1] = B[1];
			B_old[2] = B[2];
			B_old[3] = B[3];
			b_old = b;


			// calculate new position by following field line
			P[1] += B[1]/b*step_size;
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
		

			//calc 2nd point B
			magfield_elliptic(P,B,filedir);
			b = sqrt( B[1]*B[1] + B[2]*B[2] + B[3]*B[3] );

			/*
			if( P_buffer[3] > 0.4 && button == 0 ){
				cout << "B-LINE BUG: data of first 2 points in RxB with coords and Bs:" << endl;
				cout << "Cx Cy Cz, CBx CBy CBz, Dx Dy Dz, DBx DBy DBz " << endl;
				cout << P_buffer[1] << "\t" << P_buffer[2] << "\t" << P_buffer[3] << "\t" << B_old[1] << "\t" << B_old[2] << "\t" << B_old[3];
				cout << "\t" << P[1] << "\t"<< P[2] << "\t" << P[3] << "\t" << B[1] << "\t" << B[2] << "\t" << B[3] << endl;

				button ++;
			}
			*/

			////////// Drift calculation (also in other areas, where curvature can happen)
			// Bx stream line follow + analytical drift added (calc curvature)
			if ( driftOn == 1 || driftOn == 2 ){ // drift added to field line (1st or 2nd order)

				// create perpendicular normalized vector of the two b-field vectors to create plane, in which curvature radius is calculated
				planevec[1] = (B_old[2]*B[3] - B_old[3]*B[2]);
				planevec[2] = (B_old[3]*B[1] - B_old[1]*B[3]);
				planevec[3] = (B_old[1]*B[2] - B_old[2]*B[1]);
				planevec_length = sqrt( planevec[1]*planevec[1] + planevec[2]*planevec[2] + planevec[3]*planevec[3] );
				planevec[1] = planevec[1]/planevec_length;
				planevec[2] = planevec[2]/planevec_length;
				planevec[3] = planevec[3]/planevec_length;
	
				// if plane vec is very close to ZERO, B and B_old are parallel -> no curvature
				if (  planevec_length > 1.e-6 ){

					// perpend vec of 2nd point pointing to curvature center
					n[1] = (planevec[2]*B[3] - planevec[3]*B[2]);
					n[2] = (planevec[3]*B[1] - planevec[1]*B[3]);
					n[3] = (planevec[1]*B[2] - planevec[2]*B[1]);
					n_length = sqrt( n[1]*n[1] + n[2]*n[2] + n[3]*n[3] );
					n[1] = n[1]/n_length;
					n[2] = n[2]/n_length;
					n[3] = n[3]/n_length;

					// 1st perpend vec, pointing to curv center
					n_old[1] = (planevec[2]*B_old[3] - planevec[3]*B_old[2])/b_old;
					n_old[2] = (planevec[3]*B_old[1] - planevec[1]*B_old[3])/b_old;
					n_old[3] = (planevec[1]*B_old[2] - planevec[2]*B_old[1])/b_old;
					n_old_length = sqrt( n_old[1]*n_old[1] + n_old[2]*n_old[2] + n_old[3]*n_old[3] );
					n_old[1] = n_old[1]*n_old_length;
					n_old[2] = n_old[2]*n_old_length;
					n_old[3] = n_old[3]*n_old_length;

					// for 1st ORDER, we dont actually need R, but only the angle between the 2 perpend. vectors
					angle = acos( n[1]*n_old[1] + n[2]*n_old[2] + n[3]*n_old[3] );

					
					if( driftOn == 1 ){

						v_drift = (double)driftdir*momentum*1000./c/b_old *angle  *(cos(thetaEm)+1./cos(thetaEm))/2.;
	
						P[1] += planevec[1]*v_drift;
						P[2] += planevec[2]*v_drift;
						P[3] += planevec[3]*v_drift;
					}


					if ( driftOn == 2 ){

						//here we have to build in a security measure because if denominator is close to equal, R_curve -> Infinity
						denomi = n_old[2]*n[3] - n[2]*n_old[3];
						if ( fabs( denomi ) > 1.e-12 ){
							
							R_curve = (P_buffer[3] - P[3])*n[2] - (P_buffer[2] - P[2])*n[3];
							R_curve = R_curve / denomi;

							//cout << "B-LINE BUG: R_curve calc = " << R_curve << endl;
							if (R_curve < R_min ) R_min = R_curve;
	
							// now calculate RxB drift velocity 
							// we use here the velocity formula for the drift calc, and we add them with planevec
							v_drift = (double)driftdir*momentum*1000./c/b_old *step_size/R_curve  *(cos(thetaEm)+1./cos(thetaEm))/2.;
	
							P[1] += planevec[1]*v_drift;
							P[2] += planevec[2]*v_drift;
							P[3] += planevec[3]*v_drift;
	
						}
						else {
							R_curve = 1000.; // maximum radius chosen
						}
					} // 2nd drift close
				} // if curvature
			} // 1st or 2nd drift

	
			// break conditions for outside area
			if (P[3] < -13. ) {cout << "B-line: -z outside" << endl; break;}
			if (P[3] > R_1+0.5 ) {cout << "B-line: z outside" << endl; break;}
			if (P[2] > R_1+0.5 ) {cout << "B-line: y outside" << endl; break;}
			if (P[2] < -R_1-0.5 ) {cout << "B-line: -y outside" << endl; break;}
			if (P[1] < -0.5 ) {cout << "B-line: -x outside" << endl; break;}
			if (P[1] > 0.5 ) {cout << "B-line: x outside" << endl; break;}
		}
		//cout << "B-LINE: R_min = " << R_curve << endl;
		cout << "B-LINE: Iterations = " <<  counter << endl;
		fileout.close();
	}

}
