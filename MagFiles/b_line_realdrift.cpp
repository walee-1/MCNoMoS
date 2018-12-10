void b_line_real(bool NoMoSOn,double R_1,double alpha,  double start_posZ,bool onlyCornerCenter, int horilines, int vertilines, double apertsizeX, double apertsizeY, double apertYshift,double apertXshift, string inputcoil, string filedir, int prep,int driftdir,double momentum,int driftOn, double detposZ)
{
	//parameters: R1, startZ, .., .. , .. , .., inputcoil, filedir, magfield prep, maxdrift, driftOn: X Drift nach winkel dazu addieren, projection ON/OFF: X = xdrift+BX oder ohne BX, detposZ
	
	//parameters
	double ffactor(double);
	double localtheta(double, double, double);
	double c= 299792458; // light speed
	double thetaEm = 20./180.*M_PI; // emission angle for f()
	double B[4];          // B-field components [T]
	double P[4];          // Descartes coordinates [m]
	double RxBEnd[4];
	double P_old[4];   // Descartes coordinates [m]
	double B_old[4];
	double m_e = 510998.;
	double planevec[4];
	double planevec_B[4];
	double planevec_B_length, planevec_length, n_length, n_old_length, n_new_length;
	double n_old[4];
	double n[4];
	double n_new[4];
	double B_new[4];
	double P_new[4], diff_vec[4];
   double RxB_vec[4], RxB_vec_new[4];
	double v1st_A, v1st_B;
	double b, b_old, b_new, R_A, R_B, startB;
	double angle;
	double drift1st, denomi;
	int max_coord;
	int button = 0;

	
	// Adjustable:
	double step_size = 5.e-4;  // point distance in [m]
	int n_loops = 10;         // number of loops, when values are written
	int i_loop;	

	// calc startpoints
	int lines;
	if(onlyCornerCenter == true) lines =5;
	else lines = horilines*vertilines;
	double startP[2][lines];
	int counter=0;
	
	// differentiation if onlyCornerCenter == true
	if( onlyCornerCenter == true){
		//first the center point
		startP[0][counter] = 0.; startP[1][counter] = 0.;
		counter++;
		
		for(int x=-1; x <2; x=x+2){
			for(int y=-1; y<2; y=y+2){
				startP[0][counter]=x*apertsizeX/2.;
				startP[1][counter]=y*apertsizeY/2.;
				counter ++;
			}
		}
		cout << "BLINEBUG: startPs" << endl;
		for(int i=0;i<5;i++) cout << startP[0][i] << "\t" << startP[1][i] << endl;
	}
	else{
		double xstep, ystep;
        	if (horilines == 1) xstep = 0;
        	else xstep = apertsizeX/(double)(horilines-1);
		if( vertilines == 1) ystep = 0;
		else ystep = apertsizeY/(double)(vertilines-1);

		// define start points dependent on apertsize and hori/verti-lines	
		for(int x=1;x<=horilines; x++){
			for(int y=1;y<=vertilines; y++){
				
				startP[0][counter]= -apertsizeX/2. + (x-1)*xstep;
				startP[1][counter]= -apertsizeY/2. + (y-1)*ystep;
			 	counter++;
			}
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
	AnaInput << "R_1" << "\t" << "horilines" << "\t" << "vertilines" << "\t" << "apertsizeX" << "\t" << "apertsizeY" << "\t" << "driftOn" << "\t" << "onlyCornerCenter" << endl;
	AnaInput << R_1 << "\t" << horilines << "\t" << vertilines << "\t" << apertsizeX << "\t" << apertsizeY << "\t" << driftOn << "\t" << onlyCornerCenter << endl;
	AnaInput.close();
	
	cout << endl << "B-LINE: Drift = " << driftOn << endl;

	// Save the B-field lines along the beam axis
	for ( int line =0; line <lines; line++) {
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
	    	fileout <<"\t" << "B-Field" <<"\t"<< "B_x" <<"\t"<< "B_y" <<"\t"<< "B_z" << "\n";
		
		// status report
		std::cout << "B-LINE: line = " << line +1 << std::endl;
		
		// starting point
		P[1] = startP[0][line] + apertXshift;
		P[2] = R_1 + startP[1][line] + apertYshift;
		P[3] = start_posZ;
		counter = 0;
		cout << "BlineBUG: startP = " << P[1] << "\t" << P[2] << "\t" << P[3] << endl;
			
		angle =0.;
		i_loop = n_loops;   // loop counter, final value at beginning so that first point is definitely written 
		// calculate B-field values at every step on the cut line:
		bool atDet= false;
		while ( !atDet ) {
			
			// calculating the B-field
			magfield_elliptic(P,B,filedir);
			b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // claculating abs. B-field
			
			// get startB
			if( counter == 0 ) startB = b;

			// write out current data point if ..
			if (i_loop == n_loops) {
	      			// store the datapoint's position
		      		fileout << P[1] <<"\t"<< P[2]<<"\t"<< P[3];
		      		// store the B-field values in the datapoint
	      			fileout <<"\t"<< b <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< B[3] << "\n" ;
				i_loop = 0;
			} else {i_loop++;}
			
			// store old datapoint, lets say point A
			P_old[1] = P[1];
			P_old[2] = P[2];
			P_old[3] = P[3];
			B_old[1] = B[1];
			B_old[2] = B[2];
			B_old[3] = B[3];
			b_old = b;


			// calculate new position by following field line, call it point B
			P[1] += B[1]/b*step_size;
			P[2] += B[2]/b*step_size;
			P[3] += B[3]/b*step_size;

			// Detector break! //////////////////
			if ( P[3] <= detposZ && P[2] < 0) atDet = true;
			if( atDet ) {cout << "BLINE: finished at Detector" << endl;}

			//if P[3] gets over det positions, while loop is finished and last point will not be written away -> manual write away last point before that
			// but if last point of arc is already written away because of i_loop, then do not write away, because the following point is already out of arc
			if( atDet && i_loop != 0) {
	      			// store the datapoint's position
		      		fileout << P_old[1] <<"\t"<< P_old[2]<<"\t"<< P_old[3];
		      		// store the B-field values in the datapoint
	      			fileout <<"\t"<< b <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< B[3] << "\n" ;
				break;
			}
			/////////////////////////////////////
			counter = counter +1;
		
			//calc 2nd point Bfield
			magfield_elliptic(P,B,filedir);
			b = sqrt( B[1]*B[1] + B[2]*B[2] + B[3]*B[3] );



			////////////////////////////////
			////////// Drift calculation (also in other areas, where curvature can happen)
			if ( driftOn == 1 || driftOn == 2 ){ // drift added to field line (1st or 2nd order)

				// create perpendicular normalized vector of the two b-field vectors to create plane, in which curvature radius is calculated
				planevec[1] = (B_old[2]*B[3] - B_old[3]*B[2]);
				planevec[2] = (B_old[3]*B[1] - B_old[1]*B[3]);
				planevec[3] = (B_old[1]*B[2] - B_old[2]*B[1]);
				planevec_length = sqrt( planevec[1]*planevec[1] + planevec[2]*planevec[2] + planevec[3]*planevec[3] );
	
				// if plane vec is very close to ZERO, B and B_old are parallel -> no curvature
				if (  planevec_length > 1.e-9 ){

					// perpend vec of 2nd point pointing to curvature center
					n[1] = (planevec[2]*B[3] - planevec[3]*B[2]);
					n[2] = (planevec[3]*B[1] - planevec[1]*B[3]);
					n[3] = (planevec[1]*B[2] - planevec[2]*B[1]);
					n_length = sqrt( n[1]*n[1] + n[2]*n[2] + n[3]*n[3] );
					n[1] = n[1]/n_length;
					n[2] = n[2]/n_length;
					n[3] = n[3]/n_length;

					// 1st perpend vec, pointing to curv center -> is also negative unit vector of R
					n_old[1] = (planevec[2]*B_old[3] - planevec[3]*B_old[2]);
					n_old[2] = (planevec[3]*B_old[1] - planevec[1]*B_old[3]);
					n_old[3] = (planevec[1]*B_old[2] - planevec[2]*B_old[1]);
					n_old_length = sqrt( n_old[1]*n_old[1] + n_old[2]*n_old[2] + n_old[3]*n_old[3] );
					n_old[1] = n_old[1]/n_old_length;
					n_old[2] = n_old[2]/n_old_length;
					n_old[3] = n_old[3]/n_old_length;

					// for 1st ORDER, we dont actually need R, but only the angle between the 2 perpend. vectors
					angle = acos( n[1]*n_old[1] + n[2]*n_old[2] + n[3]*n_old[3] );
                    
		                        // now calc the NON normed RxB vector for the drift distance with normed R vector (-n_old) and B not normed -> additional 1/B in drift1st
		                        // minus because n_old points towards center, R away from center
		                        RxB_vec[1] = -(n_old[2]*B_old[3] - n_old[3]*B_old[2]);
		                        RxB_vec[2] = -(n_old[3]*B_old[1] - n_old[1]*B_old[3]);
		                        RxB_vec[3] = -(n_old[1]*B_old[2] - n_old[2]*B_old[1]);

					drift1st = (double)driftdir*momentum*1000./c/b_old/b_old *angle  * ffactor(localtheta(thetaEm, startB, b_old));
					//cout << "BLINEBUG P_old[3] = " << P_old[3] << endl;
					//cout << "BLine-BUG: Drift1st = " << drift1st << endl;
	
					// this is now point B'
					P[1] += RxB_vec[1]*drift1st;
					P[2] += RxB_vec[2]*drift1st;
		                        P[3] += RxB_vec[3]*drift1st;
				
					
					// 2nd Order needs R at point A and B' 
					if ( driftOn == 2 ){

                        			// Bfield of B'
                        			magfield_elliptic(P,B,filedir);
                        			b = sqrt( B[1]*B[1] + B[2]*B[2] + B[3]*B[3] );
                        
						// calc R_A
						//here we have to build in a security measure because if denominator is close to equal, R_curve -> Infinity -> v_1st to ZERO
						denomi = n_old[2]*n[3] - n[2]*n_old[3];
						if ( fabs( denomi ) > 1.e-12 ){
							
							R_A = (P_old[3] - P[3])*n[2] - (P_old[2] - P[2])*n[3];
							R_A = R_A / denomi;

							//cout << "B-LINE BUG: R_A = " << R_A << endl;
	
							// now calculate RxB drift velocity for point A
			                            // m_e is missing here because it cancels out with outer factor of D2nd as well as gamma
			                            // we devide by 2 b_old because one time from Nenner and one time from cross product where we didn't used a normed B vector
							v1st_A = (double)driftdir*momentum*1000.*momentum*1000./b_old/b_old/R_A/c/c *cos(localtheta(thetaEm, startB, b_old)) *ffactor(localtheta(thetaEm, startB, b_old));
							//cout << "BLINE BUG: v1st_A = " << v1st_A << endl;
						} // the direction of this drift velocity is still planevec, so don't overwrite



						// for R_B we still need 1stDrift DISTANCE from B'
						// so this new point is called C
						P_new[1] = P[1] + B[1]/b*step_size;
						P_new[2] = P[2] + B[2]/b*step_size;
						P_new[3] = P[3] + B[3]/b*step_size;

						// Bfield of C
						magfield_elliptic(P_new,B_new,filedir);
						b_new = sqrt( B_new[1]*B_new[1] + B_new[2]*B_new[2] + B_new[3]*B_new[3] );

						// plane vector of B fields of points B' and C
						planevec_B[1] = (B[2]*B_new[3] - B[3]*B_new[2]);
						planevec_B[2] = (B[3]*B_new[1] - B[1]*B_new[3]);
						planevec_B[3] = (B[1]*B_new[2] - B[2]*B_new[1]);
						planevec_B_length = sqrt( planevec_B[1]*planevec_B[1] + planevec_B[2]*planevec_B[2] + planevec_B[3]*planevec_B[3] );
						planevec_B[1] = planevec_B[1]/planevec_B_length;
						planevec_B[2] = planevec_B[2]/planevec_B_length;
						planevec_B[3] = planevec_B[3]/planevec_B_length;

						// if plane vec is very close to ZERO, B and B_old are parallel -> no curvature
						if (  planevec_B_length > 1.e-9 ){

							// perpend vec of point C pointing to curvature center
							n_new[1] = (planevec_B[2]*B_new[3] - planevec_B[3]*B_new[2]);
							n_new[2] = (planevec_B[3]*B_new[1] - planevec_B[1]*B_new[3]);
							n_new[3] = (planevec_B[1]*B_new[2] - planevec_B[2]*B_new[1]);
							n_new_length = sqrt( n_new[1]*n_new[1] + n_new[2]*n_new[2] + n_new[3]*n_new[3] );
							n_new[1] = n_new[1]/n_new_length;
							n_new[2] = n_new[2]/n_new_length;
							n_new[3] = n_new[3]/n_new_length;

							// perpend vec of point B' , pointing to curv center
							// here we simply use n_old because we don't need the old value anymore
							n_old[1] = (planevec_B[2]*B[3] - planevec_B[3]*B[2]);
							n_old[2] = (planevec_B[3]*B[1] - planevec_B[1]*B[3]);
							n_old[3] = (planevec_B[1]*B[2] - planevec_B[2]*B[1]);
							n_old_length = sqrt( n_old[1]*n_old[1] + n_old[2]*n_old[2] + n_old[3]*n_old[3] );
							n_old[1] = n_old[1]/n_old_length;
							n_old[2] = n_old[2]/n_old_length;
							n_old[3] = n_old[3]/n_old_length;

							// for 1st ORDER, we dont actually need R, but only the angle between the 2 perpend. vectors
							angle = acos( n_new[1]*n_old[1] + n_new[2]*n_old[2] + n_new[3]*n_old[3] );
                            
			                            // now calc the NON normed RxB vector for the drift distance with normed R vector (-n_old) and B not normed -> additional 1/B in drift1st
			                            // minus because n_old points towards center, R away from center
			                            RxB_vec_new[1] = -(n_old[2]*B[3] - n_old[3]*B[2]);
			                            RxB_vec_new[2] = -(n_old[3]*B[1] - n_old[1]*B[3]);
			                            RxB_vec_new[3] = -(n_old[1]*B[2] - n_old[2]*B[1]);

							// drift distance for point B'
							drift1st = (double)driftdir*momentum*1000./c/b/b *angle  * ffactor(localtheta(thetaEm, startB, b));
							//cout << "drift1st of B' = " << drift1st << endl;
	
							// this is now point C'
							P_new[1] += RxB_vec_new[1]*drift1st;
							P_new[2] += RxB_vec_new[2]*drift1st;
							P_new[3] += RxB_vec_new[3]*drift1st;
						
							// Bfield of C'
							magfield_elliptic(P_new,B_new,filedir);
							b_new = sqrt( B_new[1]*B_new[1] + B_new[2]*B_new[2] + B_new[3]*B_new[3] );
						}

						// now finally we can calc R_B
						//here we have to build in a security measure because if denominator is close to equal, R_curve -> Infinity -> v_1st to ZERO
						denomi = n_old[2]*n_new[3] - n_new[2]*n_old[3];
						if ( fabs( denomi ) > 1.e-12 ){
							
							R_B = (P[3] - P_new[3])*n_new[2] - (P[2] - P_new[2])*n_new[3];
							R_B = R_B / denomi;

							//cout << "B-LINE BUG: R_B = " << R_B << endl;
	
							// now calculate RxB drift velocity for point B
							v1st_B = (double)driftdir*momentum*1000.*momentum*1000./b/b/R_B/c/c *cos(localtheta(thetaEm, startB, b))  *ffactor(localtheta(thetaEm, startB, b));
							//cout << "v1st_B = " << v1st_B << endl;
						}
                        
			                        // now we need to calc the RxB cross product with normalized R and not normalized Bs for v1st_A and v1st_B
			                        // we use diff_vec here first for v1st_B and substract then 2nd term
			                        diff_vec[1] = v1st_B*RxB_vec_new[1];
			                        diff_vec[2] = v1st_B*RxB_vec_new[2];
			                        diff_vec[3] = v1st_B*RxB_vec_new[3];
                        
                        
						// now we can calculate B'' from B' from 2nd Order Drift DISTANCE
						// this is the velocity difference vectorially calculated with directions of planevec
						diff_vec[1] -= v1st_A*RxB_vec[1];
						diff_vec[2] -= v1st_A*RxB_vec[2];
						diff_vec[3] -= v1st_A*RxB_vec[3];
						
						// then calc the cross product of the diff with Bfield of point A
						diff_vec[1] = (diff_vec[2]*B_old[3] - diff_vec[3]*B_old[2]);
						diff_vec[2] = (diff_vec[3]*B_old[1] - diff_vec[1]*B_old[3]);
						diff_vec[3] = (diff_vec[1]*B_old[2] - diff_vec[2]*B_old[1]);
						
						
						// and now add the full 2nd Order drift to point B' to get B''
       				                 // division by /c/c already happened in v1st and gamma and me were left out because they cancel
						P[1] += (double)driftdir*diff_vec[1]/b_old;
						P[2] += (double)driftdir*diff_vec[2]/b_old;
						P[3] += (double)driftdir*diff_vec[3]/b_old;
						//cout << "B-LINE BUG: 2ndD (y) = " << diff_vec[2]/b_old*m_e/c/c << endl << endl;

					} // 2nd drift close
				} // if curvature
			} // 1st or 2nd drift

	
			// break conditions for outside area
			if (P[3] < -13. ) {cout << "B-line: -z outside" << endl; break;} //bline wrong direction
			if (P[3] > R_1+0.5 ) {cout << "B-line: z outside" << endl; break;} //shot over NoMoS
			if (P[2] > R_1+0.5 ) {cout << "B-line: y outside" << endl; break;} // shot down over NoMoS
			if (P[2] < -R_1-0.5 ) {cout << "B-line: -y outside" << endl; break;} // shot over up NoMoS
			if (P[1] < -0.5 ) {cout << "B-line: -x outside" << endl; break;} //x out
			if (P[1] > 0.5 ) {cout << "B-line: x outside" << endl; break;}
		}
		cout << "B-LINE: Iterations = " <<  counter << endl;
		fileout.close();
	}

}

double localtheta(double starttheta,double startB,double localB){
	return asin( sin(starttheta)* sqrt( localB / startB ) );
}


double ffactor(double localtheta){
	return (cos(localtheta)+1./cos(localtheta))/2.;
}
