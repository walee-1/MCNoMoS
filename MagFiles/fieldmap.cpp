void fieldmap(double * StartP, double * Length, double step_size, string inputcoil, string filedir, int prep){
 
       	//parameters
  double B[4];          // B-field components [T]
  double P[4];          // Descartes coordinates [m]
  double b;             // Absolute B-field [T]
    double err[4],phi,x,y,z;    // not used
    int i_x, i_y, i_z;

// Prepare simulation with magfield 3:
  if(prep==1){
    input_coils(inputcoil);
    test_coils(filedir);
    magsource(filedir);
    }

std::cout <<  endl<<  "MAP: Calculate fieldmap...  ";

// open outputfile
  std::string fileoutname;
  std::stringstream fieldstream;
  fieldstream << step_size*1000. << "mm.txt";
  fieldstream >> fileoutname;
  fileoutname = filedir + "fieldmap_" + fileoutname;
  ofstream fileout;
  fileout.open (fileoutname.c_str(), ios::out);
  // file header positions:
  fileout <<setiosflags(ios::left)<<setw(14)<< "pos_x [m]" <<setw(14)<< "pos_y [m]" <<setw(14)<< "pos_z [m]";
  // file header B-field values
  fileout <<setiosflags(ios::left)<<setw(16)<< "B-Field [T]" <<setw(16)<< "B_x [T]" <<setw(16)<< "B_y [T]" <<setw(16)<< "B_z [T]\n";

  // Save the B-field lines along the beam axis
  for( i_z = 0; i_z < Length[2]/step_size; i_z++) {       // different lines
    cout << "MAP: i_z = " << i_z << endl;
    for( i_y = 0; i_y < Length[1]/step_size; i_y++) {
   	for( i_x = 0; i_x < Length[0]/step_size; i_x++) {

	    // set coordinates of the datapoint
      		P[1] = StartP[0] + (double)i_x*step_size;
	      P[2] = StartP[1] + (double)i_y*step_size;
	      P[3] = StartP[2] + (double)i_z*step_size;

	      	magfield_elliptic(P,B,filedir);         // calculating the B-field
		b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // calculating abs. B-field
		// store the datapoint's position
      		fileout <<setiosflags(ios::left)<<setw(14)<< P[1] <<setw(14)<< P[2]<<setw(14)<< P[3]<<setw(14);
      		// store the B-field values in the datapoint
      		fileout <<setiosflags(ios::left)<<setw(16)<< b <<setw(16)<< B[1] <<setw(16)<< B[2] <<setw(16)<< B[3] << "\n" ;
      	}
    }
  }
 
  fileout.close();

	cout << "Finished" << endl;
//return 0;
}
