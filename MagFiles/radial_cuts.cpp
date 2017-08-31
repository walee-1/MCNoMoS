
void radial_cuts(int n_lines, int N_coil, double line_end, double step_size, std::string inputcoil, std::string filedir, int prep)
{
  //parameters
   double B[4];          // B-field components [T]
   double P[4];          // Descartes coordinates [m]
   double b;             // Absolute B-field [T]
   double angle;
     double err[4],phi,x,y,z;

// Simulate with magfield 3:
if(prep==1){
   input_coils(inputcoil);
   test_coils(filedir);
   magsource(filedir);
   }
//    outputcoil();

// Save the B-field lines along the beam axis

for(int i_line = 0; i_line <= n_lines; i_line++)        // different lines
   {
   P[1] = 0.;
   if(i_line==n_lines){angle = M_PI/((double)N_coil-1.)*((double)((int)(0.5*((double)N_coil-1.))));}
   else{angle = M_PI/((double)N_coil-1.)*((double)((int)((double)i_line/((double)n_lines-1.)*((double)N_coil-1.)))+0.5);}                 // calculate the angle of the cut line
   //cout << "angle: " << angle << " = " << (int)(angle*180./M_PI*100.) << " sin(angle): " << sin(angle) <<endl;
   std::string fileoutname;
   std::stringstream fieldstream;
   fieldstream << (int)(angle*180./M_PI*100.)+ 500000;
   fieldstream >> fileoutname;
   fileoutname = filedir + "radial_line_" + fileoutname + ".txt";
   ofstream fileout;
   fileout.open (fileoutname.c_str());
   fileout << "pos_x [m]\t pos_y[m]\t pos_z[m]\t pos[m]\t B-Field\tB_x\tB_y\tB_z\tB_trans\n";
   for(int k = 0; k < line_end/step_size; k++)        // follow the radius outwards
      {
      P[2] = cos(angle)*step_size*k;
      P[3] = sin(angle)*step_size*k;
      magfield_elliptic(P,B,filedir);                      // calculating the B-field
      b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // claculating abs. B-field
      fileout << P[1] << "\t" << P[2]<< "\t"  << P[3]<< "\t"  << sqrt(P[1]*P[1]+P[2]*P[2]+P[3]*P[3])<< "\t"  ;
      fileout << b<< "\t"  << B[1]<< "\t"  << B[2]<< "\t"  <<B[3]<< "\t" <<fabs(B[1])<< "\n" ;
      }
      cout << "RADIAL: ArcTan= " << atan2(P[3],P[2]) << endl;
    fileout.close();
    }

}
