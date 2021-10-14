
void coil_center_2(double R_1, double angle, double start_pos, double bend_start, double length_detect, double r_tube, double b_beam, double h_beam, string inputcoil, string filedir, int prep)
{
//parameters
  double B[4];          // B-field components [T]
  double P[4];          // Descartes coordinates [m]
  double b;             // Absolute B-field [T]

// Adjustable:
  int Nr_shift = 3;             // number of lines inside the beam tube
  double step_size = 15e-3;     // point distance in [m]
  double angle_factor = 4.;     // increasing datapoint's density inside RxB

// for calculations:
  double angle_step = step_size/R_1/angle_factor;       // size of angle steps
  double shift_side = 2.*b_beam/((double)Nr_shift-1.);  // size of radial shift
  double shift_up = 2.*h_beam/((double)Nr_shift-1.);    // size of vertical shift
  double steps_app = fabs(start_pos-bend_start)/step_size;  // number of steps inside aperture area
  double steps_RxB = fabs(M_PI/angle_step);                 // number of steps inside the RxB
  double steps_dec = fabs(length_detect)/step_size;         // number steps inside the detector
//std::cout << "APP " << steps_app << " dec " << steps_dec << " RxB " << steps_RxB << std::endl;

// Prepare simulation with magfield 3:
if (prep==1) {
  input_coils(inputcoil);
  test_coils(filedir);
  magsource(filedir);
  }

// Save the B-field lines along the beam axis
for (int i_side = 0; i_side < Nr_shift; i_side++) { // different lines shifted radial
  for (int i_up = 0; i_up < Nr_shift; i_up++) {     // different lines shifted vertically
    // open outputfile
    std::string fileoutname;
    std::stringstream fieldstream;
    fieldstream << i_side*3 + 1 + i_up << ".txt";
    fieldstream >> fileoutname;
    fileoutname = filedir + "coil_center_" + fileoutname;
    ofstream fileout;
    fileout.open (fileoutname.c_str(), ios::out);
    // file header positions:
    fileout <<setiosflags(ios::left)<<setw(14)<< "pos_x [m]" <<setw(14)<< "pos_y[m]" <<setw(14)<< "pos_z[m]" <<setw(8)<< "trash";
    // file header B-field values
    fileout <<setiosflags(ios::left)<<setw(16)<< "B-Field" <<setw(16)<< "B_x" <<setw(16)<< "B_y" <<setw(16)<< "B_z" <<setw(6)<< "trash\n";

    // radius for specific cut lines
    double radius_bend = R_1 + ((double)i_side-((double)Nr_shift-1.)/2.)*shift_side;
    // status report
    std::cout << "tangential line: " << i_up+i_side*3 << std::endl;
    // vertical component (drift direction)
    P[1] = ((double)i_up - ((double)Nr_shift - 1.) / 2.) * shift_up;
    // calculate B-field values at every step on the cut line:
    for(int j_steps = 0; j_steps < (steps_app+steps_RxB+steps_dec); j_steps++) {
      // moving along the geometric beam axis
      if (j_steps<steps_app) {                    // inside aperture area
        P[2] = radius_bend + R_1*(1-sin(M_PI/4.));
        P[3] = start_pos + step_size*(double)j_steps;
      } else if (j_steps<(steps_app+steps_RxB) && angle_step*((double)j_steps-steps_app) < M_PI/4.) { // inside RxB
        P[2] = (R_1*2. + ((double)i_side-((double)Nr_shift-1.)/2.)*shift_side)*cos(angle_step*((double)j_steps-steps_app))-R_1*sin(M_PI/4.);
        P[3] = bend_start + (R_1*2. + ((double)i_side-((double)Nr_shift-1.)/2.)*shift_side)*sin(angle_step*((double)j_steps-steps_app));
      } else if (j_steps<(steps_app+steps_RxB) && angle_step*((double)j_steps-steps_app) > 3.*M_PI/4.) { // inside RxB
        P[2] = (R_1*2. + ((double)i_side-((double)Nr_shift-1.)/2.)*shift_side)*cos(angle_step*((double)j_steps-steps_app))+R_1*sin(M_PI/4.);
        P[3] = bend_start + (R_1*2. + ((double)i_side-((double)Nr_shift-1.)/2.)*shift_side)*sin(angle_step*((double)j_steps-steps_app));
      } else if (j_steps<(steps_app+steps_RxB)) { // inside RxB
        P[2] = radius_bend*cos(angle_step*((double)j_steps-steps_app));
        P[3] = bend_start + radius_bend*sin(angle_step*((double)j_steps-steps_app))+R_1*cos(M_PI/4.);
      } else {                              // inside detector area
        P[2] = (radius_bend + R_1*(1-sin(M_PI/4.)))*cos(angle) + ((double)j_steps-(steps_app+steps_RxB))*step_size*sin(angle);
        P[3] = (radius_bend + R_1*(1-sin(M_PI/4.)))*sin(angle) + ((double)j_steps-(steps_app+steps_RxB))*step_size*cos(angle);
      }
      // calculating the B-field with Magfield 3
      magfield_elliptic(P,B,filedir);
      b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);  // claculating abs. B-field

      // store the datapoint's position
      fileout <<setiosflags(ios::left)<<setw(14)<< P[1] <<setw(14)<< P[2]<<setw(14)<< P[3]<<setw(8)<< "0";
      // store the B-field values in the datapoint
      fileout <<setiosflags(ios::left)<<setw(16)<< b <<setw(16)<< B[1] <<setw(16)<< B[2] <<setw(16)<< B[3] <<setw(8)<< "0" << "\n" ;
      }
    fileout.close();
    }
  }
}
