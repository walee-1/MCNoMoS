// different functions to creat the inputcoil.dat used for magfield3.cpp

//#include "geometrics.cpp"                // geometrics, rotations,...
void linear_coils_to_inputcoil_dat(double *coil_current, int i_coil, double *starting_point, double *direction, double coil_length, double r_coil, double coil_thickness, double coil_distance, int n_coils, string filename, int rev) {
    // write linear coils into the inputcoil.dat
/*  coil_parameters:
    coil_parameter[0] =
    coil_parameter[1] = starting_point[1]:          // x coordinate of the staring point [m]
    coil_parameter[2] = starting_point[2];          // y coordinate of the staring point [m]
    coil_parameter[3] = starting_point[3];          // z coordinate of the staring point [m]
    coil_parameter[4] = direction[1];               // x component of the coil axis direction [m]
    coil_parameter[5] = direction[2];               // y component of the coil axis direction [m]
    coil_parameter[6] = direction[3];               // z component of the coil axis direction [m]
    coil_parameter[7] = coil_length;                // length of the coil [m]
    coil_parameter[8] = coil_thickness;             // thickness of the coils [m]
    coil_parameter[9] = r_coil;                     // inner coil radius [m]
    coil_parameter[10] = coil_distance;             // distance between the coils [m]
    coil_parameter[11] = rev;                       // adding of the next coils before or after (1 or 0)
    coil_parameter[] =
*/

// normalize the direction:
normalize_vector(direction);

int j_coil;

// open the inputcoil.dat file to append new coils
ofstream inputcoil_file;
inputcoil_file.open(filename.c_str(), ios::app);

// write the coil by coil the data into the inputcoil.dat
for(int j_coil = 0; j_coil < n_coils; j_coil++){        // write more coils with the same aligment into the inputcoil.dat
    // 1. write the current density
    inputcoil_file <<setiosflags(ios::left)<<setw(14)<< coil_current[i_coil+j_coil];         // current density [A/m²]
    // 2-4 Descartes coordinates of starting point of the coil axis [m]
    int width_pos = 12;
    double coil_start = (coil_length+coil_distance)*(double)j_coil*pow(-1.,rev);
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< starting_point[1]+direction[1]*coil_start;   // x-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< starting_point[2]+direction[2]*coil_start;   // y-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< starting_point[3]+direction[3]*coil_start;   // z-coordinate of the coils beginning [m]
    // 5-7 Descartes coordinates of the  end point of the coil axis [m]
    double coil_end = coil_start + coil_length;
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< starting_point[1]+direction[1]*coil_end;     // x-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< starting_point[2]+direction[2]*coil_end;     // y-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< starting_point[3]+direction[3]*coil_end;     // z-coordinate of the coils beginning [m]
    // 8. write inner coil radius [m]
    int width_radius = 10;
    inputcoil_file <<setiosflags(ios::left)<<setw(width_radius) << r_coil;  // inner coil radius [m]
    // 9. write outer coil radius into inputcoil.dat [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_radius) << r_coil + coil_thickness; // outer coil radius [m]
    // 10. Precissions of magfield (=> 30 is proofen to be good) + End of line
    inputcoil_file <<setiosflags(ios::left)<< setw(4) <<  30 << "\n";       // Percission ?
    }

// close the inputcoil.dat file
inputcoil_file.close();
}



void circular_coils_to_inputcoil_dat(double *coil_current, int i_coil, double *starting_point, double *direction, double *grand_Radius, double r_coil, double coil_thickness, double coil_length, double alpha, int n_coils, string filename) {
// Parameters
int width_pos = 12;
int width_radius = 10;

// calculated parameters
double R_1 = vector_length(grand_Radius);           // absoult value of the grand_Radius R_1 of the RxB [m]
double centerpoint[4];          // center point of the coil circle
for (int i = 1; i <= 3; i++) {centerpoint[i] = starting_point[i] - grand_Radius[i];} // calculate the center point
//cout  << "centerpoint: ( " << centerpoint[1] << " , " << centerpoint[2] << " , " << centerpoint[3] << " )" << endl; // test if the centerpoint is right
ofstream inputcoil_file;                            // open inputcoil.dat
inputcoil_file.open(filename.c_str(), ios::app);

// write the coil by coil the data into the inputcoil.dat
for(int j_coil = 0; j_coil < n_coils; j_coil++){        // write more coils with the same alignment into the inputcoil.dat
    // 1. write the current density
    inputcoil_file <<setiosflags(ios::left)<<setw(14)<< coil_current[i_coil+j_coil];         // current density [A/m²]
    // 2-4 Descartes coordinates of starting point of the coil axis [m]
    double coil_angle = alpha/((double)n_coils-1.)*((double)j_coil);        // opening angle between to neighboring coils
    //cout << "coil number: " << j_coil + 1 << "  angle in degree:" << coil_angle/M_PI*180. << endl;      // test the calculation of the angle between first and current coil in RxB
    // calculate the coil center position (it is on the Radius rotated for the coil angle in the plan spanned by the Radius and the direction)
    double Rxd[4];
    crossproduct(grand_Radius,direction,Rxd);
    double new_R_1[4];
    for (int i = 1; i <= 3; i++) {new_R_1[i]= grand_Radius[i];}
    //rotate_vector_in_3d(grand_Radius,centerpoint,Rxd,coil_angle,new_R_1);             // not working because of signs now step by step rotation
    double pos_x = centerpoint[1] + new_R_1[1];
    double pos_y = centerpoint[2] + new_R_1[2];
    double pos_z = centerpoint[3] + new_R_1[3];
    //cout << "coil number " << j_coil+1 << " center pos: ( " << pos_x << " , " << pos_y << " , " << pos_z << " )" << endl; // test if the coil center was calculated right

    // write the position into the inpucoil.dat file
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_x - coil_length/2.*direction[1];   // x-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_y - coil_length/2.*direction[2];   // y-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_z - coil_length/2.*direction[3];   // z-coordinate of the coils beginning [m]
    // 5-7 Descartes coordinates of the  end point of the coil axis [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_x + coil_length/2.*direction[1];     // x-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_y + coil_length/2.*direction[2];     // y-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_z + coil_length/2.*direction[3];     // z-coordinate of the coils beginning [m]

    // 8. write inner coil radius [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_radius) << r_coil;  // inner coil radius [m]
    // 9. write outer coil radius into inputcoil.dat [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_radius) << r_coil + coil_thickness; // outer coil radius [m]
    // 10. Precissions of magfield (=> 30 is proofen to be good) + End of line
    inputcoil_file <<setiosflags(ios::left)<< setw(4) <<  30 << "\n";       // Percission ?

    // rotated the direction for the next step
    rotate_vector_in_3d(direction,centerpoint,Rxd,alpha/((double)n_coils-1.));
    rotate_vector_in_3d(grand_Radius,centerpoint,Rxd,alpha/((double)n_coils-1.));
    }
inputcoil_file.close();     // close inputcoil.dat
}



void circular_coils_to_inputcoil_dat(double *coil_current, int i_coil, double *starting_point, double *direction, double *grand_Radius, double r_coil, double coil_thickness, double coil_length, double alpha, int n_coils, string filename, double shift_angle) {
// Parameters
int width_pos = 12;
int width_radius = 10;

// calculated parameters
double R_1 = vector_length(grand_Radius);           // absoult value of the grand_Radius R_1 of the RxB [m]
double centerpoint[4];          // center point of the coil circle
for (int i = 1; i <= 3; i++) {centerpoint[i] = starting_point[i] - grand_Radius[i];} // calculate the center point
//cout  << "centerpoint: ( " << centerpoint[1] << " , " << centerpoint[2] << " , " << centerpoint[3] << " )" << endl; // test if the centerpoint is right
ofstream inputcoil_file;                            // open inputcoil.dat
inputcoil_file.open(filename.c_str(), ios::app);


// rotate the hole setup
//    rotate_vector_in_3d(direction,centerpoint,Rxd,shift_angle);
  //  rotate_vector_in_3d(grand_Radius,centerpoint,Rxd,shift_angle);

// write the coil by coil the data into the inputcoil.dat
for(int j_coil = 0; j_coil < n_coils; j_coil++){        // write more coils with the same alignment into the inputcoil.dat
    // 1. write the current density
    inputcoil_file <<setiosflags(ios::left)<<setw(14)<< coil_current[i_coil+j_coil];         // current density [A/m²]
    // 2-4 Descartes coordinates of starting point of the coil axis [m]
    double coil_angle = alpha/((double)n_coils-1.)*((double)j_coil);        // opening angle between to neighboring coils
    //cout << "coil number: " << j_coil + 1 << "  angle in degree:" << coil_angle/M_PI*180. << endl;      // test the calculation of the angle between first and current coil in RxB
    if (coil_angle > alpha) {break;}
    // calculate the coil center position (it is on the Radius rotated for the coil angle in the plan spanned by the Radius and the direction)
    double Rxd[4];
    crossproduct(grand_Radius,direction,Rxd);
    double new_R_1[4];
    for (int i = 1; i <= 3; i++) {new_R_1[i]= grand_Radius[i];}
    //rotate_vector_in_3d(grand_Radius,centerpoint,Rxd,coil_angle,new_R_1);             // not working because of signs now step by step rotation
    double pos_x = centerpoint[1] + new_R_1[1];
    double pos_y = centerpoint[2] + new_R_1[2];
    double pos_z = centerpoint[3] + new_R_1[3];
    //cout << "coil number " << j_coil+1 << " center pos: ( " << pos_x << " , " << pos_y << " , " << pos_z << " )" << endl; // test if the coil center was calculated right

    // write the position into the inpucoil.dat file
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_x - coil_length/2.*direction[1];   // x-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_y - coil_length/2.*direction[2];   // y-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_z - coil_length/2.*direction[3];   // z-coordinate of the coils beginning [m]
    // 5-7 Descartes coordinates of the  end point of the coil axis [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_x + coil_length/2.*direction[1];     // x-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_y + coil_length/2.*direction[2];     // y-coordinate of the coils beginning [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_pos)<< pos_z + coil_length/2.*direction[3];     // z-coordinate of the coils beginning [m]

    // 8. write inner coil radius [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_radius) << r_coil;  // inner coil radius [m]
    // 9. write outer coil radius into inputcoil.dat [m]
    inputcoil_file <<setiosflags(ios::left)<<setw(width_radius) << r_coil + coil_thickness; // outer coil radius [m]
    // 10. Precissions of magfield (=> 30 is proofen to be good) + End of line
    inputcoil_file <<setiosflags(ios::left)<< setw(4) <<  30 << "\n";       // Percission ?

    // rotated the direction for the next step
    rotate_vector_in_3d(direction,centerpoint,Rxd,alpha/((double)n_coils-1.));
    rotate_vector_in_3d(grand_Radius,centerpoint,Rxd,alpha/((double)n_coils-1.));
    }
inputcoil_file.close();     // close inputcoil.dat
}



void klothoid_coils_to_inputcoil_dat(double *coil_current, int i_coil, double *starting_point, double *direction, double *grand_Radius, double r_coil, double coil_thickness, double coil_distance, double alpha, int n_coils, string filename, double shift_angle) {

}



void PERC_coils_to_inputcoil_dat(double *starting_point, double *direction,double global_scale,double sol_scale,double filter_scale,double connec_scale, string filedir, string filename, int Perc_coordinatsystem){
    // Perc Parameterrs:
int Perc_N = 13;                         // number of coils in Perc
string Perc_text_name = filedir + "Perc_coils_2.txt";     // path to store the Perccoil configuration with the PERC's end in the origin and the PERC orientation
// Parameters of the PERC coils
  // solenoid
    double Perc_Solenoid_length = 8.000;          // length of the Solenoid [m]
    double Perc_Solenoid_radius = 0.200;              // inner radius of the solenoid [m]
    double Perc_Solenoid_thickness = 0.0098;         // coil thickness of the Solenoid [m]
    double Perc_Solenoid_j = 1067./Perc_Solenoid_length/Perc_Solenoid_thickness*8943.2*sol_scale*global_scale;  // maximum current density of the Solenoid and its correction coils to reach B0=0.5-1.5 T [A/m²]
    double Perc_Solenoid_center = -4.0;
    double Perc_Solenoid_cor_1_length = 0.1474;    // length of the Solenoid's first correction coil [m]
    double Perc_Solenoid_cor_1_radius = 0.2098;               // inner radius of the solenoid's first correction coil [m]
    double Perc_Solenoid_cor_1_thickness = 0.0049;  // coil thickness of the Solenoid's 1st correction coil [m]
    double Perc_Solenoid_cor_1_center = -7.9263;
    double Perc_Solenoid_cor_1_j = 1067./Perc_Solenoid_cor_1_length/Perc_Solenoid_cor_1_thickness*80.4*sol_scale*global_scale;
    double Perc_Solenoid_cor_2_length = 0.08;  // length of the Solenoid's 2nd correction coil [m]
    double Perc_Solenoid_cor_2_center = -0.1328;
    double Perc_Solenoid_cor_23_radius = 0.2073;   // inner radius of the solenoid's 2nd & 3rd correction coil [m]
    double Perc_Solenoid_cor_23_thickness = 0.0024;// coil thickness of the Solenoid's 2nd & 3rd correction coil [m]
    double Perc_Solenoid_cor_2_j = -1067./Perc_Solenoid_cor_2_length/Perc_Solenoid_cor_23_thickness*22.4*sol_scale*global_scale;   // current density of the Solenoid's 2nd &3rd correction coil to reach B0=0.5-1.5 T [A/m²]
    double Perc_Solenoid_cor_3_length = 0.023;  // length of the Solenoid's 3rd correction coil [m]
    double Perc_Solenoid_cor_3_center = -0.3144;
    double Perc_Solenoid_cor_3_j = -1067./Perc_Solenoid_cor_3_length/Perc_Solenoid_cor_23_thickness*6.4*sol_scale*global_scale;
  // Bending coils 1 - 3
    double Perc_Bender_1_length = 0.1723;          // length of the Bender_1 [m]
    double Perc_Bender_1_radius = 0.256;              // inner radius of the Bender_1 [m]
    double Perc_Bender_1_shift = -0.031;                // shift of the Bender_1-center  in the vertical axis away from the neutronbeam [m]
    double Perc_Bender_1_angle = 10.5/180.*M_PI;                  // angle of the first bendercoil [rad]
    double Perc_Bender_12_thickness = 0.0308;   // coil thickness of the Benders 1 & 2 [m]
    double Perc_Bender_1_j = 555./Perc_Bender_1_length/Perc_Bender_12_thickness*664.3*global_scale;    // current density of the Benders to reach B1=3-6 T [A/m²]
    double Perc_Bender_1_center = 0.1678;
    double Perc_Bender_2_length = 0.2456;          // length of the Bender_2 [m]
    double Perc_Bender_2_radius = 0.314;              // inner radius of the Bender_2 [m]
    double Perc_Bender_2_shift = 0.047;                 // shift of the Bender_2-center  in the vertical axis away from the neutronbeam[m]
    double Perc_Bender_2_angle = 21.5/180.*M_PI;                  // angle of the second bendercoil [rad]
    double Perc_Bender_2_j = 555./Perc_Bender_2_length/Perc_Bender_12_thickness*951.6*global_scale;
    double Perc_Bender_2_center = 0.4579;
    double Perc_Backbender_length = 0.2393;        // length of the Bender_2 [m]
    double Perc_Backbender_radius = 0.3118;            // inner radius of the Backbender [m]
    double Perc_Backbender_shift = 0.034;               // shift of the Backbender-center  in the vertical axis away from the neutronbeam[m]
    double Perc_Backbender_angle = -24./180.*M_PI;                 // angle of the backbendercoil [rad]
    double Perc_Backbender_thickness = 0.0356;  // coil thickness of the Backbender [m]
    double Perc_Backbender_j = 555./Perc_Backbender_length/Perc_Backbender_thickness*1069.5*global_scale;
    double Perc_Backbender_center = 1.9091;
  // Seperator
    double Perc_Selector_length = 0.8773;          // length of the Selector [m]
    double Perc_Selector_radius = 0.300;              // inner radius of the Selector [m]
    double Perc_Selector_thickness = 0.0356;     // coil thickness of the Selector [m]
    double Perc_Selector_shift = 0.08;                  // shift of the selectorcenter  in the vertical axis away from the neutronbeam [m]
    double Perc_Selector_j = 555./Perc_Selector_length/Perc_Selector_thickness*3961.5*global_scale;         // current density of the Selector to reach B1=3-6 T [A/m²]
    double Perc_Selector_center = 1.1787;
    double Perc_Selector_cor_12_length = 0.103;   // length of Selector's correction coils  [m]
    double Perc_Selector_cor_12_radius = 0.3368;       // inner radius of the Selector's correction coils [m]
    double Perc_Selector_cor_12_thickness = 0.0332;         // coil thickness of the Selector's first correction coil [m]
    double Perc_Selector_cor_12_shift = 0.08;
    double Perc_Selector_cor_12_j = 555./Perc_Selector_cor_12_length/Perc_Selector_cor_12_thickness*421.4*global_scale;
    double Perc_Selector_cor_1_center = 0.7916;
    double Perc_Selector_cor_2_center = 1.5659;
  // Filter Coils
    double Perc_Selector_Filter_12_length = 0.1914;// length of Selector's filter coils  [m]
    double Perc_Selector_Filter_12_radius = 0.3368;    // inner radius of the Selector's filter coils [m]
    double Perc_Selector_Filter_12_thickness = 0.0924;      // coil thickness of the Selector's filter coils [m]
    double Perc_Selector_Filter_12_j = 591./Perc_Selector_Filter_12_length/Perc_Selector_Filter_12_thickness*2215.2*filter_scale*global_scale;        // current density of the Selector's filter coils to reach B1=3-6 T [A/m²]
    double Perc_Selector_Filter_1_center = 0.98;
    double Perc_Selector_Filter_2_center = 1.3775;
    double Perc_Selector_Filter_12_shift = 0.08;
  // Detector/Connector Coils
    double Perc_Connector_length = 0.9965;          // length of the Connector [m]
    double Perc_Connector_radius = 0.250;              // inner radius of the Connector [m]
    double Perc_Connector_thickness = 0.0195;        // coil thickness of the Connector [m]
    double Perc_Connector_j = 560./Perc_Connector_length/Perc_Connector_thickness*2220.8*global_scale*connec_scale;   // current density of the Connector to reach B2=0,5-1 T [A/m²]
    double Perc_Connector_center = 2.6825;

// save the coil configuration in a text file
    ofstream Perc_text_file;
    Perc_text_file.open(Perc_text_name.c_str(), ios::out);

  // Solenoid
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Solenoid_j <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Solenoid_center - Perc_Solenoid_length/2. <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Solenoid_center+Perc_Solenoid_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Solenoid_radius <<setw(10)<< Perc_Solenoid_radius + Perc_Solenoid_thickness <<setw(4)<< 30 << endl;
  // Solenoid first correction coil
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Solenoid_cor_1_j <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Solenoid_cor_1_center -Perc_Solenoid_cor_1_length/2. <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Solenoid_cor_1_center +Perc_Solenoid_cor_1_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Solenoid_cor_1_radius <<setw(10)<< Perc_Solenoid_cor_1_radius + Perc_Solenoid_cor_1_thickness <<setw(4)<< 30 << endl;
  // Solenoid seconde correction coil
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Solenoid_cor_2_j <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Solenoid_cor_2_center - Perc_Solenoid_cor_2_length/2. <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Solenoid_cor_2_center + Perc_Solenoid_cor_2_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Solenoid_cor_23_radius <<setw(10)<< Perc_Solenoid_cor_23_radius + Perc_Solenoid_cor_23_thickness <<setw(4)<< 30 << endl;
  // Solenoid third correction coil
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Solenoid_cor_3_j <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Solenoid_cor_3_center - Perc_Solenoid_cor_3_length/2. <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Solenoid_cor_3_center + Perc_Solenoid_cor_3_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Solenoid_cor_23_radius <<setw(10)<< Perc_Solenoid_cor_23_radius + Perc_Solenoid_cor_23_thickness <<setw(4)<< 30 << endl;

  // Bender coil 1
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Bender_1_j <<setw(12)<< Perc_Bender_1_shift - Perc_Bender_1_length/2.*sin(Perc_Bender_1_angle) <<setw(12)<< 0. <<setw(12)<< Perc_Bender_1_center - Perc_Bender_1_length/2.*cos(Perc_Bender_1_angle);
    Perc_text_file <<setiosflags(ios::left)<<setw(12)<< Perc_Bender_1_shift + Perc_Bender_1_length/2.*sin(Perc_Bender_1_angle) <<setw(12)<< 0. <<setw(12)<< Perc_Bender_1_center + Perc_Bender_1_length/2.*cos(Perc_Bender_1_angle);
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Bender_1_radius <<setw(10)<< Perc_Bender_1_radius + Perc_Bender_12_thickness <<setw(4)<< 30 << endl;
  // Bender coil 2
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Bender_2_j <<setw(12)<< Perc_Bender_2_shift - Perc_Bender_2_length/2.*sin(Perc_Bender_2_angle) <<setw(12)<< 0. <<setw(12)<< Perc_Bender_2_center - Perc_Bender_2_length/2.*cos(Perc_Bender_2_angle);
    Perc_text_file <<setiosflags(ios::left)<<setw(12)<< Perc_Bender_2_shift + Perc_Bender_2_length/2.*sin(Perc_Bender_2_angle) <<setw(12)<< 0. <<setw(12)<< Perc_Bender_2_center + Perc_Bender_2_length/2.*cos(Perc_Bender_2_angle);
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Bender_2_radius <<setw(10)<< Perc_Bender_2_radius + Perc_Bender_12_thickness <<setw(4)<< 30 << endl;

  // Selector coil
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Selector_j <<setw(12)<< Perc_Selector_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_center- Perc_Selector_length/2. <<setw(12)<< Perc_Selector_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_center+Perc_Selector_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Selector_radius <<setw(10)<< Perc_Selector_radius + Perc_Selector_thickness <<setw(4)<< 30 << endl;
  // Selector first correction coil
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Selector_cor_12_j <<setw(12)<< Perc_Selector_cor_12_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_cor_1_center - Perc_Selector_cor_12_length/2. <<setw(12)<< Perc_Selector_cor_12_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_cor_1_center+Perc_Selector_cor_12_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Selector_cor_12_radius <<setw(10)<< Perc_Selector_cor_12_radius + Perc_Selector_cor_12_thickness <<setw(4)<< 30 << endl;
  // Selector second correction coil
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Selector_cor_12_j <<setw(12)<< Perc_Selector_cor_12_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_cor_2_center - Perc_Selector_cor_12_length/2. <<setw(12)<< Perc_Selector_cor_12_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_cor_2_center+Perc_Selector_cor_12_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Selector_cor_12_radius <<setw(10)<< Perc_Selector_cor_12_radius + Perc_Selector_cor_12_thickness <<setw(4)<< 30 << endl;

  // first filter coil
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Selector_Filter_12_j <<setw(12)<< Perc_Selector_Filter_12_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_Filter_1_center - Perc_Selector_Filter_12_length/2. <<setw(12)<< Perc_Selector_Filter_12_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_Filter_1_center+Perc_Selector_Filter_12_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Selector_Filter_12_radius <<setw(10)<< Perc_Selector_Filter_12_radius + Perc_Selector_Filter_12_thickness <<setw(4)<< 30 << endl;
  // second filter coil
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Selector_Filter_12_j <<setw(12)<< Perc_Selector_Filter_12_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_Filter_2_center - Perc_Selector_Filter_12_length/2. <<setw(12)<< Perc_Selector_Filter_12_shift <<setw(12)<< 0. <<setw(12)<< Perc_Selector_Filter_2_center+Perc_Selector_Filter_12_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Selector_Filter_12_radius <<setw(10)<< Perc_Selector_Filter_12_radius + Perc_Selector_Filter_12_thickness <<setw(4)<< 30 << endl;

  //Backbender
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Backbender_j <<setw(12) << Perc_Backbender_shift - Perc_Backbender_length/2.*sin(Perc_Backbender_angle) <<setw(12)<< 0. <<setw(12)<< Perc_Backbender_center- Perc_Backbender_length/2.*cos(Perc_Backbender_angle);
    Perc_text_file <<setiosflags(ios::left)<<setw(12)<< Perc_Backbender_shift + Perc_Backbender_length/2.*sin(Perc_Backbender_angle) <<setw(12)<< 0. <<setw(12)<< Perc_Backbender_center + Perc_Backbender_length/2.*cos(Perc_Backbender_angle);
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Backbender_radius <<setw(10)<< Perc_Backbender_radius + Perc_Backbender_thickness <<setw(4)<< 30 << endl;

  // Connector
    Perc_text_file <<setiosflags(ios::left)<<setw(16)<< Perc_Connector_j <<setw(12) << 0. <<setw(12)<< 0. <<setw(12)<< Perc_Connector_center- Perc_Connector_length/2. <<setw(12)<< 0. <<setw(12)<< 0. <<setw(12)<< Perc_Connector_center+Perc_Connector_length/2.;
    Perc_text_file <<setiosflags(ios::left)<<setw(10)<< Perc_Connector_radius <<setw(10)<< Perc_Connector_radius + Perc_Connector_thickness <<setw(4)<< 30 << endl;

  Perc_text_file.close();

  // store the Perc-coils into the inputcoil.dat file
double PercEndinPercSystem = Perc_Connector_center+Perc_Connector_length/2.;
ifstream Perc_coils;
Perc_coils.open(Perc_text_name.c_str(), ios::in);   // open the File with the information of all PERC-coils
double buffer_p, buffer_xp, buffer_yp, buffer_zp;   // buffer for read out
ofstream inputcoil_file;
inputcoil_file.open(filename.c_str(),ios::app);     // open the inputcoil.dat file
for (int Perc_i = 0; Perc_i <Perc_N; Perc_i++) {    // copy and relocate coil by coil
    // 1 copy the current density [a/m²]
    Perc_coils >> buffer_p;
    inputcoil_file <<setiosflags(ios::left)<<setw(14)<< buffer_p;   // current density [A/m²]
    // 2-7 copy the starting&end point of the coils [m]
    for (int i = 0; i < 2; i++){
        Perc_coils >> buffer_xp >> buffer_yp >> buffer_zp;              // load the Descartes coordinates of beginn of the coil [m]
        if (Perc_coordinatsystem == 1) {          // store in PERC coordinate system: xp=>xp, yp=yp, zp=>zp (x upwards orientation)
            inputcoil_file <<setiosflags(ios::left)<<setw(12)<< buffer_xp+starting_point[1] <<setw(12)<< buffer_yp+starting_point[2] <<setw(12)<< buffer_zp -PercEndinPercSystem +starting_point[3];
        } else if (Perc_coordinatsystem == 0) {          // Normal c0ordinate system: xp=>z, yp=-y, zp=>x (z upwards orientation)
            inputcoil_file <<setiosflags(ios::left)<<setw(12)<< buffer_zp+starting_point[1] <<setw(12)<< -buffer_yp+starting_point[2] <<setw(12)<< buffer_xp+starting_point[3];
        } else if (Perc_coordinatsystem == 2) {          // PERC cordinate system: xp=>-yp, yp=xp, zp=>zp but lying PERC (-y upwards orientation)
            inputcoil_file <<setiosflags(ios::left)<<setw(12)<< buffer_yp+starting_point[1] <<setw(12)<< -buffer_xp+starting_point[2] <<setw(12)<< buffer_zp -PercEndinPercSystem +starting_point[3];
        }
      }
     // 8-10 copy radius and precision
     for (int i = 0; i < 3; i++){
        Perc_coils >> buffer_p;
        inputcoil_file <<setiosflags(ios::left)<<setw(10)<< buffer_p;
        }
        inputcoil_file << "\n";     // End the line
    }
    inputcoil_file.close();
}
