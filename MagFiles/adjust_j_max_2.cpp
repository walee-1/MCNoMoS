void adjust_j_max_2(double R_1, int N_coil, double d_coil, int N_dummy, string inputcoil, string filedir,int iterations){

// Magfield3 variables
double B[4];          // B-field components [T]
double P[4];          // Descartes coordinates [m]
double b;             // Absolute B-field [T]

// parameters
double b_goal = 0.2;  // goal magnetic field strength [T]
double b_goal_medial; // goal magnetic field strength at the beginning [T]
double b_goal_lateral;// goal magnetic field strength at the end of RxB [T]
double radial_shift = 0.02;   // shift of magnetic lines in radial direction [m]
double stepsize = 0.1;// stepsize of the integration
double j_old[100];           // old current setting of one coil [A/m²]
double j_new[100];           // new current setting of one coil [A/m²]
int i_detector = 6;     // until which coil the B-field should be uniform
int i_aperture = 0;     // number of coils in aperture in which the B-field should be uniform
int i_recoil = 0;       // number of first not adjusted recoil coil
int n_recoil = 0;       // number of coils on the other side of the pumpport
double B_buffer;        // buffer for the B-fiedl [T]
double NoMoS_start = -1.63307;  // beginning of NoMoS [T]
double  d_half = 0.05;  // half of l_coil
double fraction = 0.1;  // neigbor incluence fraction on the B-field
int n_steps = 5;           // number of points taking for averageing



// used variables for the file paths
    string input_old_path;          // old values
    string input_new_path;          // new values
// iterations
for (int it = 0; it < iterations; it++) {       // iterations to do
    // define filename of the old configuration
    if (it==0) {input_old_path=inputcoil;
    } else {input_old_path=input_new_path;}
    // define filename of the new configuration
    stringstream ssname;
    ssname << it+1;
    ssname >> input_new_path;
    input_new_path = filedir + "inputcoil_RxB_" + input_new_path + ".dat";

    // perpare Magfield
    input_coils(input_old_path);
    test_coils(filedir);
    magsource(filedir);

    //Writing a new Inputcoil.dat
    string buffer;
    double coil_number;
    ifstream inputcoil_old;
    inputcoil_old.open(input_old_path.c_str(),ios::in);
    ofstream inputcoil_new;
    inputcoil_new.open(input_new_path.c_str(),ios::out);
    inputcoil_old >> coil_number;
    inputcoil_new << coil_number << "\n";//(N_coil+2*N_dummy);

    // read in the old inputcoil file
    for (int i_coil=0; i_coil < coil_number;i_coil++) {
      inputcoil_old >> j_old[i_coil];
        for(int j_1=0; j_1 < 9; j_1++){
            inputcoil_old >> buffer;}
            //cout << i_coil +1 << "\t" << j_old[i_coil]<<endl;
      }

        inputcoil_old.close();

    // sweep through all coils and calculate the new j_old_levels
    for (int i_coil=0; i_coil < coil_number;i_coil++) {
            // RxB coils
      if (i_coil < N_coil && (double)i_coil*M_PI/((double)N_coil-1.) < 2.5) {
        P[1] = 0.;
        P[2] = R_1*cos((double)i_coil*M_PI/((double)N_coil-1.));
        P[3] = R_1*sin((double)i_coil*M_PI/((double)N_coil-1.));
        magfield_elliptic(P,B,filedir);                 // calculating the B-field
        b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);          // claculating abs. B-field
        if (i_coil>0) {j_new[i_coil-1]=j_new[i_coil-1]*(1.-fraction+fraction/b*b_goal);}      // change the previous coil
        j_new[i_coil]=j_old[i_coil]*(fraction*2.+(1.-fraction*2.)/b*b_goal);        // change the coil itseld
        j_old[i_coil+1]=j_old[i_coil+1]*(1.-fraction+fraction/b*b_goal);        // change the next coils
        if ( i_coil == 0) {j_new[i_coil]=j_old[i_coil]*(fraction*2.+(1.-fraction*2.)/b*b_goal*1.02); }       // change the coil itseld
        // find the values for the end and the beginning
        if (i_coil >= (N_coil - 1) / 2 && i_coil < (N_coil + 1) / 2) {
          P[2] = (R_1-radial_shift)*cos((double)i_coil*M_PI/((double)N_coil-1.));
          P[3] = (R_1-radial_shift)*sin((double)i_coil*M_PI/((double)N_coil-1.));
          magfield_elliptic(P,B,filedir);                   // calculating the B-field
          b_goal_medial=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);// claculating abs. B-field
          P[2] = (R_1+radial_shift)*cos((double)i_coil*M_PI/((double)N_coil-1.));
          P[3] = (R_1+radial_shift)*sin((double)i_coil*M_PI/((double)N_coil-1.));
          magfield_elliptic(P,B,filedir);                   // calculating the B-field
          b_goal_lateral=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);// claculating abs. B-field
          }
          // adjust the end coils to lower values to have a small gradient in the end
      } else if (i_coil < N_coil && (double)i_coil*M_PI/((double)N_coil-1.) >= 2.9) {
        P[1] = 0.;
        P[2] = (R_1+radial_shift)*cos((double)i_coil*M_PI/((double)N_coil-1.));
        P[3] = (R_1+radial_shift)*sin((double)i_coil*M_PI/((double)N_coil-1.));
        magfield_elliptic(P,B,filedir);                 // calculating the B-field
        b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);          // claculating abs. B-field
        j_new[i_coil-1]=j_new[i_coil-1]*(1.-fraction+fraction/b*b_goal_lateral);        // change the previous coil
        j_new[i_coil]=j_old[i_coil]*(fraction*2.+(1.-fraction*2.)/b*b_goal_lateral);    // change the coil itseld
        j_old[i_coil+1]=j_old[i_coil+1]*(1.-fraction+fraction/b*b_goal_lateral);        // change the next coils

            // Detector coils
      } else if (i_coil >= N_coil && i_coil < N_coil + i_detector) {    // homgenize the first coils in the detector
        b=0.;
        P[2] = -R_1;
        for (int i_half = 1; i_half < n_steps+1; i_half++) {
          P[3] = -(double)(i_coil-N_coil+1.)*d_coil+d_half/((double)n_steps-1.)*(double)i_half;
          magfield_elliptic(P,B,filedir);                 // calculating the B-field
          b+=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);}          // claculating abs. B-field
        b /= (double)n_steps;
        j_new[i_coil]=j_old[i_coil]/b*b_goal_lateral;
        B_buffer = b;
      } else if (i_coil >= N_coil+ i_detector && i_coil < N_coil + N_dummy) { // delete the maximas inside detector
        P[2] = -R_1;
        b=0.;
        for (int i_half = 1; i_half < n_steps+1; i_half++) {
          P[3] = -(double)(i_coil-N_coil+1.)*d_coil+d_half/((double)n_steps-1.)*(double)i_half;
          magfield_elliptic(P,B,filedir);                 // calculating the B-field
          b+=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);          // claculating abs. B-field
          }
        b /= (double)n_steps;
        cout << "b: " << b << "  B_buffer: "   << B_buffer << "  B_goal: " << b_goal_lateral << endl;
        if (b > B_buffer) {
            if (B_buffer < b_goal_lateral) {j_new[i_coil]=j_old[i_coil]/b*B_buffer;
            } else {j_new[i_coil]=j_old[i_coil]/b*b_goal_lateral;}
        }
        B_buffer = b;
        // adjust the last Aperture coils
        cout << "B_goal_medial: " << b_goal_medial << endl;
      } else if (i_coil >= N_coil+N_dummy && i_coil < N_coil+N_dummy + i_aperture) {
        P[2] = +R_1;
        b=0.;
        double B_buf[4];
        for (int i_half = 1; i_half < n_steps+1; i_half++) {
          P[3] = -(double)(i_coil-N_coil-N_dummy+1.)*d_coil+d_half/((double)n_steps-1.)*(double)i_half;
          magfield_elliptic(P,B,filedir);                 // calculating the B-field
          //b+=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);          // claculating abs. B-field
          B_buf[i_half] = sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);
          }
        for (int i_half = 1; i_half < n_steps+1; i_half++) {
          //if (b != 0 && b > B_buf[i_half]) b = B_buf[i_half];
          b= (B_buf[1]+B_buf[2])/2.;
        }
        j_new[i_coil]=j_old[i_coil]/b*b_goal_medial;
        //adjust the recoil coils
      } else if (i_coil >= coil_number - n_recoil && i_coil < coil_number -i_recoil) {
        P[2] = +R_1;
        b=0.;
        for (int i_half = 1; i_half < n_steps+1; i_half++){
          P[3] = NoMoS_start+(double)(coil_number-i_coil)*d_coil-d_half/((double)n_steps-1.)*(double)i_half;
          magfield_elliptic(P,B,filedir);                 // calculating the B-field
          b+=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);          // claculating abs. B-field
          }
        b /= (double)n_steps;
        j_new[i_coil]=j_old[i_coil]/b*b_goal;
        B_buffer = b;
      } else if (i_coil >= coil_number - i_recoil && i_coil < coil_number) {
        P[2] = +R_1;
        b=0.;
        for (int i_half = 1; i_half < n_steps+1; i_half++){
          P[3] = NoMoS_start+(double)(coil_number-i_coil)*d_coil-d_half/((double)n_steps-1.)*(double)i_half;
        //cout << P[3] << endl;
          magfield_elliptic(P,B,filedir);                 // calculating the B-field
          b+=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);          // claculating abs. B-field
          }
        b /= (double)n_steps;
        if (b > B_buffer) {
            if (B_buffer < b_goal) {j_new[i_coil]=j_old[i_coil]/b*B_buffer;
            } else {j_new[i_coil]=j_old[i_coil]/b*b_goal;}}
        B_buffer = b;
      } else {j_new[i_coil]=j_old[i_coil];}/*
        if (it==0){j_new[i_coil]=j_old/b*b_goal;}
        else {j_new[i_coil]=j_new[i_coil]/b*b_goal;}
        }
    for (int i=N_coil; i < N_coil+N_dummy; i++) {
        j_new[i]=j_old*j_lin;}
    for(int i=N_coil+N_dummy; i < N_coil+2*N_dummy; i++){
        j_new[i]= -j_old*j_lin;}
    /*for(int i=0; i < N_dummy;i++){                  // aperature area
        P[1] = 0.;
        P[2] = R_1;
        for(int j=0; j<= 1./stepsize; j++){
            P[3] = -d_coil*((double)i+(double)j*stepsize);
            magfield_elliptic(P,B,filedir);                 // calculating the B-field
            b = (b*(double)j + sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]))/((double)j+1.);          // claculating abs. B-field
            }
        if(it==0){j_new[N_coil+i]=j_old*j_lin/b*b_goal;}
        else{j_new[N_coil+i]=j_new[N_coil+i]/b*b_goal;}
        }
    for(int i=0; i < N_dummy;i++){                      // detector area
        P[1] = 0.;
        P[2] = -R_1;
        P[3] = -d_coil*((double)i+0.3);
        magfield_elliptic(P,B,filedir);                 // calculating the B-field
        b=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);          // claculating abs. B-field
        if(it==0){j_new[N_coil+N_dummy+i]=-j_old*j_lin/b*b_goal;}
        else{j_new[N_coil+N_dummy+i]=j_new[N_coil+N_dummy+i]/b*b_goal;}
        }*/
    }
    // Reopen the old inputfile for copying
    inputcoil_old.open(input_old_path.c_str(),ios::in);
    inputcoil_old >> buffer;    // read out coilnumber


 for (int i_coil=0; i_coil < coil_number;i_coil++) {
        // store the new currentvalues
        inputcoil_old >> buffer;
        if (fabs(j_new[i_coil]) <=5e6 && fabs(j_new[i_coil]) >=1e4 || i_coil > 70) {
            inputcoil_new <<setiosflags(ios::left)<<setw(14)<< j_new[i_coil];
        } else if (fabs(j_new[i_coil]) <=1e4) {
            inputcoil_new <<setiosflags(ios::left)<<setw(14)<< 1e4;
        } else {
            inputcoil_new <<setiosflags(ios::left)<<setw(14)<< 5e6;}
        for(int j_1=0; j_1 < 9; j_1++){ // copy all outer values
            inputcoil_old >> buffer;
            inputcoil_new <<setiosflags(ios::left)<<setw(12)<< buffer;
            }
        inputcoil_new << "\n";
        }
    /*inputcoil_old >> buffer;
    while(!inputcoil_old.eof()){
        //if(inputcoil_old.eof()){break;}
        inputcoil_new << j_new[N_coil];
        for(int j_1=0; j_1 < 9; j_1++){
            inputcoil_old >> buffer;
            inputcoil_new << buffer<< "\t";
            }
        inputcoil_new << "\n";
        inputcoil_old >> buffer;
        }*/
        inputcoil_new.close();
        inputcoil_old.close();
    }
}
