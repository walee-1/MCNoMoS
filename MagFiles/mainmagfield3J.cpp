#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <string.h>
#include <iostream>                     // for c streams (cout <<)
#include <fstream>                      // ofr filstreams (input&output)
#include <sstream>                      // for stringstream
#include <iomanip>                      // for setw()
#include <string>                       // using strings
using namespace std;                    // makes std::cout , std::string unnecessary


// for magfield3
#define Ncoilmax 100
#define nmaxmag 4000
#define Ncenmax 5000
#define Ncenaxisymmmax 5000

//other files
#include "magfield3.cpp"                // Ference's programm already adapted for C++
#include "geometrics.cpp"                // geometrics, rotations,...

#include "input_coil_creat.cpp"         // write the coil parameters into the inputcoil.dat for different configuration (linear, circular, clothoide,...
#include "adjust_currents.cpp"          // adjust the currents in many different ways
#include "adjust_j_max_2.cpp"           // old adjustment function

#include "b_line.cpp"                   // calculating the B-fields lines starting in the aperture (foreward and backward)
#include "b_line_realdrift.cpp"
#include "coil_center.cpp"              // calculates the B-field following the expected beam size following the coils central axis
#include "radial_cuts.cpp"              // calculates the B-field on radial lines
#include "vertical_cuts.cpp"            // calculates the B-field on vertical lines
#include "fieldmap.cpp"                 // calculates the B-field on the central plane
#include "coil_center_2.cpp"
#include "geoline+anadrift.cpp"
#include "cooling.cpp"
#include "elecPower.cpp"

#include "../ConfigFiles/config.h"


//main braucht envp für parser, weil da error logs weiter gegeben werden
int main(int argc,char** argv, char* envp[] )
{/* mainmagfield3.cpp:
	    Global Parameters
	        General Settings for the simulation
	    NoMoS Parameter
	        All parameters needed for the simulation (Geometrie,...)
	    Sweep adjustments
	        Possibility to make several simulations with varied parameters.
	    Write the inputcoil.dat file
	        first clothoid
	        RxB
	        double layer
	        detector
	        Aperture
	        Recoil
	        Adapter
	        PERC
	    Adjust current densities
	    Simulations
	
	*/

	//directory from program call
	if (argc <2) {cout << "MAG: no directory given!" << endl; return 1;}
	string filedir =argv[1];
	cout <<"MAG: dir given: " <<  filedir << endl;


/////////////Config parser////////////// actually, name of config file should be given in bash file!
	if(argc <3) {cout << "MAG: no config file given!" << endl; return 1;}
	string configname = argv[2];
	cout << "MAG: config given: " << configname << endl;

	Config myconfig(configname,envp);// dir, where executable will be located, so without "../"

//////////////////////////////////////////////////////////////////////////////////
//-------------------------Parameter setting
//////////////////////////////////////////////////////////////////////////////////

	int Perc_coordinatsystem = 2;  // = 0 normal descard system: beam direction x, transversal horizontal y, vertical z
	                            // = 1 PERC cordinate system: x=>z, y=-y, z=>x
	                            // = 2 PERC coordinate system but lying PERC (reality standing NoMoS)
	int Wang_on = 0;            // if = 1 then Wang's settings are applied
	int doublecoil = 0;         // if = 1: additional coils in the RxB inbetween the normal coils
	bool coilcenterbool= myconfig.pBool("coilcenter");
	bool geodriftbool = myconfig.pBool("geo+drift");
	bool blinebool = myconfig.pBool("bline");
	bool fieldmapbool = myconfig.pBool("fieldmap");
	double apertYshift = myconfig.pDouble("apertYshift");
	double apertXshift = myconfig.pDouble("apertXshift");
	int horilines = myconfig.pInt("horilines");
	int vertilines = myconfig.pInt("vertilines");
	double ApertX = myconfig.pDouble("ApertX");
	double ApertY = myconfig.pDouble("ApertY");
	
	int Perc_on = myconfig.pBool("PercOn");
	double perc_global_scale = myconfig.pDouble("perc_global_scale");
	double perc_sol_scale = myconfig.pDouble("perc_sol_scale");
	double perc_filter_scale = myconfig.pDouble("perc_filter_scale");
	int nomos_on = myconfig.pBool("NoMoSOn");
	string Conductor = myconfig.pString("Conductor");
	double raisin = myconfig.pDouble("Raisin");
	double water_dia = myconfig.pDouble("water_dia");
	double i_max[100];                                          // real current densities for all coils.
	double PERClastcoil_to_flansch = myconfig.pDouble("PERClastcoil_to_flansch");
	double flansch_perc_thick = myconfig.pDouble("flansch_perc_thick");
	double fieldmap_step = myconfig.pDouble("fieldmap_step");

	// RxB Parameter
	double RxB_conduct = myconfig.pDouble("RxBconduct");
	int RxB_l = myconfig.pInt("RxBlTurns");
	int RxB_b = myconfig.pInt("RxBbTurns");
	double d_space = myconfig.pDouble("RxBgap");//simulation_parameters[5];       // minimal spacing between coils inside the RxB [m]
	int N_coil = myconfig.pInt("NCoils");//simulation_parameters[2];                 // number of coils in the RxB
	double alpha = myconfig.pDouble("alpha")/180.*M_PI;//simulation_parameters[6];   // opening angle of the RxB [rad]
	double RxBCurrent = myconfig.pDouble("RxBCurrent");
	double RxB_second_scale = myconfig.pDouble("RxB_second_scale");
	double l_coil_RxB =  RxB_conduct*RxB_l+raisin*(RxB_l-1);//length of RxB coil
	double t_coil_RxB = RxB_b*RxB_conduct+(RxB_b-1)*raisin;//simulation_parameters[4];      // thickness of the RxB coils [m]
	double RxBCurrDens = RxBCurrent*RxB_l*RxB_b/(l_coil_RxB * t_coil_RxB);//cross section of conductors includes raisin inside, not around [A/m²]
	// thickness calculated is without outer raisin, because thats how coils are defined. for curvature Radius, we include outer raisin aswell
	double r_coil = myconfig.pDouble("RxBInner");//simulation_parameters[3];      // inner coil radius of all coils [m]
	double R_1 = r_coil + t_coil_RxB + 2.*raisin + ( (l_coil_RxB + 2*raisin)/2. + d_space ) / tan ( alpha/(2*(N_coil -1))  );
	double FirstLastRxB_scale = myconfig.pDouble("FirstLastRxB_scale");

	// Outlet Coils (Detector)
	int outlet_l = myconfig.pInt("outlet_l");//simulation_parameters[18];            // linear coil length [m]
	int outlet_b = myconfig.pInt("outlet_b");//r2_wire;//simulation_parameters[19]; // linear coil thickness [m]
	double d_outlet = myconfig.pDouble("d_outlet");//0.005;//simulation_parameters[20];   // distance between linear coils [m]
	int n_outlet = myconfig.pInt("n_outlet");                // number of coils in the detector region (without correction coils)
	double l_coil_outlet = outlet_l*RxB_conduct + (outlet_l -1)*raisin;
	double t_coil_outlet = outlet_b*RxB_conduct + (outlet_b -1)*raisin;
	double outlet_r_shift = myconfig.pDouble("outlet_r_shift");             // s_coil_dec radial shift of all detector coils
	double outlet_x_shift = myconfig.pDouble("outlet_x_shift");
	double outlet_scale = myconfig.pDouble("outlet_curr_scale");
	double outlet_inner = myconfig.pDouble("outlet_inner");

	// Outlet Correction
	int n_outletCorr = myconfig.pInt("n_outletCorr");
	int outletCorr_l;
	int outletCorr_b;
	double d_outletCorr;
	double l_coil_outletCorr;
	double t_coil_outletCorr;
	double outletCorr_r_shift;
	double outletCorr_x_shift;
	double outletCorr_scale;
	double outletCorr_inner_gap;
	double outletCorr_inner;
	double outletCorr_fromoutEnd;
	if( n_outletCorr > 0 ){
		outletCorr_l = myconfig.pInt("outletCorr_l");
		outletCorr_b = myconfig.pInt("outletCorr_b");
		d_outletCorr = myconfig.pDouble("d_outletCorr");
		l_coil_outletCorr = outletCorr_l*RxB_conduct + (outletCorr_l -1)*raisin;
		t_coil_outletCorr = outletCorr_b*RxB_conduct + (outletCorr_b -1)*raisin;
		outletCorr_r_shift = myconfig.pDouble("outletCorr_r_shift");             // s_coil_dec radial shift of all detector coils
		outletCorr_x_shift = myconfig.pDouble("outletCorr_x_shift");
		outletCorr_scale = myconfig.pDouble("outletCorr_curr_scale");
		outletCorr_inner_gap = myconfig.pDouble("outletCorr_inner_gap"); // gap between outer radius of underlying coil and helmholtz coil
		outletCorr_inner = outlet_inner + t_coil_outlet + outletCorr_inner_gap;
		outletCorr_fromoutEnd = myconfig.pDouble("outletCorr_fromoutEnd");
	}

	// INLET	
	// for all Inlet coils
	double inlet_r_shift = myconfig.pDouble("inlet_r_shift");           // shift of all coils in the aperture region in the y direction
	double inlet_x_shift = myconfig.pDouble("inlet_x_shift");
	double inlet_inner = myconfig.pDouble("inlet_inner");

	// AF Inlet Coils
	int af_inlet_l = myconfig.pInt("af_inlet_l");
	int af_inlet_b = myconfig.pInt("af_inlet_b");
	double d_af_inlet = myconfig.pDouble("d_af_inlet");
	int n_af_inlet = myconfig.pInt("n_af_inlet");
	double l_coil_af_inlet = af_inlet_l*RxB_conduct + (af_inlet_l -1)*raisin;
	double t_coil_af_inlet = af_inlet_b*RxB_conduct + (af_inlet_b -1)*raisin;
	double af_inlet_scale = myconfig.pDouble("af_inlet_curr_scale");
	double af_inlet_inner = inlet_inner;
	
	// BF Inlet Coils
	int bf_inlet_l = myconfig.pInt("bf_inlet_l");
	int bf_inlet_b = myconfig.pInt("bf_inlet_b");
	double d_bf_inlet = myconfig.pDouble("d_bf_inlet");
	int n_bf_inlet = myconfig.pInt("n_bf_inlet");
	double l_coil_bf_inlet = bf_inlet_l*RxB_conduct + (bf_inlet_l -1)*raisin;
	double t_coil_bf_inlet = bf_inlet_b*RxB_conduct + (bf_inlet_b -1)*raisin;
	double bf_inlet_scale = myconfig.pDouble("bf_inlet_curr_scale");
	double bf_inlet_inner = inlet_inner;

	double inlet_dv_scale = myconfig.pDouble("inlet_dv_coil_scale");
	int inlet_dv_coil_start = myconfig.pInt("inlet_dv_coil_start");
	int inlet_dv_offN = myconfig.pInt("inlet_dv_offN");
	
	// filter Inlet coil
	int filter_inlet_l = myconfig.pInt("filter_inlet_l");
	int filter_inlet_b = myconfig.pInt("filter_inlet_b");
	double filter_inlet_scale = myconfig.pDouble("filter_inlet_curr_scale");
	double l_coil_filter_inlet = filter_inlet_l*RxB_conduct + (filter_inlet_l -1)*raisin;
	double t_coil_filter_inlet = filter_inlet_b*RxB_conduct + (filter_inlet_b -1)*raisin;
	double filter_inlet_inner = inlet_inner;

	// INNER Inlet coil
	int inner_inlet_l = myconfig.pInt("inner_inlet_l");
	int inner_inlet_b = myconfig.pInt("inner_inlet_b");
	double l_coil_inner_inlet = inner_inlet_l*RxB_conduct + (inner_inlet_l -1)*raisin;
	double t_coil_inner_inlet = inner_inlet_b*RxB_conduct + (inner_inlet_b -1)*raisin;
	double inner_inlet_r_shift = myconfig.pDouble("inner_inlet_r_shift");           // shift of all coils in the aperture region in the y direction
	double inner_inlet_x_shift = myconfig.pDouble("inner_inlet_x_shift");
	double inner_inlet_scale = myconfig.pDouble("inner_inlet_curr_scale");
	double inner_inlet_inner = myconfig.pDouble("inner_inlet_inner");
	if(inner_inlet_inner + t_coil_inner_inlet >= inlet_inner && inner_inlet_scale != 0.){
		cout << "MAG: INNER Inlet coil doesn't fit" << endl;
		return 1;
	}

	// filter coils parameters
	int n_filter = myconfig.pInt("n_filter");               // number of filter coils
	int filter_l = myconfig.pInt("filter_l");
	int filter_b = myconfig.pInt("filter_b");
	double d_filter = myconfig.pDouble("d_filter");
	double filter_inner_gap = myconfig.pDouble("filter_inner_gap"); // gap between outer radius of underlying coil and helmholtz coil
	double l_coil_filter = filter_l*RxB_conduct + (filter_l -1)*raisin;
	double t_coil_filter = filter_b*RxB_conduct + (filter_b -1)*raisin;
	double filter_r_shift = myconfig.pDouble("filter_r_shift");           // shift of all coils in the aperture region in the y direction
	double filter_x_shift = myconfig.pDouble("filter_x_shift");
	double filter_scale = myconfig.pDouble("filter_curr_scale");
	double filter_inner = inlet_inner + t_coil_filter_inlet + filter_inner_gap;
	
	// connector coils
	// after gate
	int n_af_connec = myconfig.pInt("af_connec_n");    
	int af_connec_l = myconfig.pInt("af_connec_l");
	int af_connec_b = myconfig.pInt("af_connec_b");
	double d_af_connec = myconfig.pDouble("d_af_connec");
	double l_coil_af_connec = af_connec_l*RxB_conduct + (af_connec_l -1)*raisin;
	double t_coil_af_connec = af_connec_b*RxB_conduct + (af_connec_b -1)*raisin;
	double af_connec_r_shift = myconfig.pDouble("af_connec_r_shift");           // shift of all coils in the aperture region in the y direction
	double af_connec_x_shift = myconfig.pDouble("af_connec_x_shift");
	double af_connec_scale = myconfig.pDouble("af_connec_curr_scale");
	double af_connec_inner = myconfig.pDouble("af_connec_inner");
	if(af_connec_inner != inlet_inner ) {cout << "MAG: InnerR of af_connec not equal to Inlet InnerR" << endl;}
	double af_connec_savedist = myconfig.pDouble("af_connec_savedist");
	
	// before gate
	int n_bf_connec = myconfig.pInt("bf_connec_n");    
	int bf_connec_l = myconfig.pInt("bf_connec_l");
	int bf_connec_b = myconfig.pInt("bf_connec_b");
	double d_bf_connec = myconfig.pDouble("d_bf_connec");
	double l_coil_bf_connec = bf_connec_l*RxB_conduct + (bf_connec_l -1)*raisin;
	double t_coil_bf_connec = bf_connec_b*RxB_conduct + (bf_connec_b -1)*raisin;
	double bf_connec_r_shift = myconfig.pDouble("bf_connec_r_shift");           // shift of all coils in the aperture region in the y direction
	double bf_connec_x_shift = myconfig.pDouble("bf_connec_x_shift");
	double bf_connec_scale = myconfig.pDouble("bf_connec_curr_scale");
	double bf_connec_inner = myconfig.pDouble("bf_connec_inner");
	if(bf_connec_inner != inlet_inner ) {
		cout << "MAG: InnerR of bf_connec not equal to af InnerR" << endl;
		if(bf_connec_x_shift != 0.) {cout << "MAG: bf_connec coils go around PERC Flansch. On purpose?" << endl;}
	}
	double bf_connec_savedist = myconfig.pDouble("bf_connec_savedist");
	double gate_dist = myconfig.pDouble("gate_dist");
	
	// Helmholtz coils around Pumpport
	int helm_l = myconfig.pInt("helm_l");
	int helm_b = myconfig.pInt("helm_b");
	double helm_inner_gap = myconfig.pDouble("helm_inner_gap"); // gap between outer radius of underlying coil and helmholtz coil
	double l_coil_helm = helm_l*RxB_conduct + (helm_l -1)*raisin;
	double t_coil_helm = helm_b*RxB_conduct + (helm_b -1)*raisin;
	double helm_r_shift = myconfig.pDouble("helm_r_shift");           // shift of all coils in the aperture region in the y direction
	double helm_x_shift = myconfig.pDouble("helm_x_shift");
	double pumpport_distance = myconfig.pDouble("pumpport_distance");//simulation_parameters[29];    // Distance between adapter coils and aperture coils of NoMoS [m]
	double helm_scale = myconfig.pDouble("helm_curr_scale");
	double helm_inner;
	char gateinfo;
	if(inlet_inner + t_coil_bf_inlet < af_connec_inner + t_coil_af_connec) { 
	      helm_inner = af_connec_inner + t_coil_af_connec + helm_inner_gap;
	}
	else{helm_inner = inlet_inner + t_coil_bf_inlet + helm_inner_gap;}	
	if(l_coil_helm > l_coil_af_connec){
		cout << "MAG: l_helm > l_af_connec, abort? (else, l_coil_helm gets overwritten by af_connec)" << endl;
		cin >> gateinfo;
		if(gateinfo == 'y'){return 1;}
		else{l_coil_helm = l_coil_af_connec;}
	}


	// Perc Parameterrs:
	int Perc_N = 13;                // number of coils of the PERC system
	double perclastcoil_toflansch = myconfig.pDouble("PERClastcoil_to_flansch");
	


	
	///////////////////////////////////////////////////////////////
////////////////////////////// Supra leiter coil parameters /////////////////////
	////////////////////////////////////////////////////////////////// 
	
	double SL_length, SL_thick;
	double SL_innerR = 0.37/2. + 0.01;
	if ( Conductor == "HTSL" ){
		SL_length = 0.0045;
		SL_thick = 0.00065;
		raisin = 0.; // HTSL sizes are already with isolation, therefore raisin = 0
	}
	else if( Conductor == "TTSL" ){
		SL_length = 0.0005;
		SL_thick = 0.0005;
		raisin = 0.0019 - SL_length;
	}

	double RxBtoOut_d = d_outlet;
	double RxBtoIn_d = d_af_inlet;
	double FIntoAFIn_d = d_af_inlet;
	double FIntoBFIn_d = d_bf_inlet;
	double HelmToBF_off = 0.;
	double HelmToConnec = 0.;

	if(Conductor == "HTSL" || Conductor == "TTSL" ){   // here, we have to adapt the coil data to super conductor values


		// RxB coils 
		l_coil_RxB =  SL_length*RxB_l + (RxB_l - 1)*raisin;
		t_coil_RxB = SL_thick*RxB_b + (RxB_b - 1)*raisin;
		r_coil = SL_innerR;
		RxBCurrDens = RxBCurrent*RxB_l*RxB_b/(l_coil_RxB * t_coil_RxB); // same current * turns but on smaller area -> much higher density
		R_1 = r_coil + t_coil_RxB + 2.*raisin + ( (l_coil_RxB + 2*raisin)/2. + d_space ) / tan ( alpha/(2*(N_coil -1))  );


		// Outlet
		l_coil_outlet = SL_length*outlet_l + (outlet_l -1)*raisin;
		t_coil_outlet = SL_thick*outlet_b + (outlet_b -1)*raisin;
		outlet_inner = SL_innerR + 0.05;
		
		// Outlet Corr
		l_coil_outletCorr = SL_length*outletCorr_l + (outletCorr_l -1)*raisin;
		t_coil_outletCorr = SL_thick*outletCorr_b + (outletCorr_b -1)*raisin;
		outletCorr_inner = outlet_inner + t_coil_outlet + outletCorr_inner_gap; // with new t, get new inner R


		// AF Inlet
		l_coil_af_inlet = SL_length*af_inlet_l + (af_inlet_l -1)*raisin;
		t_coil_af_inlet = SL_thick*af_inlet_b + (af_inlet_b -1)*raisin;
		af_inlet_inner = SL_innerR;

		// BF Inlet
		l_coil_bf_inlet = SL_length*bf_inlet_l +(bf_inlet_l -1)*raisin;
		t_coil_bf_inlet = SL_thick*bf_inlet_b +(bf_inlet_b -1)*raisin;
		bf_inlet_inner = SL_innerR;

		// filter Inlet
		l_coil_filter_inlet = SL_length*filter_inlet_l +(filter_inlet_l -1)*raisin;
		t_coil_filter_inlet = SL_thick*filter_inlet_b +(filter_inlet_b -1)*raisin;
		filter_inlet_inner = SL_innerR;

		/*
		// inner Inlet
		average_radius = inner_inlet_inner + t_coil_inner_inlet/2.;
		Delta_l = SL_length*inner_inlet_l; // we use delta l inbetween so we can calc delta l before changing l_coil
		Delta_l = l_coil_inner_inlet - Delta_l;
		l_coil_inner_inlet = SL_length*inner_inlet_l;
		t_coil_inner_inlet = SL_thick*inner_inlet_b;
		inner_inlet_inner = average_radius - t_coil_inner_inlet/2.; // with new t, get new inner R
		*/

		// filter
		l_coil_filter = SL_length*filter_l +(filter_l -1)*raisin;
		t_coil_filter = SL_thick*filter_b +(filter_b -1)*raisin;
		filter_inner = filter_inlet_inner + t_coil_filter_inlet + filter_inner_gap;


		// AF Connector
		l_coil_af_connec = SL_length*af_connec_l +(af_connec_l -1)*raisin;
		t_coil_af_connec = SL_thick*af_connec_b +(af_connec_b -1)*raisin;
		af_connec_inner = SL_innerR;

		// BF Connector
		l_coil_bf_connec = SL_length*bf_connec_l +(bf_connec_l -1)*raisin;
		t_coil_bf_connec = SL_thick*bf_connec_b +(bf_connec_b -1)*raisin;
		bf_connec_inner = SL_innerR;

		// helmholtz
		l_coil_helm = SL_length*helm_l +(helm_l -1)*raisin;
		t_coil_helm = SL_thick*helm_b +(helm_b -1)*raisin;
		if( bf_inlet_inner + t_coil_bf_inlet > af_connec_inner + t_coil_af_connec ){
			helm_inner = bf_inlet_inner + t_coil_bf_inlet + helm_inner_gap;
		}
		else{
			helm_inner = af_connec_inner + t_coil_af_connec + helm_inner_gap;
		}


	}

	cout << "MAG: R_1 =" << R_1 << endl; 
	cout << "MAG: RxBCurrDens = " << RxBCurrDens << endl;	


//////////////////////////////	// current set /////////////////////////////
	
	int counter = 0;
	
	// RxB
	for (int i = 0; i < N_coil; i++) {
		if(i == 0 || i == N_coil -1){ // First or Last RxB Coil
			i_max[i] = RxBCurrDens * FirstLastRxB_scale;
		}
		else if( i >= (N_coil-1)/2 ){
			i_max[i] = RxBCurrDens * RxB_second_scale;
		}
		else{
			i_max[i] = RxBCurrDens;
		}
	}	
	
	// Outlet
	for (int i = N_coil; i < N_coil + n_outlet; i++) {i_max[i] = RxBCurrDens*outlet_scale;}
	
	// Inlet
	// AF inlet
	for (int i = N_coil + n_outlet; i < N_coil + n_outlet + n_af_inlet; i++) {
		i_max[i] = RxBCurrDens*af_inlet_scale;
		counter ++;
	}	
	// Inlet Filter coil
	i_max[ N_coil + n_outlet + n_af_inlet ] = RxBCurrDens * filter_inlet_scale;
	counter ++;
	// BF inlet
	for (int i = N_coil + n_outlet + n_af_inlet + 1; i < N_coil + n_outlet + n_af_inlet + 1 + n_bf_inlet; i++) {
		if( counter >= n_af_inlet+1+n_bf_inlet-1 - (inlet_dv_coil_start-1) - (inlet_dv_offN-1) && counter <= n_af_inlet+1+n_bf_inlet-1 - (inlet_dv_coil_start-1) ){
			//if counter is at coils, that are downscaled for DV
			i_max[i] = RxBCurrDens*inlet_dv_scale;
		}
		else{
			i_max[i] = RxBCurrDens*bf_inlet_scale;
		}
		counter ++;
	}	

	int n_inlet = n_af_inlet + 1 + n_bf_inlet;
	
	// Helmholtz
	for (int i = N_coil + n_outlet + n_inlet; i < N_coil + n_outlet + n_inlet + 2; i++) {i_max[i] = RxBCurrDens*helm_scale;}	
	
	// Filter
	for (int i = N_coil + n_outlet + n_inlet + 2; i < N_coil + n_outlet + n_inlet + 2 + n_filter; i++) {i_max[i] = RxBCurrDens*filter_scale;}	
	
	// Connector after gate
	for (int i = N_coil + n_outlet + n_inlet + 2 + n_filter; i < N_coil + n_outlet + n_inlet + 2 + n_filter + n_af_connec; i++) {i_max[i] = RxBCurrDens*af_connec_scale;}	
	
	// Connector before gate
	for (int i = N_coil + n_outlet + n_inlet + 2 + n_filter + n_af_connec; i < N_coil + n_outlet + n_inlet + 2 + n_filter + n_af_connec + n_bf_connec; i++) {i_max[i] = RxBCurrDens*bf_connec_scale;}
	
	// Outlet Corr
	for (int i = N_coil + n_outlet + n_inlet + 2 + n_filter + n_af_connec + n_bf_connec; i < N_coil + n_outlet + n_inlet + 2 + n_filter + n_af_connec + n_bf_connec + n_outletCorr; i++) {
		i_max[i] = RxBCurrDens*outletCorr_scale;
	}	
	
	// INNER Inlet
	i_max[ N_coil + n_outlet + n_inlet + 2 + n_filter + n_af_connec + n_bf_connec + n_outletCorr ] = RxBCurrDens*inner_inlet_scale;
	
	

//////////////////////////// calc coil number

	double n_coil_nomos = N_coil + n_inlet + n_outlet + 2 + n_filter+ n_af_connec + n_bf_connec + n_outletCorr; // RxB, inlet, outlet, helm, filter, af, bf, outletCorr, innerInlet
	if (inner_inlet_scale != 0.) n_coil_nomos ++;
	
	cout << "MAG: Number of all NoMoS coils is: " << n_coil_nomos << endl;



/////////////////////////////////////////////////////////////////////
	//elec power + cooling calc
/////////////////////////////////////////////////////////////////////
	double copper_edge_radius = 0.00045; //edge radius in meter
	double edge_area = copper_edge_radius*copper_edge_radius*(4. - M_PI);
	double copper_area = RxB_conduct*RxB_conduct - water_dia*water_dia/4*M_PI - edge_area;
	if( Conductor == "SL" ){
		copper_area = SL_length*SL_thick;
	}
	//cout << "EL-POWER: Copper area = " << copper_area << endl;
	//cooling(r_coil,t_coil_RxB,l_coil_RxB,RxB_l*RxB_b,N_coil,water_dia,copper_area);
	double elecpowersum[11];
	
	// RxB
	elecpowersum[0] =(double)N_coil * elecPower(r_coil,r_coil+t_coil_RxB,RxBCurrent,RxB_l*RxB_b,copper_area);
	cout << "EL-POWER: 0. RxB-Coil-power = " << elecpowersum[0]/(double)N_coil << " Watt, total = " << elecpowersum[0] << " Watt" << endl;
	
	// Outlet
	elecpowersum[1]=(double)n_outlet * elecPower(outlet_inner,outlet_inner+t_coil_outlet,RxBCurrent*outlet_scale,outlet_l*outlet_b,copper_area);
	cout << "EL-POWER: 1. Outlet-Coil-power = " << elecpowersum[1]/(double)n_outlet << " Watt, total = " << elecpowersum[1] << " Watt" << endl;
	
	// AF Inlet
	elecpowersum[2]=((double)n_af_inlet) * elecPower(af_inlet_inner,af_inlet_inner+t_coil_af_inlet,RxBCurrent*af_inlet_scale,af_inlet_l*af_inlet_b,copper_area); 
	cout << "EL-POWER: 2. AF-Inlet-Coil-power = " << elecpowersum[2]/(double)n_af_inlet << " Watt, total = " << elecpowersum[2] << " Watt" << endl;
	// filter Inlet
	elecpowersum[3]= elecPower(filter_inlet_inner,filter_inlet_inner+t_coil_filter_inlet,RxBCurrent*filter_inlet_scale,filter_inlet_l*filter_inlet_b,copper_area); 
	cout << "EL-POWER: 3. Filter-Inlet-Coil-power = " << elecpowersum[3] << " Watt, total = " << elecpowersum[3] << " Watt" << endl;
	// BF Inlet
	elecpowersum[4]= (double)n_bf_inlet * elecPower(bf_inlet_inner,bf_inlet_inner+t_coil_bf_inlet,RxBCurrent*bf_inlet_scale,bf_inlet_l*bf_inlet_b,copper_area); 
	cout << "EL-POWER: 4. BF-Inlet-Coil-power = " << elecpowersum[4]/(double)n_bf_inlet << " Watt, total = " << elecpowersum[4] << " Watt" << endl;
	
	// Helm
	elecpowersum[5]=2. * elecPower(helm_inner,helm_inner+t_coil_helm,RxBCurrent*helm_scale,helm_l*helm_b,copper_area);
	cout << "EL-POWER: 5. Helmholtz-Coil-power = " << elecpowersum[5]/2. << " Watt, total = " << elecpowersum[5] << " Watt" << endl;
	
	// Filter
	if(filter_scale == 0.){elecpowersum[6] = 0.;}
	else{
	elecpowersum[6]=(double)n_filter * elecPower(filter_inner,filter_inner+t_coil_filter,RxBCurrent*filter_scale,filter_l*filter_b,copper_area);
	cout << "EL-POWER: 6. Filter-Coil-power = " << elecpowersum[6]/(double)n_filter << " Watt, total = " << elecpowersum[6] << " Watt" << endl;
	}
	
	// AF Connec
	elecpowersum[7]=(double)n_af_connec * elecPower(af_connec_inner,af_connec_inner+t_coil_af_connec,RxBCurrent*af_connec_scale,af_connec_l*af_connec_b,copper_area);
	cout << "EL-POWER: 7. AF-Connec-Coil-power = " << elecpowersum[7]/(double)n_af_connec << " Watt, total = " << elecpowersum[7] << " Watt" << endl;
	
	// BF Connec
	elecpowersum[8]=(double)n_bf_connec * elecPower(bf_connec_inner,bf_connec_inner+t_coil_bf_connec,RxBCurrent*bf_connec_scale,bf_connec_l*bf_connec_b,copper_area);
	cout << "EL-POWER: 8. BF-Connec-Coil-power = " << elecpowersum[8]/(double)n_bf_connec << " Watt, total = " << elecpowersum[8] << " Watt" << endl;
	
	// Outlet Corr
	elecpowersum[9]=(double)n_outletCorr * elecPower(outletCorr_inner,outletCorr_inner+t_coil_outletCorr,RxBCurrent*outletCorr_scale,outletCorr_l*outletCorr_b,copper_area);
	cout << "EL-POWER: 9. OutletCorr-Coil-power = " << elecpowersum[9]/(double)n_outletCorr << " Watt, total = " << elecpowersum[9] << " Watt" << endl;
	
	// INNER Inlet
	elecpowersum[10]= elecPower(inner_inlet_inner,inner_inlet_inner+t_coil_inner_inlet,RxBCurrent*inner_inlet_scale,inner_inlet_l*inner_inlet_b,copper_area);
	cout << "EL-POWER: 10. Inner-Inlet_Coil-power = " << elecpowersum[10] << " Watt, total = " << elecpowersum[10] << " Watt" << endl;
	

	double total=0.;
	for(int i = 0; i<11;i++){total+=elecpowersum[i];}
	cout << "EL-POWER: Total el. power consumption of NoMoS: " << total << " Watt" << endl;

	double currents[11] = {1.,outlet_scale,af_inlet_scale,filter_inlet_scale,bf_inlet_scale,helm_scale,filter_scale,af_connec_scale,bf_connec_scale,outletCorr_scale,inner_inlet_scale};
	for(int i=0; i<11; i++){ currents[i] = currents[i]* RxBCurrent;}
	double highestcurr=0.;
	int curr_position;
	for(int i=0;i<11;i++){
		if(highestcurr < currents[i]){highestcurr = currents[i]; curr_position = i;}
	}
	cout << "CURRDENS: Highest Curr. at pos: " << curr_position << ", CurrDens in COPPER: " << highestcurr/copper_area << "A/m*m"<< endl;
	cout << "CURRDENS: ... and CurrDens over coil (COOLING): " << highestcurr/(RxB_conduct*RxB_conduct) <<"A/m*m" <<  endl;



//////////////////// WANG Setting ////////////////////////////////
/*
	if (Wang_on == 1) 
		{
	    // Wangs parameters for his proposed RxB design (some)
		RxB = 0.06296;
	 	R_1 = 0.4;                              // RxB grand radius [m]
	    	d_j_max = 140./(0.0072*0.06296)*18.*3.;   // constant coil current in [A/m^2]
	    	i_N_coil = 21;                          // number of coils in the RxB
	    	r_coil = 0.15;                          // inner coil radius of all coils [m]
	   	r2_wire = 0.0072;                       // thickness of the RxB coils [m]
	   	d_space = 0;                            // minimal spacing between coils inside the RxB [m]
	   	n_filter = 0;                           // number of filter coils
	   	j_lin = 127./l_coil/t_coil*30.*3.;
	   	for (i_coil = 0; i_coil < i_N_coil; ++i_coil) {
	   	    if (i_coil == 0) {i_max[i_coil] = j_max/140.*132.;       // current density [A/m²]
	   	      } else if (i_coil == 1) {i_max[i_coil] = 133./140.*j_max;
	   	      } else if (i_coil == 2) {i_max[i_coil] =  135./140.*j_max;
	   	      } else if (i_coil >= 3 && i_coil <= 4) {i_max[i_coil] = 137./140.*j_max;
	   	      } else if (i_coil >= 5 && i_coil <= 7) {i_max[i_coil] = 138./140.*j_max;
	   	      }
	   	    }
	   	for(i_coil = i_N_coil; i_coil < i_N_coil + n_detector; i_coil++){
	   	    i_max[i_coil] = j_lin;} // current density
	   	for(i_coil = i_N_coil + n_detector; i_coil < i_N_coil + n_detector+n_apert; i_coil++) {
	   	    if(i_coil == i_N_coil + 3){i_max[i_coil] = 115./127.*j_lin;}
	   	    else if(i_coil == i_N_coil + 2){i_max[i_coil] = 119./127.*j_lin;}
	   	    else if(i_coil >= i_N_coil && i_coil < i_N_coil+2){i_max[i_coil] = 95./127.*j_lin;}
	   	    else{i_max[i_coil] = d_j_max;}         // current density [A/m²]
	        	}
		}
*/
	


// quickly write out file with R1 in it, so trajectory can read in, also aperture grid points for bline VS traj comparison

	ofstream configinfo;
	string infoname = filedir + "info.txt";
	configinfo.open(infoname.c_str(),ios::out);
	configinfo << "R1" << "\t" << R_1 << endl;
	double blineStartZ;
       	if( Perc_on == 0 ){
		blineStartZ = -l_coil_RxB/2. - RxBtoIn_d - (d_af_inlet+l_coil_af_inlet)*((double)n_af_inlet-1.) - l_coil_af_inlet;
	       blineStartZ = blineStartZ - FIntoAFIn_d - l_coil_filter_inlet - FIntoBFIn_d - (d_bf_inlet+l_coil_bf_inlet)*((double)n_bf_inlet-1.) - l_coil_bf_inlet - pumpport_distance/2.;
		configinfo << "blineStartZ" << "\t" << blineStartZ << endl;
	}// if on, written later, after PERC is built


//////////////////////////////////////////////////////////////////////
/////////////// inputcoil.dat creation
////////////////////////////////////////////////////////////////////////////

	string inputcoil_RxB = filedir + "inputcoil_RxB.dat";
	ofstream inputcoil_file;
	inputcoil_file.open(inputcoil_RxB.c_str(), ios::out);
	inputcoil_file << n_coil_nomos + Perc_N*Perc_on << "\n";
	inputcoil_file.close();
	

	
	    double starting_point[4];
	    double direction[4];
	double gate_begin;
	double connec_end;
	double filterinletstartZ = -l_coil_RxB/2. - RxBtoIn_d - (d_af_inlet+l_coil_af_inlet)*((double)n_af_inlet-1.) - l_coil_af_inlet - FIntoAFIn_d - l_coil_filter_inlet;
	double pump_end;

	if(nomos_on == 1)
	{

	    // RxB 
	    //double offset_angle= alpha/((double)N_coil-1.)*(double)n_cloth;       // angle of the first coil
	    starting_point[1] = 0.;
	    starting_point[2] = R_1;
	    starting_point[3] = 0.;
	    direction[1] = 0.;
	    direction[2] = 0.;
	    direction[3] = 1;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
		// input: current (vector), start coil number, starting point, direction, curv. Radius (vector), inner R, thickness,coil_length, curve angle, N_coil; filedir,shiftangle
	    circular_coils_to_inputcoil_dat(i_max, 0, starting_point,direction,starting_point,r_coil,t_coil_RxB,l_coil_RxB,alpha,N_coil, inputcoil_RxB,0.);

	    
		    // OUTLET coils
	    starting_point[1] = outlet_x_shift;
	    starting_point[2] = (R_1+outlet_r_shift)*cos(alpha)-(l_coil_RxB/2.+RxBtoOut_d)*sin(alpha);
	    starting_point[3] = (R_1+outlet_r_shift)*sin(alpha)+(l_coil_RxB/2.+RxBtoOut_d)*cos(alpha);
	    direction[1] = 0.;
	    direction[2] = -sin(alpha);
	    direction[3] = cos(alpha);
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
		// input: current, starting coil number, start P, direction, length, innerR, thick, gap, n_outlet, filedir, reverse direction yes =1
	    linear_coils_to_inputcoil_dat(i_max, N_coil, starting_point, direction, l_coil_outlet, outlet_inner, t_coil_outlet, d_outlet, n_outlet, inputcoil_RxB,0);


		    // OUTLET CORR coils
	    starting_point[1] = outletCorr_x_shift;
	    starting_point[2] = starting_point[2] - (n_outlet*l_coil_outlet+(n_outlet-1)*d_outlet)*sin(alpha) +(n_outletCorr*l_coil_outletCorr+(n_outletCorr-1)*d_outletCorr+outletCorr_fromoutEnd)*sin(alpha);
	    starting_point[3] = starting_point[3] + (n_outlet*l_coil_outlet+(n_outlet-1)*d_outlet)*cos(alpha) -(n_outletCorr*l_coil_outletCorr+(n_outletCorr-1)*d_outletCorr+outletCorr_fromoutEnd)*cos(alpha);
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
		// input: current, starting coil number, start P, direction, length, innerR, thick, gap, n_outlet, filedir, reverse direction yes =1
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+n_af_connec+n_bf_connec, starting_point, direction, l_coil_outletCorr, outletCorr_inner, t_coil_outletCorr, d_outletCorr, n_outletCorr, inputcoil_RxB,0);
	
	
		    // INLET
	    
		    // AF inlet
	    starting_point[1] = inlet_x_shift;
	    starting_point[2] = R_1+inlet_r_shift;
	    starting_point[3] = -l_coil_RxB/2.-RxBtoIn_d-l_coil_af_inlet;
	    direction[1] = 0.;
	    direction[2] = 0.;
	    direction[3] = 1.;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet, starting_point, direction, l_coil_af_inlet, bf_inlet_inner, t_coil_af_inlet, d_af_inlet, n_af_inlet, inputcoil_RxB,1);
		//reversed direction..	
	    

		// filter inlet
	    starting_point[3] = filterinletstartZ;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_af_inlet, starting_point, direction, l_coil_filter_inlet, filter_inlet_inner, t_coil_filter_inlet, d_af_inlet, 1, inputcoil_RxB,1);
	    

		    // BF inlet
	    starting_point[3] = filterinletstartZ - (FIntoBFIn_d + l_coil_bf_inlet);
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_af_inlet+1, starting_point, direction, l_coil_bf_inlet, bf_inlet_inner, t_coil_bf_inlet, d_bf_inlet, n_bf_inlet, inputcoil_RxB,1);
		//reversed direction..	

	
	    	// HELMHOLTZ coils
	    starting_point[1] = helm_x_shift;
	    starting_point[2] = R_1+helm_r_shift;
	    pump_end = filterinletstartZ - FIntoBFIn_d - ((double)n_bf_inlet-1.)*(d_bf_inlet+l_coil_bf_inlet) - l_coil_bf_inlet - HelmToBF_off; 
	    starting_point[3] = pump_end;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_inlet, starting_point, direction, l_coil_helm, helm_inner, t_coil_helm, pumpport_distance, 2, inputcoil_RxB,1);
	
	    
		// FILTER coils
	    starting_point[1] = filter_x_shift;
	    starting_point[2] = R_1 + filter_r_shift;
	    starting_point[3] = filterinletstartZ + l_coil_filter_inlet/2. - l_coil_filter/2.;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_inlet +2 , starting_point,direction, l_coil_filter, filter_inner, t_coil_filter, d_filter, n_filter, inputcoil_RxB,0);
	
	
	    // CONNECTOR coils AFTER gate
	    // before save gap for screwing
	    connec_end = pump_end - pumpport_distance - HelmToConnec;        
	    starting_point[1] = af_connec_x_shift;
	    starting_point[2] = R_1 + af_connec_r_shift;
	    starting_point[3] = connec_end - l_coil_af_connec;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter,starting_point,direction,l_coil_af_connec,af_connec_inner,t_coil_af_connec,d_af_connec,n_af_connec-1,inputcoil_RxB,1);

	    // first AF CONNEC coil with save screw gap
	    starting_point[3] = connec_end - l_coil_af_connec*(double)(n_af_connec -1) - (double)(n_af_connec -2)*d_af_connec - af_connec_savedist - l_coil_af_connec;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+n_af_connec-1,starting_point,direction,l_coil_af_connec,af_connec_inner,t_coil_af_connec,d_af_connec,1,inputcoil_RxB,1);



	    // CONNECTOR coils BEFORE gate
	    gate_begin = connec_end - l_coil_af_connec*(double)(n_af_connec -1) - (double)(n_af_connec -2)*d_af_connec - af_connec_savedist - l_coil_af_connec - gate_dist; 
	    starting_point[1] = bf_connec_x_shift;
	    starting_point[2] = R_1 + bf_connec_r_shift;
	    starting_point[3] = gate_begin - l_coil_bf_connec;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+n_af_connec,starting_point,direction,l_coil_bf_connec,bf_connec_inner,t_coil_bf_connec,bf_connec_savedist, 2, inputcoil_RxB,1);
	
	    
	    cout << "MAG: last startPz (BF_connec) = " << starting_point[3] << endl;	
	    	
	    // INNER INLET
	    if(inner_inlet_scale != 0.){
	    	starting_point[1] = inner_inlet_x_shift;
	    	starting_point[2] = R_1 + inner_inlet_r_shift;
	    	starting_point[3] = filterinletstartZ +l_coil_filter_inlet/2. - l_coil_inner_inlet/2.;
	    	// here, we use the starting point of filter coil, but calc the center of the filter coils, and go back half the inner inlet coil length so that the center is aligned with filter coils
	    	if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    	linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+n_af_connec+n_bf_connec+n_outletCorr, starting_point, direction, l_coil_inner_inlet, inner_inlet_inner, t_coil_inner_inlet, 0., 1, inputcoil_RxB,0);
	    }


	}


	if (Perc_on == 1) {
	        // PERC coils 
		direction[1] = 0.;
		direction[2] = 0.;
		direction[3] = 1.;
	        pump_end = filterinletstartZ - FIntoBFIn_d - ((double)n_bf_inlet-1.)*(d_bf_inlet+l_coil_bf_inlet) - l_coil_bf_inlet - HelmToBF_off; 
	        connec_end = pump_end - pumpport_distance - HelmToConnec;        // position of the upstream Pump port z-position
	        gate_begin = connec_end - l_coil_af_connec*(double)(n_af_connec -1) - (double)(n_af_connec -2)*d_af_connec - af_connec_savedist - l_coil_af_connec - gate_dist; 
	        starting_point[1] = 0.;
	        starting_point[2] = R_1+ bf_connec_r_shift;
	        starting_point[3] = gate_begin - l_coil_bf_connec*n_bf_connec - bf_connec_savedist - flansch_perc_thick - PERClastcoil_to_flansch;
	        if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	        PERC_coils_to_inputcoil_dat(starting_point, direction,perc_global_scale, perc_sol_scale, perc_filter_scale, filedir, inputcoil_RxB, Perc_coordinatsystem);
                
		configinfo << "PercEnd" << "\t" << starting_point[3] << endl;
		blineStartZ = starting_point[3] - (2.6825+0.9965/2.) -4.;
	}
	
	configinfo.close();
	

		///////////////////////////////////////////////////
	    // Run the simulation and calculate values of the B-field
		////////////////////////////////////////////////////
	input_coils(inputcoil_RxB);
	test_coils(filedir);
	magsource(filedir);



	
	if(coilcenterbool){	
		coil_center(R_1, alpha, -2., 0., 1.,horilines,vertilines,ApertX,ApertY, apertYshift,apertXshift, inputcoil_RxB, filedir, 0);     // R1, alpha,startPz,RxBStart,outlet-length,horil,vertil,apX,apY,inputc,filed,prep
	}

	if(geodriftbool){
		//fieldmap(2.*R_1, fieldmap_step, inputcoil_RxB, filedir, 0);       // map size [m], step_size [m], inputdatafile , initialize magfield (0=no)
		//radial_cuts(5, N_coil, R_1*2., 5e-3, inputcoil_RxB, filedir, 0); // lines, shift[rad], line_end[m], step_size [m], inputdatafile , initialize magfield (0=no)
		//vertical_cuts(5, N_coil, 1., 5e-3, R_1, inputcoil_RxB, filedir, 0); // lines, shift[rad], line_end[m], step_size [m], inputdatafile , initialize magfield (0=no)
		geoline_anadrift(R_1,alpha,horilines,vertilines,ApertX,ApertY,apertYshift,filedir,inputcoil_RxB,0,1187.29,true); // R1, alpha, horilines,vertilines,apertsizeX,apertsizeY,..,..,prep,momentum kev,drifton/off
		geoline_anadrift(R_1,alpha,horilines,vertilines,ApertX,ApertY,apertYshift,filedir,inputcoil_RxB,0,1187.29,false); // R1, alpha, horilines,vertilines,apertsizeX,apertsizeY,..,..,prep,momentum keV,drifton/off
	}

	if(blinebool){
		b_line_real(R_1,alpha,blineStartZ,horilines,vertilines,ApertX,ApertY,apertYshift,apertXshift,inputcoil_RxB,filedir,0,-1,1187.29,1,-1.); //..,..,Z,horil,vertil,apertx,aperty,..,..,prep,driftdir,momentum keV,driftOn, detposZ
		b_line_real(R_1,alpha,blineStartZ,horilines,vertilines,ApertX,ApertY,apertYshift,apertXshift,inputcoil_RxB,filedir,0,-1,1187.29,0,-1.); //..,..,Z,horil,vertil,apertx,aperty,..,..,prep,driftdir,mom,driftOn, detposZ
	}
	

	    // Run the simulation and calculate map the B-field, precise
	if(fieldmapbool){
		double fieldmapStartP[3];
		double fieldmapLength[3];
		if(Perc_on == 1) {
			fieldmapStartP[0] = -0.25;
			fieldmapStartP[1] = -R_1 - 0.4;
			fieldmapStartP[2] = -13.;
			fieldmapLength[0] = 0.5;
			fieldmapLength[1] = 2*R_1+2*0.4;
			fieldmapLength[2] = 13 + R_1 + 0.4;
		}
		else{
			fieldmapStartP[0] = -0.25;
			fieldmapStartP[1] = -R_1 - 0.4;
			fieldmapStartP[2] = -2.;
			fieldmapLength[0] = 0.5;
			fieldmapLength[1] = 2*R_1+2*0.4;
			fieldmapLength[2] = 2. + R_1 + 0.4;
		}
		fieldmap(fieldmapStartP,fieldmapLength, fieldmap_step, inputcoil_RxB, filedir, 0);       // map size [m], step_size [m], inputdatafile , initialize magfield (0=no)
	}



	
/*
	////// calculate B-field for each coil: In Root the coil currents will be adjusted and all solution merged
	    if(Fitting_on == 1){
	        ifstream inputcoil_file;
	        inputcoil_file.open(inputcoil_RxB.c_str(), ios::in);
	        double i_coil_max;
	        inputcoil_file >> i_coil_max;
	        for(i_coil = 1; i_coil <= i_coil_max-(Perc_N-2)*Perc_on; i_coil++){    // extract the coil datas in single input coils
	            // find the output path
	            stringstream outputcoil_stream;
	            outputcoil_stream << i_coil;
	            string outputcoil_name;
	            outputcoil_stream >> outputcoil_name;
	            if(i_coil == i_coil_max-(Perc_N-1)*Perc_on && Perc_on == 1){outputcoil_name="PERC";}
	            if(i_coil == i_coil_max-(Perc_N-2)*Perc_on && Perc_on == 1){outputcoil_name="PERConnect";}
	            string outputcoil_dir;
	            outputcoil_dir = filedir + "coils/" + outputcoil_name + "/";
	            outputcoil_name = outputcoil_dir + "inputcoil_" + outputcoil_name + ".dat";
	            // creat output file
	            ofstream outputcoil_file;
	            outputcoil_file.open(outputcoil_name.c_str(), ios::out);
	            double buffer;
	            // write first line
	            if(i_coil <= i_coil_max-Perc_N*Perc_on){
	                outputcoil_file << "1" << "\n";
	                if(i_coil <= N_coil+N_dummy || i_coil > N_coil + 2*N_dummy){outputcoil_file << "1000000" << "\t";}
	                else{outputcoil_file << "1000000" << "\t";}
	                inputcoil_file >> buffer;
	                }
	            else if (i_coil == i_coil_max-Perc_N*Perc_on + 1 && Perc_on == 1){outputcoil_file << Perc_N - 1  << "\n";}
	            else{outputcoil_file << "1" << "\n";}
	            // copy line of one coil
	            if(i_coil <= i_coil_max-Perc_N*Perc_on){
	                for(int j = 1; j <= 9; j++){
	                    inputcoil_file >> buffer;
	                    outputcoil_file << buffer << "\t";
	                    }
	                }
	            else if (i_coil == i_coil_max-Perc_N*Perc_on + 1 && Perc_on == 1){
	                for(int Perc_i = 0; Perc_i < Perc_N-1; Perc_i++){
	                    for(int j = 0; j <= 9; j++){
	                        inputcoil_file >> buffer;
	                        outputcoil_file << buffer << "\t";
	                        }
	                    outputcoil_file << "\n";
	                    }
	                }
	            else{                           // PERC connector
	                for(int Perc_i = 0; Perc_i < 1; Perc_i++){
	                    for(int j = 0; j <= 9; j++){
	                        inputcoil_file >> buffer;
	                        outputcoil_file << buffer << "\t";
	                        }
	                    outputcoil_file << "\n";
	                    }
	                }
	            outputcoil_file.close();
	            // calculate the B-field on the tangential lines to later on fit the currents
	            if(Fitting_on == 1 && i_coil >= 1){//N_coil+n_detector+n_apert){//
	                coil_center(d_R_1, alpha, -3, 0.0, 1., r_tube, b_beam, h_beam, outputcoil_name, outputcoil_dir, 1);     // Radius[m], Angle[rad], Start_pos, Bend_Start, r1_wire_detect
	                }
	            }
	        inputcoil_file.close();
	        }
	
	}}} }}} }}}       // close all possible sweep loops
*/

	return 0;

}
