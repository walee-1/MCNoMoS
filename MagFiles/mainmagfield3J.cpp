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
#define Ncoilmax 130
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
{

	//directory from program call
	if (argc <2) {cout << "MAG: no directory given!" << endl; return 1;}
	string filedir =argv[1];
	cout <<"MAG: dir given: " <<  filedir << endl;


/////////////Config parser////////////// 
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
	bool coilcenterbool= myconfig.pBool("coilcenter");
	bool geodriftbool = myconfig.pBool("geo+drift");
	bool blinebool = myconfig.pBool("bline");
	bool startinPERC = myconfig.pBool("startinPERC");
	bool fieldmapbool = myconfig.pBool("fieldmap");
	double apertYshift = myconfig.pDouble("apertYshift");
	double apertXshift = myconfig.pDouble("apertXshift");
	bool onlyCornerCenter = myconfig.pBool("onlyCornerCenter");
	int horilines = myconfig.pInt("horilines");
	int vertilines = myconfig.pInt("vertilines");
	double ApertX = myconfig.pDouble("ApertX");
	double ApertY = myconfig.pDouble("ApertY");
	double detpos = myconfig.pDouble("detpos");
	
	bool Perc_on = myconfig.pBool("PercOn");
	double perc_global_scale = myconfig.pDouble("perc_global_scale");
	double perc_sol_scale = myconfig.pDouble("perc_sol_scale");
	double perc_filter_scale = myconfig.pDouble("perc_filter_scale");
	double perc_connec_scale = myconfig.pDouble("perc_connec_scale");
	int NoMoSOn = myconfig.pBool("NoMoSOn");
	string Conductor = myconfig.pString("Conductor");
	double conduct_s = myconfig.pDouble("conduct_s");
	double Raisin = myconfig.pDouble("Raisin");
	double water_dia = myconfig.pDouble("water_dia");
	double i_max[130];                                          // real current densities for all coils.
	double PERClastcoil_to_flansch = myconfig.pDouble("PERClastcoil_to_flansch");
	double flansch_perc_thick = myconfig.pDouble("flansch_perc_thick");
	double fieldmap_step = myconfig.pDouble("fieldmap_step");
	double manualStartZ = myconfig.pDouble("manualStartZ");
	bool boolManualStartZ = myconfig.pBool("boolManualStartZ");

	// RxB Parameter
	int RxB_l = myconfig.pInt("RxB_l");
	int RxB_b = myconfig.pInt("RxB_b");
	double d_space = myconfig.pDouble("d_RxB");//simulation_parameters[5];       // minimal spacing between coils inside the RxB [m]
	int N_coil = myconfig.pInt("n_RxB");//simulation_parameters[2];                 // number of coils in the RxB
	double alpha = myconfig.pDouble("alpha")/180.*M_PI;//simulation_parameters[6];   // opening angle of the RxB [rad]
	double RxBCurrent = myconfig.pDouble("RxBCurrent");
	double RxB_second_scale = myconfig.pDouble("RxB_second_scale");
	double r_coil = myconfig.pDouble("RxBInner");//simulation_parameters[3];      // inner coil radius of all coils [m]
	double FirstRxB_scale = myconfig.pDouble("FirstRxB_scale");
	double LastRxB_scale = myconfig.pDouble("LastRxB_scale");
	// outer RxB
	int oRxB_l = myconfig.pInt("oRxB_l");
	int oRxB_b = myconfig.pInt("oRxB_b");
	double R_shift = myconfig.pDouble("R_shift");
	int No_coil = N_coil -1;
	double alphastart = alpha/(N_coil-1)/2.;
	double alpha2 = alpha-alphastart*2;
	double oRxB_scale = myconfig.pDouble("oRxB_curr_scale");
	double ro_coil = myconfig.pDouble("oRxBInner");//simulation_parameters[3];      // inner coil radius of all coils [m]
	bool oRxB_on = myconfig.pBool("oRxB_on");

	// trans Corr
	int trans_l = myconfig.pInt("trans_l");//simulation_parameters[18];            // linear coil length [m]
	int trans_b = myconfig.pInt("trans_b");//r2_wire;//simulation_parameters[19]; // linear coil thickness [m]
	double trans_r_shift = myconfig.pDouble("trans_r_shift");             // s_coil_dec radial shift of all detector coils
	double trans_scale = myconfig.pDouble("trans_curr_scale");
	double trans_inner = myconfig.pDouble("trans_inner");
	double trans_angle = myconfig.pDouble("trans_angle");
	trans_angle = trans_angle/180.*M_PI;

	// Outlet Coils (Detector)
	int outlet_l = myconfig.pInt("outlet_l");//simulation_parameters[18];            // linear coil length [m]
	int outlet_b = myconfig.pInt("outlet_b");//r2_wire;//simulation_parameters[19]; // linear coil thickness [m]
	double d_outlet = myconfig.pDouble("d_outlet");//0.005;//simulation_parameters[20];   // distance between linear coils [m]
	int n_outlet = myconfig.pInt("n_outlet");                // number of coils in the detector region (without correction coils)
	double outlet_r_shift = myconfig.pDouble("outlet_r_shift");             // s_coil_dec radial shift of all detector coils
	double outlet_x_shift = myconfig.pDouble("outlet_x_shift");
	double outlet_scale = myconfig.pDouble("outlet_curr_scale");
	double outlet_inner = myconfig.pDouble("outlet_inner");
	double outlet_first_scale = myconfig.pDouble("outlet_first_scale");
	double outlet_screw_gap = myconfig.pDouble("outlet_screw_gap");

	// Outlet Correction
	int n_outletCorr = myconfig.pInt("n_outletCorr");
	double outletCorr_l = myconfig.pInt("outletCorr_l");
	double outletCorr_b = myconfig.pInt("outletCorr_b");
	double d_outletCorr = myconfig.pDouble("d_outletCorr");
	double outletCorr_r_shift = myconfig.pDouble("outletCorr_r_shift");             // s_coil_dec radial shift of all detector coils
	double outletCorr_x_shift = myconfig.pDouble("outletCorr_x_shift");
	double outletCorr_scale = myconfig.pDouble("outletCorr_curr_scale");
	double outletCorr_inner_gap = myconfig.pDouble("outletCorr_inner_gap"); // gap between outer radius of underlying coil and helmholtz coil
	double outletCorr_fromoutEnd = myconfig.pDouble("outletCorr_fromoutEnd");
	
	// Exit Corr
	int exit_l = myconfig.pInt("exit_l");
	int exit_b = myconfig.pInt("exit_b");
	double exit_r_shift = myconfig.pDouble("exit_r_shift");
	double exit_x_shift = myconfig.pDouble("exit_x_shift");
	double exit_scale = myconfig.pDouble("exit_curr_scale");
	double exit_inner = myconfig.pDouble("exit_inner");
	double exit_angle = myconfig.pDouble("exit_angle")/180.*M_PI;
	double exit_z_shift = myconfig.pDouble("exit_z_shift");


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
	double af_inlet_scale = myconfig.pDouble("af_inlet_curr_scale");
	double inlet_last_scale = myconfig.pDouble("inlet_last_scale");
	double af_screw_gap = myconfig.pDouble("af_screw_gap");
	
	// BF Inlet Coils
	int bf_inlet_l = myconfig.pInt("bf_inlet_l");
	int bf_inlet_b = myconfig.pInt("bf_inlet_b");
	double d_bf_inlet = myconfig.pDouble("d_bf_inlet");
	int n_bf_inlet = myconfig.pInt("n_bf_inlet");
	double bf_inlet_scale = myconfig.pDouble("bf_inlet_curr_scale");
	double bf_screw_gap = myconfig.pDouble("bf_screw_gap");

	double inlet_dv_scale = myconfig.pDouble("inlet_dv_coil_scale");
	int inlet_dv_coil_start = myconfig.pInt("inlet_dv_coil_start");
	int inlet_dv_offN = myconfig.pInt("inlet_dv_offN");
	
	// filter Inlet coil
	int filter_inlet_l = myconfig.pInt("filter_inlet_l");
	int filter_inlet_b = myconfig.pInt("filter_inlet_b");
	double filter_inlet_scale = myconfig.pDouble("filter_inlet_curr_scale");

	// enter Corr
	int enter_l = myconfig.pInt("enter_l");
	int enter_b = myconfig.pInt("enter_b");
	double enter_r_shift = myconfig.pDouble("enter_r_shift");
	double enter_x_shift = myconfig.pDouble("enter_x_shift");
	double enter_scale = myconfig.pDouble("enter_curr_scale");
	double enter_inner = myconfig.pDouble("enter_inner");
	double enter_angle = myconfig.pDouble("enter_angle")/180.*M_PI;
	double enter_z_shift = myconfig.pDouble("enter_z_shift");


	// filter coils parameters
	int n_filter = myconfig.pInt("n_filter");               // number of filter coils
	int filter_l = myconfig.pInt("filter_l");
	int filter_b = myconfig.pInt("filter_b");
	double d_filter = myconfig.pDouble("d_filter");
	double filter_inner_gap = myconfig.pDouble("filter_inner_gap"); // gap between outer radius of underlying coil and helmholtz coil
	double filter_r_shift = myconfig.pDouble("filter_r_shift");           // shift of all coils in the aperture region in the y direction
	double filter_x_shift = myconfig.pDouble("filter_x_shift");
	double filter_scale = myconfig.pDouble("filter_curr_scale");
	double filter_z_shift = myconfig.pDouble("filter_z_shift");

	// Helmholtz coils 1 and 2
	int helm1_l = myconfig.pInt("helm1_l");
	int helm1_b = myconfig.pInt("helm1_b");
	double helm1_inner = myconfig.pDouble("helm1_inner"); // gap between outer radius of underlying coil and helmholtz coil
	double helm_r_shift = myconfig.pDouble("helm_r_shift");           // shift of all coils in the aperture region in the y direction
	double helm_x_shift = myconfig.pDouble("helm_x_shift");
	double helm1_scale = myconfig.pDouble("helm1_curr_scale");
	double helm1_bf_gap = myconfig.pDouble("helm1_bf_gap");
	
	double NC_Curr = myconfig.pDouble("NC_Curr");
	int helm2_l = myconfig.pInt("helm2_l");
	int helm2_b = myconfig.pInt("helm2_b");
	double helm2_inner = myconfig.pDouble("helm2_inner"); // gap between outer radius of underlying coil and helmholtz coil
	double helm2_scale = myconfig.pDouble("helm2_curr_scale");



	// New Connector coils
	// GC connec
	double pumpport_distance = myconfig.pDouble("pumpport_dist");//simulation_parameters[29];    // Distance between adapter coils and aperture coils of NoMoS [m]
	int gc_connec_l = myconfig.pInt("gc_connec_l");
	int gc_connec_b = myconfig.pInt("gc_connec_b");
	double gc_connec_r_shift = myconfig.pDouble("gc_connec_r_shift");           // shift of all coils in the aperture region in the y direction
	double gc_connec_x_shift = myconfig.pDouble("gc_connec_x_shift");
	double gc_connec_scale = myconfig.pDouble("gc_connec_curr_scale");
	double gc_connec_inner = myconfig.pDouble("gc_connec_inner");
	double gate_dist = myconfig.pDouble("gate_dist");

	// PreGate connec
	int pg_connec_l = myconfig.pInt("pg_connec_l");
	int pg_connec_b = myconfig.pInt("pg_connec_b");
	double pg_connec_r_shift = myconfig.pDouble("pg_connec_r_shift");           // shift of all coils in the aperture region in the y direction
	double pg_connec_x_shift = myconfig.pDouble("pg_connec_x_shift");
	double pg_connec_scale = myconfig.pDouble("pg_connec_curr_scale");
	double pg_connec_inner = myconfig.pDouble("pg_connec_inner");
	double pipe_dist = myconfig.pDouble("pipe_dist");

	// PrePipe connec
	double screw_dist = myconfig.pDouble("screw_dist");
	int pp_connec_l = myconfig.pInt("pp_connec_l");
	int pp_connec_b = myconfig.pInt("pp_connec_b");
	double pp_connec_r_shift = myconfig.pDouble("pp_connec_r_shift");           // shift of all coils in the aperture region in the y direction
	double pp_connec_x_shift = myconfig.pDouble("pp_connec_x_shift");
	double pp_connec_scale = myconfig.pDouble("pp_connec_curr_scale");
	double pp_connec_inner = myconfig.pDouble("pp_connec_inner");

	// Perc Parameterrs:
	int Perc_N = 13;                // number of coils of the PERC system
	


	
	///////////////////////////////////////////////////////////////
////////////////////////////// Supra leiter coil parameters ///////////////////// or parameters depending on technology
	////////////////////////////////////////////////////////////////// 
	
	// some distance parameters that are normally at standard values
	double RxBtoOut_d = d_outlet; // gap between last RxB coil and outlet coil
	double RxBtoIn_d = d_af_inlet; // same for Inlet
	double FIntoAFIn_d = d_af_inlet;
	double FIntoBFIn_d = d_bf_inlet;
	double HelmToBF_off = 0.;
	double HelmToConnec = 0.;
	double raisin;	
	
	double SL_length, SL_thick;
	
	if( Conductor == "NL" ){
		SL_length = conduct_s;
		SL_thick = conduct_s;
		raisin = Raisin; // get value from config file
	}
	if ( Conductor == "HTSL" ){
		SL_length = 0.0045;
		SL_thick = 0.00065;
		raisin = 0.; // HTSL sizes are already with isolation, therefore raisin = 0
	}
	if( Conductor == "TTSL" ){
		SL_length = 0.0005;
		SL_thick = 0.0005;
		raisin = 0.0019 - SL_length;
	}

	// now we set the sizes of the coils according to the previously set parameters for the conductor type

	// RxB coils 
	double l_coil_RxB =  SL_length*RxB_l + (RxB_l - 1)*raisin;
	double t_coil_RxB = SL_thick*RxB_b + (RxB_b - 1)*raisin;
	double RxBCurrDens = RxBCurrent*RxB_l*RxB_b/(l_coil_RxB * t_coil_RxB); // same current * turns but on smaller area -> much higher density
	double R_1 = r_coil + t_coil_RxB + 2.*raisin + ( (l_coil_RxB + 2*raisin)/2. + d_space ) / tan ( alpha/(2*(N_coil -1))  );
	// outer RxB
	double l_coil_oRxB =  SL_length*oRxB_l + (oRxB_l - 1)*raisin;
	double t_coil_oRxB = SL_thick*oRxB_b + (oRxB_b - 1)*raisin;
	double oRxBCurrDens = RxBCurrent*oRxB_l*oRxB_b/(l_coil_oRxB * t_coil_oRxB); 
	double R_2 = R_1 + R_shift;

	// RxB trans Corr
	double l_coil_trans =  SL_length*trans_l + (trans_l - 1)*raisin;
	double t_coil_trans = SL_thick*trans_b + (trans_b - 1)*raisin;
	double transCurrDens = RxBCurrent*trans_l*trans_b/(l_coil_trans * t_coil_trans); 

	// Outlet
	double l_coil_outlet = SL_length*outlet_l + (outlet_l -1)*raisin;
	double t_coil_outlet = SL_thick*outlet_b + (outlet_b -1)*raisin;
	double outletCurrDens = RxBCurrent*outlet_l*outlet_b/(l_coil_outlet * t_coil_outlet); 

	// Outlet Corr
	double l_coil_outletCorr = SL_length*outletCorr_l + (outletCorr_l -1)*raisin;
	double t_coil_outletCorr = SL_thick*outletCorr_b + (outletCorr_b -1)*raisin;
	double outletCorr_inner = outlet_inner + t_coil_outlet + outletCorr_inner_gap; // with new t, get new inner R
	double outletCorrCurrDens = RxBCurrent*outletCorr_l*outletCorr_b/(l_coil_outletCorr * t_coil_outletCorr); 

	// Exit
	double l_coil_exit = SL_length*exit_l + (exit_l -1)*raisin;
	double t_coil_exit = SL_thick*exit_b + (exit_b -1)*raisin;
	double exitCurrDens = RxBCurrent*exit_l*exit_b/(l_coil_exit * t_coil_exit); 


	// AF Inlet
	double l_coil_af_inlet = SL_length*af_inlet_l + (af_inlet_l -1)*raisin;
	double t_coil_af_inlet = SL_thick*af_inlet_b + (af_inlet_b -1)*raisin;
	double af_inletCurrDens = RxBCurrent*af_inlet_l*af_inlet_b/(l_coil_af_inlet * t_coil_af_inlet); 

	// BF Inlet
	double l_coil_bf_inlet = SL_length*bf_inlet_l +(bf_inlet_l -1)*raisin;
	double t_coil_bf_inlet = SL_thick*bf_inlet_b +(bf_inlet_b -1)*raisin;
	double bf_inletCurrDens = RxBCurrent*bf_inlet_l*bf_inlet_b/(l_coil_bf_inlet * t_coil_bf_inlet); 

	// Enter
	double l_coil_enter = SL_length*enter_l +(enter_l -1)*raisin;
	double t_coil_enter = SL_thick*enter_b +(enter_b -1)*raisin;
	double enterCurrDens = RxBCurrent*enter_l*enter_b/(l_coil_enter * t_coil_enter); 

	// filter Inlet
	double l_coil_filter_inlet = SL_length*filter_inlet_l +(filter_inlet_l -1)*raisin;
	double t_coil_filter_inlet = SL_thick*filter_inlet_b +(filter_inlet_b -1)*raisin;
	double filter_inletCurrDens = RxBCurrent*filter_inlet_l*filter_inlet_b/(l_coil_filter_inlet * t_coil_filter_inlet); 

	// filter
	double l_coil_filter = SL_length*filter_l +(filter_l -1)*raisin;
	double t_coil_filter = SL_thick*filter_b +(filter_b -1)*raisin;
	double filter_inner = inlet_inner + t_coil_filter_inlet + filter_inner_gap;
	double filterCurrDens = RxBCurrent*filter_l*filter_b/(l_coil_filter * t_coil_filter); 

	// helmholtz
	double l_coil_helm1 = SL_length*helm1_l +(helm1_l -1)*raisin;
	double t_coil_helm1 = SL_thick*helm1_b +(helm1_b -1)*raisin;
	double helm1CurrDens = RxBCurrent*helm1_l*helm1_b/(l_coil_helm1 * t_coil_helm1); 

	
	// if we use Super conductor, the following coils use normal conductor sizes still
	if( Conductor == "TTSL" || Conductor == "HTSL" ){
		SL_length = conduct_s;
		SL_thick = conduct_s;
		raisin = Raisin;
	}

	double l_coil_helm2 = SL_length*helm2_l +(helm2_l -1)*raisin;
	double t_coil_helm2 = SL_thick*helm2_b +(helm2_b -1)*raisin;
	double helm2_dens = NC_Curr*helm2_b*helm2_l/(l_coil_helm2*t_coil_helm2);

	//  Connector 1
	double l_coil_gc_connec = SL_length*gc_connec_l +(gc_connec_l -1)*raisin;
	double t_coil_gc_connec = SL_thick*gc_connec_b +(gc_connec_b -1)*raisin;
	double gc_dens = NC_Curr*gc_connec_b*gc_connec_l/(l_coil_gc_connec*t_coil_gc_connec);
	//  Connector 2
	double l_coil_pg_connec = SL_length*pg_connec_l +(pg_connec_l -1)*raisin;
	double t_coil_pg_connec = SL_thick*pg_connec_b +(pg_connec_b -1)*raisin;
	double pg_dens = NC_Curr*pg_connec_b*pg_connec_l/(l_coil_pg_connec*t_coil_pg_connec);
	//  Connector 3
	double l_coil_pp_connec = SL_length*pp_connec_l +(pp_connec_l -1)*raisin;
	double t_coil_pp_connec = SL_thick*pp_connec_b +(pp_connec_b -1)*raisin;
	double pp_dens = NC_Curr*pp_connec_b*pp_connec_l/(l_coil_pp_connec*t_coil_pp_connec);


	cout << "MAG: R_1 = " << R_1 << " m" << endl; 
	cout << "MAG: RxBCurrDens = " << RxBCurrDens << " A/m*m" << endl;	

	
/////////////////////////////////////////////////////////////////////
//////////////////////////////	// current set /////////////////////////////
/////////////////////////////////////////////////////////////////////
	
	int counter = 0;

	// RxB
	for (int i = 0; i < N_coil; i++) {
		if(i == 0 ){ // First or Last RxB Coil
			i_max[i] = RxBCurrDens * FirstRxB_scale;
		}
		else if(i == N_coil -1) i_max[i] = RxBCurrDens * LastRxB_scale;
		else if( i >= (N_coil-1)/2 ){
			i_max[i] = RxBCurrDens * RxB_second_scale;
		}
		else{
			i_max[i] = RxBCurrDens;
		}
	}	
	
	// Outlet
	for (int i = N_coil; i < N_coil + n_outlet; i++) {
		
		if( i == N_coil) i_max[i] = outletCurrDens * outlet_first_scale;
		else i_max[i] = outletCurrDens*outlet_scale;
	
	}
	
	// Inlet
	// AF inlet
	for (int i = N_coil + n_outlet; i < N_coil + n_outlet + n_af_inlet; i++) {
				if( i == N_coil + n_outlet ) i_max[i] = af_inletCurrDens * inlet_last_scale;
		else i_max[i] = af_inletCurrDens*af_inlet_scale;
		counter ++;
	}	

	// Inlet Filter coil
	i_max[ N_coil + n_outlet + n_af_inlet ] = filter_inletCurrDens*filter_inlet_scale;
	counter ++;
	// BF inlet
	for (int i = N_coil + n_outlet + n_af_inlet + 1; i < N_coil + n_outlet + n_af_inlet + 1 + n_bf_inlet; i++) {
		if( counter >= n_af_inlet+1+n_bf_inlet-1 - (inlet_dv_coil_start-1) - (inlet_dv_offN-1) && counter <= n_af_inlet+1+n_bf_inlet-1 - (inlet_dv_coil_start-1) ){
			//if counter is at coils, that are downscaled for DV
			i_max[i] = bf_inletCurrDens*inlet_dv_scale;
		}
		else{
			i_max[i] = bf_inletCurrDens*bf_inlet_scale;
		}
		counter ++;
	}	

	int n_inlet = n_af_inlet + 1 + n_bf_inlet;
	
	// Helmholtz
	i_max[N_coil + n_outlet + n_inlet] = helm1CurrDens*helm1_scale;
	i_max[N_coil + n_outlet + n_inlet + 1] = helm2_dens*helm2_scale;
	
	// Filter
	if(n_filter > 0){
		for (int i = N_coil + n_outlet + n_inlet + 2; i < N_coil + n_outlet + n_inlet + 2 + n_filter; i++) {i_max[i] = filterCurrDens*filter_scale;}	
	}
	
// here we adapted quickly the numbers so that we can have 3 NC coils as connectors

	// Connector 1
	for (int i = N_coil + n_outlet + n_inlet + 2 + n_filter ; i < N_coil + n_outlet + n_inlet + 2 + n_filter + 1; i++) {i_max[i] = gc_dens*gc_connec_scale;}	
	
	// Connector 2
	for (int i = N_coil + n_outlet + n_inlet + 2 + n_filter + 1; i < N_coil + n_outlet + n_inlet + 2 + n_filter + 2; i++) {i_max[i] = pg_dens*pg_connec_scale;}
	
	// Connector 3
	for (int i = N_coil + n_outlet + n_inlet + 2 + n_filter + 2; i < N_coil + n_outlet + n_inlet + 2 + n_filter + 3; i++) {i_max[i] = pp_dens*pp_connec_scale;}

	// Outlet Corr
	for (int i = N_coil + n_outlet + n_inlet + 2 + n_filter + 3; i < N_coil + n_outlet + n_inlet + 2 + n_filter + 3 + n_outletCorr; i++) {
		i_max[i] = outletCorrCurrDens*outletCorr_scale;
	}	
	
	// enter Corr
	i_max[ N_coil + n_outlet + n_inlet + 2 + n_filter + 3 + n_outletCorr ] = enterCurrDens*enter_scale;
	
	// exit Corr
	i_max[ N_coil + n_outlet + n_inlet + 2 + n_filter + 3 + n_outletCorr + 1 ] = exitCurrDens*exit_scale;
	
	// trans Corr
	i_max[ N_coil + n_outlet + n_inlet + 2 + n_filter + 3 + n_outletCorr + 2 ] = transCurrDens*trans_scale;

	for(int i= N_coil + n_outlet + n_inlet + 2 + n_filter + 3 + n_outletCorr + 2 +1; i < N_coil + n_outlet + n_inlet + 2 + n_filter + 3 + n_outletCorr + 2 +1 +No_coil; i++){
		i_max[i]= oRxB_scale * oRxBCurrDens;
	}



//////////////////////////// calc coil number
// RxB, inlet, outlet, helm, filter, hhelm, outletCorr, first connector
	double n_coil_nomos = N_coil + n_inlet + n_outlet + 2 + n_filter + n_outletCorr +1;
	// first two connectors off or on
 	if( pg_connec_scale != 0. && pp_connec_scale != 0. ) n_coil_nomos += 2;
	// enterCorr, exit Corr, trans Corr
	if( trans_scale == 0. ) n_coil_nomos += 2;
	else n_coil_nomos += 3;	
	if(oRxB_on) n_coil_nomos += No_coil;
	cout << "MAG: Number of all NoMoS coils is: " << n_coil_nomos << endl;


///////////////////////////////////////////////////////////////////
      //elec power + cooling calc
///////////////////////////////////////////////////////////////////
	double total = 0.;
	if( Conductor == "NL" ){
		double copper_edge_radius = 0.001; //edge radius in meter
		double edge_area = copper_edge_radius*copper_edge_radius*(4. - M_PI);
		double copper_area = conduct_s*conduct_s - water_dia*water_dia/4*M_PI - edge_area;
		
		//cout << "PowerBug: copper area: " << copper_area << endl;
		//cout << "EL-POWER: Copper area = " << copper_area << endl;
		//cooling(r_coil,t_coil_RxB,l_coil_RxB,RxB_l*RxB_b,N_coil,water_dia,copper_area);
		double elecpowersum[13];
		
		// RxB
		elecpowersum[0] =((double)N_coil-2.) * elecPower(r_coil,r_coil+t_coil_RxB,RxBCurrent,RxB_l*RxB_b,copper_area);
	       	elecpowersum[0] += elecPower(r_coil,r_coil+t_coil_RxB,FirstRxB_scale*RxBCurrent,RxB_l*RxB_b,copper_area);
	       	elecpowersum[0] += elecPower(r_coil,r_coil+t_coil_RxB,LastRxB_scale*RxBCurrent,RxB_l*RxB_b,copper_area);
		cout << "EL-POWER: 0. RxB-Coil-power = " << elecpowersum[0]/(double)N_coil << " Watt, total = " << elecpowersum[0] << " Watt" << endl;
		
		// Outlet
		elecpowersum[1]=((double)n_outlet-1.) * elecPower(outlet_inner,outlet_inner+t_coil_outlet,RxBCurrent*outlet_scale,outlet_l*outlet_b,copper_area);
		elecpowersum[1] += elecPower(outlet_inner,outlet_inner+t_coil_outlet,RxBCurrent*outlet_first_scale,outlet_l*outlet_b,copper_area);
		cout << "EL-POWER: 1. Outlet-Coil-power = " << elecpowersum[1]/(double)n_outlet << " Watt, total = " << elecpowersum[1] << " Watt" << endl;
		
		// AF Inlet
		// same procedure as with RxB, last coil factors in
		elecpowersum[2]=((double)n_af_inlet -1.) * elecPower(inlet_inner,inlet_inner+t_coil_af_inlet,RxBCurrent*af_inlet_scale,af_inlet_l*af_inlet_b,copper_area); 
		elecpowersum[2] += elecPower(inlet_inner,inlet_inner+t_coil_af_inlet,RxBCurrent*inlet_last_scale, af_inlet_l*af_inlet_b,copper_area);
		cout << "EL-POWER: 2. AF-Inlet-Coil-power = " << elecpowersum[2]/(double)n_af_inlet << " Watt, total = " << elecpowersum[2] << " Watt" << endl;
		// filter Inlet
		elecpowersum[3]= elecPower(inlet_inner,inlet_inner+t_coil_filter_inlet,RxBCurrent*filter_inlet_scale,filter_inlet_l*filter_inlet_b,copper_area); 
		cout << "EL-POWER: 3. Filter-Inlet-Coil-power = " << elecpowersum[3] << " Watt, total = " << elecpowersum[3] << " Watt" << endl;
		// BF Inlet
		elecpowersum[4]= ((double)n_bf_inlet - (double)inlet_dv_offN) * elecPower(inlet_inner,inlet_inner+t_coil_bf_inlet,RxBCurrent*bf_inlet_scale,bf_inlet_l*bf_inlet_b,copper_area); 
		elecpowersum[4] += (double)inlet_dv_offN * elecPower(inlet_inner,inlet_inner+t_coil_bf_inlet,RxBCurrent*inlet_dv_scale,bf_inlet_l*bf_inlet_b,copper_area);
		cout << "EL-POWER: 4. BF-Inlet-Coil-power = " << elecpowersum[4]/(double)n_bf_inlet << " Watt, total = " << elecpowersum[4] << " Watt" << endl;
		
		// Helm
		elecpowersum[5]= elecPower(helm1_inner,helm1_inner+t_coil_helm1,RxBCurrent*helm1_scale,helm1_l*helm1_b,copper_area);
		elecpowersum[5] += elecPower(helm2_inner,helm2_inner+t_coil_helm2,NC_Curr*helm2_scale,helm2_l*helm2_b,copper_area);
		cout << "EL-POWER: 5. Helmholtz-Coil-power = " << elecpowersum[5]/2. << " Watt, total = " << elecpowersum[5] << " Watt" << endl;
		
		// Filter
		if(filter_scale == 0.){elecpowersum[6] = 0.;}
		else{
		elecpowersum[6]=(double)n_filter * elecPower(filter_inner,filter_inner+t_coil_filter,RxBCurrent*filter_scale,filter_l*filter_b,copper_area);
		cout << "EL-POWER: 6. Filter-Coil-power = " << elecpowersum[6]/(double)n_filter << " Watt, total = " << elecpowersum[6] << " Watt" << endl;
		}
		
		
		// All Connec MEAN
		elecpowersum[7]= elecPower(gc_connec_inner,gc_connec_inner+t_coil_gc_connec,NC_Curr*gc_connec_scale,gc_connec_l*gc_connec_b,copper_area);
		elecpowersum[7]+= elecPower(pg_connec_inner,pg_connec_inner+t_coil_pg_connec,NC_Curr*pg_connec_scale,pg_connec_l*pg_connec_b,copper_area);
		elecpowersum[7]+= elecPower(pp_connec_inner,pp_connec_inner+t_coil_pp_connec,NC_Curr*pp_connec_scale,pp_connec_l*pp_connec_b,copper_area);
		cout << "EL-POWER: 7. Mean Connec-Coil-power = " << elecpowersum[7]/3. << " Watt, total = " << elecpowersum[7] << " Watt" << endl;
		
		// Outlet Corr
		elecpowersum[8]=(double)n_outletCorr * elecPower(outletCorr_inner,outletCorr_inner+t_coil_outletCorr,RxBCurrent*outletCorr_scale,outletCorr_l*outletCorr_b,copper_area);
		cout << "EL-POWER: 8. OutletCorr-Coil-power = " << elecpowersum[8]/(double)n_outletCorr << " Watt, total = " << elecpowersum[8] << " Watt" << endl;
		
		// enter corr
		elecpowersum[9]= elecPower(enter_inner,enter_inner+t_coil_enter,RxBCurrent*enter_scale,enter_l*enter_b,copper_area);
		cout << "EL-POWER: 9. Enter_Coil-power = " << elecpowersum[9] << " Watt, total = " << elecpowersum[9] << " Watt" << endl;
		
		// exit corr
		elecpowersum[10]= elecPower(exit_inner,exit_inner+t_coil_exit,RxBCurrent*exit_scale,exit_l*exit_b,copper_area);
		cout << "EL-POWER: 10. Exit_Coil-power = " << elecpowersum[10] << " Watt, total = " << elecpowersum[10] << " Watt" << endl;
	
		// trans corr
		elecpowersum[11]= elecPower(trans_inner,trans_inner+t_coil_trans,RxBCurrent*trans_scale,trans_l*trans_b,copper_area);
		cout << "EL-POWER: 11. Trans_Coil-power = " << elecpowersum[11] << " Watt, total = " << elecpowersum[11] << " Watt" << endl;

		// outer RxB
		if( oRxB_on ){
			elecpowersum[12]= (double)No_coil* elecPower(ro_coil,ro_coil+t_coil_oRxB,RxBCurrent*oRxB_scale,oRxB_l*oRxB_b,copper_area);
		}
		else elecpowersum[12] = 0.;
		cout << "EL-POWER: 12. oRxB_Coil-power = " << elecpowersum[12]/(double)No_coil << " Watt, total = " << elecpowersum[12] << " Watt" << endl;

		for(int i = 0; i<13;i++){total+=elecpowersum[i];}
		cout << "EL-POWER: Total el. power consumption of NoMoS: " << total << " Watt" << endl;
	
		double currents[13] = {RxBCurrDens,outletCurrDens*outlet_scale,af_inletCurrDens*af_inlet_scale,filter_inletCurrDens*filter_inlet_scale,bf_inletCurrDens*bf_inlet_scale,helm1CurrDens*helm1_scale,filterCurrDens*filter_scale,gc_dens*gc_connec_scale,outletCorrCurrDens*outletCorr_scale,enterCurrDens*enter_scale,exitCurrDens*exit_scale,transCurrDens*trans_scale,oRxBCurrDens*oRxB_scale};
		double highestcurr=0.;
		int curr_position;
		for(int i=0;i<13;i++){
			if(highestcurr < currents[i]){highestcurr = currents[i]; curr_position = i;}
		}
		cout << "CURRDENS: Highest CurrDens. at pos: " << curr_position << ", CurrDens in coil: " << highestcurr << " A/m*m"<< endl;
		//cout << "CURRDENS: ... and CurrDens over coil (COOLING): " << highestcurr/(conduct_s*conduct_s) <<" A/m*m" <<  endl;

	}


	// in case of SC, we check the power of the NC coils
	if( Conductor == "TTSL" ){
		double copper_edge_radius = 0.00045; //edge radius in meter
		double edge_area = copper_edge_radius*copper_edge_radius*(4. - M_PI);
		double copper_area = conduct_s*conduct_s - water_dia*water_dia/4*M_PI - edge_area;

		double NC_power;
		//helm 2
		NC_power = elecPower(helm2_inner,helm2_inner+t_coil_helm2,NC_Curr,helm2_l*helm2_b,copper_area);
		NC_power += elecPower(gc_connec_inner,gc_connec_inner+t_coil_gc_connec,NC_Curr,gc_connec_l*gc_connec_b,copper_area);
		if( pg_connec_scale != 0. && pp_connec_scale != 0. ){

			NC_power += elecPower(pg_connec_inner,pg_connec_inner+t_coil_pg_connec,NC_Curr,pg_connec_l*pg_connec_b,copper_area);
			NC_power += elecPower(pp_connec_inner,pp_connec_inner+t_coil_pp_connec,NC_Curr,pp_connec_l*pp_connec_b,copper_area);
		}

		cout << "MAG: Power Consumption of NC coils is: " << NC_power << " Watts" << endl;
	}



// quickly write out file with R1 in it, so trajectory can read in, also aperture grid points for bline VS traj comparison

	ofstream configinfo;
	string infoname = filedir + "info.txt";
	configinfo.open(infoname.c_str(),ios::out);
	configinfo << "R1" << "\t" << R_1 << endl;
	double blineStartZ;

       	if( ! Perc_on ){
		if( Conductor == "NL" || Conductor == "TTSL" ){

		blineStartZ = -l_coil_RxB/2. - RxBtoIn_d - l_coil_af_inlet*2 - d_af_inlet - af_screw_gap - (d_af_inlet+l_coil_af_inlet)*((double)n_af_inlet-3.) - l_coil_af_inlet;
	       	blineStartZ = blineStartZ - FIntoAFIn_d - l_coil_filter_inlet - FIntoBFIn_d - (d_bf_inlet+l_coil_bf_inlet)*((double)n_bf_inlet-2.) - l_coil_bf_inlet - bf_screw_gap;
	      	blineStartZ = blineStartZ - l_coil_bf_inlet -helm1_bf_gap - l_coil_helm1 -  pumpport_distance/2.;
		}


		//now we overwrite this, in case we want a manual Z start
		if( boolManualStartZ ) blineStartZ = manualStartZ;

		configinfo << "blineStartZ" << "\t" << blineStartZ << endl;
	}// if on, written later, after PERC is built


//////////////////////////////////////////////////////////////////////
/////////////// inputcoil.dat creation
////////////////////////////////////////////////////////////////////////////

	string inputcoil_RxB = filedir + "inputcoil_RxB.dat";
	ofstream inputcoil_file;
	inputcoil_file.open(inputcoil_RxB.c_str(), ios::out);
	double Perc_switch;
	if( Perc_on ) Perc_switch = 1.;
	else Perc_switch = 0.;
	double totalN = n_coil_nomos*NoMoSOn + Perc_N*Perc_switch;
	inputcoil_file << totalN << "\n";
	inputcoil_file.close();
	

	double starting_point[4];
	double direction[4];
	double gate_begin;
	double connec_end;
	double filterinletstartZ = -l_coil_RxB/2. - RxBtoIn_d - l_coil_af_inlet*2 - d_af_inlet - af_screw_gap - (d_af_inlet+l_coil_af_inlet)*((double)n_af_inlet-3.) - l_coil_af_inlet;
	filterinletstartZ = filterinletstartZ - FIntoAFIn_d - l_coil_filter_inlet;
	double pump_end = filterinletstartZ - FIntoBFIn_d - (d_bf_inlet+l_coil_bf_inlet)*((double)n_bf_inlet-2.) - l_coil_bf_inlet - bf_screw_gap - l_coil_bf_inlet - helm1_bf_gap - l_coil_helm1;

	if(NoMoSOn == 1)
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
	    
		    // OUTLET coil before screw gap
	    starting_point[1] = outlet_x_shift;
	    starting_point[2] = (R_1+outlet_r_shift)*cos(alpha)-(l_coil_RxB/2.+RxBtoOut_d)*sin(alpha);
	    starting_point[3] = (R_1+outlet_r_shift)*sin(alpha)+(l_coil_RxB/2.+RxBtoOut_d)*cos(alpha);
	    direction[1] = 0.;
	    direction[2] = -sin(alpha);
	    direction[3] = cos(alpha);
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
		// input: current, starting coil number, start P, direction, length, innerR, thick, gap, n_outlet, filedir, reverse direction yes =1
	    linear_coils_to_inputcoil_dat(i_max, N_coil, starting_point, direction, l_coil_outlet, outlet_inner, t_coil_outlet, d_outlet, 1, inputcoil_RxB,0);

		    // OUTLET coils after screw gap
	    starting_point[1] = outlet_x_shift;
	    starting_point[2] = (R_1+outlet_r_shift)*cos(alpha)-(l_coil_RxB/2.+RxBtoOut_d+l_coil_outlet+outlet_screw_gap)*sin(alpha);
	    starting_point[3] = (R_1+outlet_r_shift)*sin(alpha)+(l_coil_RxB/2.+RxBtoOut_d+l_coil_outlet+outlet_screw_gap)*cos(alpha);
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
		// input: current, starting coil number, start P, direction, length, innerR, thick, gap, n_outlet, filedir, reverse direction yes =1
	    linear_coils_to_inputcoil_dat(i_max, N_coil+1, starting_point, direction, l_coil_outlet, outlet_inner, t_coil_outlet, d_outlet, n_outlet-1, inputcoil_RxB,0);

		    // OUTLET CORR coils
	    starting_point[1] = outletCorr_x_shift;
	    starting_point[2] = starting_point[2] - ((n_outlet-1)*l_coil_outlet+(n_outlet-2)*d_outlet)*sin(alpha) +(n_outletCorr*l_coil_outletCorr+(n_outletCorr-1)*d_outletCorr+outletCorr_fromoutEnd)*sin(alpha);
	    starting_point[3] = starting_point[3] + ((n_outlet-1)*l_coil_outlet+(n_outlet-2)*d_outlet)*cos(alpha) -(n_outletCorr*l_coil_outletCorr+(n_outletCorr-1)*d_outletCorr+outletCorr_fromoutEnd)*cos(alpha);
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
		// input: current, starting coil number, start P, direction, length, innerR, thick, gap, n_outlet, filedir, reverse direction yes =1
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+3, starting_point, direction, l_coil_outletCorr, outletCorr_inner, t_coil_outletCorr, d_outletCorr, n_outletCorr, inputcoil_RxB,0);
	
	
		    // INLET
		    // AF inlet before gap
	    starting_point[1] = inlet_x_shift;
	    starting_point[2] = R_1+inlet_r_shift;
	    starting_point[3] = -l_coil_RxB/2.-RxBtoIn_d-l_coil_af_inlet;
	    direction[1] = 0.;
	    direction[2] = 0.;
	    direction[3] = 1.;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet, starting_point, direction, l_coil_af_inlet, inlet_inner, t_coil_af_inlet, d_af_inlet, 2, inputcoil_RxB,1);
		//reversed direction..	
		//
		//AF inlet after gap
	    starting_point[3] = -l_coil_RxB/2.-RxBtoIn_d - l_coil_af_inlet*2 - d_af_inlet - af_screw_gap - l_coil_af_inlet;
	    direction[1] = 0.;
	    direction[2] = 0.;
	    direction[3] = 1.;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+2, starting_point, direction, l_coil_af_inlet, inlet_inner, t_coil_af_inlet, d_af_inlet, n_af_inlet-2, inputcoil_RxB,1);
	    
		// filter inlet
	    starting_point[3] = filterinletstartZ;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_af_inlet, starting_point, direction, l_coil_filter_inlet, inlet_inner, t_coil_filter_inlet, d_af_inlet, 1, inputcoil_RxB,1);

		    // BF inlet before screw gap
	    starting_point[3] = filterinletstartZ - (FIntoBFIn_d + l_coil_bf_inlet);
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_af_inlet+1, starting_point, direction, l_coil_bf_inlet, inlet_inner, t_coil_bf_inlet, d_bf_inlet, n_bf_inlet-1, inputcoil_RxB,1);
		//reversed direction..
		//after screw gap
	    starting_point[3] = pump_end + l_coil_helm1 + helm1_bf_gap;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_af_inlet+1+n_bf_inlet-1, starting_point, direction, l_coil_bf_inlet, inlet_inner, t_coil_bf_inlet, d_bf_inlet, 1, inputcoil_RxB,1);
	

	    	// HELMHOLTZ 1
	    starting_point[1] = helm_x_shift;
	    starting_point[2] = R_1+helm_r_shift;
	    starting_point[3] = pump_end;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_inlet, starting_point, direction, l_coil_helm1, helm1_inner, t_coil_helm1, 0., 1, inputcoil_RxB,0);
	
		// helm holtz 2
	    starting_point[1] = helm_x_shift;
	    starting_point[2] = R_1+helm_r_shift;
	    starting_point[3] = pump_end - pumpport_distance - l_coil_helm2;
	    connec_end = starting_point[3] - gate_dist;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_inlet +1, starting_point, direction, l_coil_helm2, helm2_inner, t_coil_helm2,0., 1, inputcoil_RxB,0);
	
	if(n_filter > 0){    
		// FILTER coils
	    starting_point[1] = filter_x_shift;
	    starting_point[2] = R_1 + filter_r_shift;
	    starting_point[3] = filterinletstartZ + l_coil_filter_inlet/2. - l_coil_filter/2.+filter_z_shift;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max, N_coil+n_outlet+n_inlet +2 , starting_point,direction, l_coil_filter, filter_inner, t_coil_filter, d_filter, n_filter, inputcoil_RxB,0);
	}
	

	    // Connector 1
	    starting_point[3] = connec_end - l_coil_gc_connec;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter,starting_point,direction,l_coil_gc_connec,gc_connec_inner,t_coil_gc_connec,0.,1,inputcoil_RxB,0);
	    
	if( pg_connec_scale != 0. && pp_connec_scale != 0. ){
	    // Connector 2
	    starting_point[3] = starting_point[3] - pipe_dist - l_coil_pg_connec;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+1,starting_point,direction,l_coil_pg_connec,pg_connec_inner,t_coil_pg_connec,0.,1,inputcoil_RxB,0);

	    // Connector 3
	    starting_point[3] = starting_point[3] - screw_dist - l_coil_pp_connec;
	    if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+2,starting_point,direction,l_coil_pp_connec,pp_connec_inner,t_coil_pp_connec,0.,1,inputcoil_RxB,0);
	}

	    cout << "MAG: End of NoMoS = " << starting_point[3] << endl;	
	    


	    // ENTER Corr
	    starting_point[1] = enter_x_shift;
	    starting_point[2] = R_1 + enter_r_shift + sin(enter_angle)*l_coil_enter/2.;
	    starting_point[3] = enter_z_shift - l_coil_RxB/2. - cos(enter_angle)*l_coil_enter/2.;
	    direction[1] = 0.;
	    direction[2] = -sin(enter_angle);
	    direction[3] = cos(enter_angle);
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+3+n_outletCorr,starting_point,direction,l_coil_enter,enter_inner,t_coil_enter,0.,1,inputcoil_RxB,0);

	    //EXIT Corr
	    starting_point[1] = exit_x_shift;
	    starting_point[2] = (R_1+outlet_r_shift)*cos(alpha)-(l_coil_RxB/2.+RxBtoOut_d-exit_z_shift)*sin(alpha);
	    starting_point[3] = (R_1+outlet_r_shift)*sin(alpha)+(l_coil_RxB/2.+RxBtoOut_d-exit_z_shift)*cos(alpha);
	    direction[1] = 0.;
	    direction[2] = -sin(alpha+exit_angle);
	    direction[3] = cos(alpha+exit_angle);
	    linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+3+n_outletCorr+1,starting_point,direction,l_coil_exit,exit_inner,t_coil_exit,0,1,inputcoil_RxB,0);
	    
	    //TRANS Corr
	    if( !trans_scale == 0. ){
	    
	        starting_point[1] = 0.;
	    	starting_point[2] = l_coil_trans/2. + l_coil_trans/2.*cos(trans_angle);
	    	starting_point[3] = R_1 + trans_r_shift + l_coil_trans/2*sin(trans_angle);
	    	direction[1] = 0.;
	    	direction[2] = -cos(trans_angle);
	    	direction[3] = -sin(trans_angle);
	    	linear_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+3+n_outletCorr+2,starting_point,direction,l_coil_trans,trans_inner,t_coil_trans,0,1,inputcoil_RxB,0);

	    }

	    // outer RxB 
	    //double offset_angle= alpha/((double)N_coil-1.)*(double)n_cloth;       // angle of the first coil
	    if(oRxB_on){
		starting_point[1] = 0.;
	    	starting_point[2] = R_2*cos(alphastart);
	    	starting_point[3] = R_2*sin(alphastart);
	    	direction[1] = 0.;
	    	direction[2] = -sin(alphastart);
	    	direction[3] = cos(alphastart);
	    	if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point, direction);}
	    	    // input: current (vector), start coil number, starting point, direction, curv. Radius (vector), inner R, thickness,coil_length, curve angle, N_coil; filedir,shiftangle
	    	circular_coils_to_inputcoil_dat(i_max,N_coil+n_outlet+n_inlet+2+n_filter+3+n_outletCorr+2+1, starting_point,direction,starting_point,ro_coil,t_coil_oRxB,l_coil_oRxB,alpha2,No_coil, inputcoil_RxB,0.);
	    }


	}


	if ( Perc_on ) {
	        // PERC coils 
		direction[1] = 0.;
		direction[2] = 0.;
		direction[3] = 1.;
		pump_end = filterinletstartZ - FIntoBFIn_d - (d_bf_inlet+l_coil_bf_inlet)*((double)n_bf_inlet-2.) - l_coil_bf_inlet - bf_screw_gap - l_coil_bf_inlet - helm1_bf_gap - l_coil_helm1;
	        connec_end = pump_end - pumpport_distance - l_coil_helm2 - gate_dist;        // position of the upstream Pump port z-position
		gate_begin = connec_end; 
	        starting_point[1] = 0.;
	        starting_point[2] = R_1+ helm_r_shift;
	        starting_point[3] = gate_begin -l_coil_gc_connec - pipe_dist- l_coil_pg_connec - screw_dist -l_coil_pp_connec - PERClastcoil_to_flansch- flansch_perc_thick;


	        if (Perc_coordinatsystem == 0) {PERC_to_normal_coordinat_transform(starting_point);}
	        PERC_coils_to_inputcoil_dat(starting_point, direction,perc_global_scale, perc_sol_scale, perc_filter_scale,perc_connec_scale, filedir, inputcoil_RxB, Perc_coordinatsystem);
                
		if( startinPERC ) {
			blineStartZ = starting_point[3] - (2.6825+0.9965/2.) -4.;
		}
		else{
			blineStartZ = manualStartZ;
			cout << "MAG: Manual Z start is set!" << endl;
		}
		configinfo << "startZ" << "\t" << blineStartZ << endl;
	}
	
	//here we write the total power consumption to the info file
	configinfo << "power" << "\t" << total << endl;
	configinfo.close();
	
	cout << "MAG: All coils built" << endl;

	///////////////////////////////////////////////////
	// Run the simulation and calculate values of the B-field
	////////////////////////////////////////////////////
	input_coils(inputcoil_RxB);
	test_coils(filedir);
	magsource(filedir);

	if(coilcenterbool){	
		coil_center(R_1, alpha, -3., 0., 1.,horilines,vertilines,ApertX,ApertY, apertYshift,apertXshift, inputcoil_RxB, filedir, 0);     // R1, alpha,startPz,RxBStart,outlet-length,horil,vertil,apX,apY,inputc,filed,prep
	}

	if(geodriftbool){
		//fieldmap(2.*R_1, fieldmap_step, inputcoil_RxB, filedir, 0);       // map size [m], step_size [m], inputdatafile , initialize magfield (0=no)
		//radial_cuts(5, N_coil, R_1*2., 5e-3, inputcoil_RxB, filedir, 0); // lines, shift[rad], line_end[m], step_size [m], inputdatafile , initialize magfield (0=no)
		//vertical_cuts(5, N_coil, 1., 5e-3, R_1, inputcoil_RxB, filedir, 0); // lines, shift[rad], line_end[m], step_size [m], inputdatafile , initialize magfield (0=no)
		geoline_anadrift(R_1,alpha,horilines,vertilines,ApertX,ApertY,apertYshift,filedir,inputcoil_RxB,0,1187.29,true); // R1, alpha, horilines,vertilines,apertsizeX,apertsizeY,..,..,prep,momentum kev,drifton/off
		geoline_anadrift(R_1,alpha,horilines,vertilines,ApertX,ApertY,apertYshift,filedir,inputcoil_RxB,0,1187.29,false); // R1, alpha, horilines,vertilines,apertsizeX,apertsizeY,..,..,prep,momentum keV,drifton/off
	}

	if(blinebool){
		b_line_real(NoMoSOn,R_1,alpha,blineStartZ,onlyCornerCenter,horilines,vertilines,ApertX,ApertY,apertYshift,apertXshift,inputcoil_RxB,filedir,0,-1,1187.29,2,detpos); //..,..,Z,horil,vertil,apertx,aperty,..,..,prep,driftdir,mom,driftOn, detposZ
		b_line_real(NoMoSOn,R_1,alpha,blineStartZ,onlyCornerCenter,horilines,vertilines,ApertX,ApertY,apertYshift,apertXshift,inputcoil_RxB,filedir,0,-1,1187.29,1,detpos); //..,..,Z,horil,vertil,apertx,aperty,..,..,prep,driftdir,momentum keV,driftOn, detposZ
		b_line_real(NoMoSOn,R_1,alpha,blineStartZ,onlyCornerCenter,horilines,vertilines,ApertX,ApertY,apertYshift,apertXshift,inputcoil_RxB,filedir,0,-1,1187.29,0,detpos); //..,..,Z,horil,vertil,apertx,aperty,..,..,prep,driftdir,mom,driftOn, detposZ
	}
	

	    // Run the simulation and calculate map the B-field, precise
	if(fieldmapbool){
		double fieldmapStartP[3];
		double fieldmapLength[3];
		if( Perc_on ) {
			fieldmapStartP[0] = -0.2;
			fieldmapStartP[1] = -0.9 - 0.2;
			fieldmapStartP[2] = -13.;
			fieldmapLength[0] = 0.4;
			fieldmapLength[1] = 2*0.9+2*0.2;
			fieldmapLength[2] = 13 + 0.9 + 0.2;
		}
		else{
			fieldmapStartP[0] = -0.2;
			fieldmapStartP[1] = -0.9 - 0.2;
			fieldmapStartP[2] = -1.7;
			fieldmapLength[0] = 0.4;
			fieldmapLength[1] = 2*0.9+2*0.2;
			fieldmapLength[2] = 1.7 + 0.9 + 0.2;
		}
		fieldmap(fieldmapStartP,fieldmapLength, fieldmap_step, inputcoil_RxB, filedir, 0);       // map size [m], step_size [m], inputdatafile , initialize magfield (0=no)
	}


	return 0;

}
