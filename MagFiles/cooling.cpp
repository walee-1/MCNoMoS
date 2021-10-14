double iterations(double delta_p, double meanTemp, double water_dia, double conduct_length, double& Rey, double& resist){

	double wDens = 999.8; //kg/m³
	double visco = (-0.0028*pow(meanTemp,3) + 0.5968*pow(meanTemp,2) - 46.822*meanTemp + 1733.6)/(1.e6);
	cout << "COOLING: VISCO = " << visco << endl;
	double v = delta_p*water_dia*water_dia/(32.*visco*conduct_length);

	Rey = v*wDens*water_dia/visco;
	
	if ( Rey <= 2100. )
		{resist = 64./Rey;}
	else if (Rey <= 3000.)
		{resist = (64/Rey+3.164/sqrt(Rey))/2;}
	else if (Rey <= 100000) 
		{resist = 3.164/sqrt(Rey);}
	else {resist = 0.0054 + 0.3964/pow(Rey,0.3);}
	
	double diff = 1.;
	double vold;
	int counter=0;

	while( (fabs(diff) > 0.01) ){
		vold = v;
		v = (sqrt(water_dia * delta_p *2. / (resist * conduct_length * 998.))+ v)/2.;
		
		Rey = v*wDens*water_dia/visco;
		
		if ( Rey <= 2100. )
			{resist = 64./Rey;}
		else if (Rey <= 3000.)
			{resist = (64/Rey+3.164/sqrt(Rey))/2;}
		else if (Rey <= 100000) 
			{resist = 3.164/sqrt(Rey);}
		else {resist = 0.0054 + 0.3964/pow(Rey,0.3);}
		
		diff = (v - vold)/v;
		counter ++;
	}
		
	cout << "ITERATOR: Iterations = " << counter << endl;

	return v;
}


int cooling(double innerR, double realthick, double realwidth, int turns, int N_coil, double water_dia,double copper_area){

	double copper_res = 1.989e-8; // Ohm * m
	double wDens = 999.8; //kg/m³

	double cond_length = (innerR + innerR+realthick)*M_PI * turns;
	cout << endl << "COOLING: Conductor length = " << cond_length << endl;
	double deltaP= 6.5e5;
	double Tin = 20.;
	double dT = 50.;
	double Rey, resist;

	double velo = iterations(deltaP,Tin+dT/2, water_dia, cond_length,Rey, resist);
	cout << "COOLING: velo = " << velo << " m/s" <<  endl;
	cout << "COOLING: Rey = " << Rey << endl;
	double flow = velo * water_dia*water_dia /4 *M_PI * wDens; //in kg/s

	double specW = 4180.5;

	double power = dT * flow * specW;
	cout << "COOLING: Cooling power per coil = " << power << " W" << endl;
	cout << "COOLING: Cooling power in RxB = " << power*N_coil << " W" << endl;

	double elRes = cond_length/copper_area*copper_res;
	double current = sqrt(power/elRes);
	cout << "COOLING: Current from Cooling = " << current << " A" <<  endl;
	cout << "COOLING: CurrDens in copper = " << current/copper_area*1.0e-6 << " A/mm*mm" << endl;
	cout << "COOLING: CurrDens over whole coil = " << current/(realthick*realwidth/turns)*1.0e-6 << " A/mm*mm" << endl;

	return 0;

}
