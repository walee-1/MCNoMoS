double elecPower(double innerR, double outerR, double current, int turns, double copper_area){

	double conductor_length = (innerR+outerR)*M_PI*turns;
	double copper_res = 1.989e-8; // Ohm * m
	double resistance = copper_res*conductor_length/copper_area;
	double power = resistance*current*current;
	
	return power;
}
