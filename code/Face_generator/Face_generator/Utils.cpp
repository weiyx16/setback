#include "Utils.h"


double bicubic_weight(double distance) {
	double distance_abs = fabs(distance);
	double distance_l = 1.0;
	double distance_h = 2.0;
	if (distance_abs <= distance_l) {
		return 1.0 - 2.0*pow(distance_abs, 2.0) + pow(distance_abs, 3.0);
	}
	else {
		if (distance_l < distance_abs < distance_h) {
			return 4.0 - 8 * distance_abs + 5 * pow(distance_abs, 2.0) - pow(distance_abs, 3.0);
		}
		else return 0.0;
	}
}

double U_calcu(double r2)
{
	if (fabs(r2) < 0.0001) {
		return 0.0;
	}
	else {
		return (r2*log(r2));
	}
}