#include "Utils.h"


double bicubic_weight(double distance) {
	double distance_abs = fabs(distance);
	double distance_l = 1.0;
	double distance_h = 2.0;
	if (distance_abs <= distance_l) {
		return 1.0 - 2.0*pow(distance_abs, 2.0) + pow(distance_abs, 3.0);
	}
	else {
		if (distance_l < distance_abs && distance_abs < distance_h) {
			return 4.0 - 8.0 * distance_abs + 5.0 * pow(distance_abs, 2.0) - pow(distance_abs, 3.0);
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

double B_spline_weight(double distance, int order)
{
	switch(order){
		case 0:{
			return -1.0/6.0*pow(distance, 3.0) + 1.0/2.0*pow(distance, 2.0) - 1.0/2.0*distance + 1.0/6.0;
			break;
		}
		case 1:{
			return 1.0/2.0*pow(distance, 3.0) - pow(distance, 2.0) + 2.0/3.0;
			break;
		}
		case 2:{
			return -1.0/2.0*pow(distance, 3.0) + 1.0/2.0*pow(distance, 2.0) + 1.0/2.0*distance + 1.0/6.0;
			break;
		}
		case 3:{
			return 1.0/6.0*pow(distance, 3.0);
			break;
		}
		default:{
			return 0.0;
			break;
		}
	}
}

void grid_init(std::vector<double> D_locx, std::vector<double> S_locx, std::vector<double> D_locy, std::vector<double> S_locy, cv::Mat & Graph_tran_grid_x, cv::Mat & Graph_tran_grid_y, 
	int grid_size, double max_s_locx, double min_s_locx, double max_s_locy, double min_s_locy)
{
	// initial the morphing grid

	int grid_height = Graph_tran_grid_x.size().height;
	int grid_width = Graph_tran_grid_x.size().width;

	int min_grid_x = floor(min_s_locx / double(grid_size)) - 1;
	min_grid_x = min_grid_x > 0 ? min_grid_x : 0;
	int max_grid_x = floor(max_s_locx / double(grid_size)) + 2;
	max_grid_x = max_grid_x > grid_height ? grid_height : max_grid_x;

	int min_grid_y = floor(min_s_locy / double(grid_size)) - 1;
	min_grid_y = min_grid_y > 0 ? min_grid_y : 0;
	int max_grid_y = floor(max_s_locy / double(grid_size)) + 2;
	max_grid_y = max_grid_y > grid_height ? grid_height : max_grid_y;

	for (int i = min_grid_x; i < max_grid_x; i++) {
		for (int j = min_grid_y; j < max_grid_y; j++) {
			std::vector<int> local_fea; // save feature index which is near ij
			for (int fea = 0; fea < S_locx.size(); fea++) {
				// the feature point should in the two grid distance in i,j and get the index here
				if (abs(S_locx[fea] / double(grid_size) - i) < 2 && abs(S_locy[fea] / double(grid_size) - j) < 2) {
					local_fea.push_back(fea);
				}
			}
			if (local_fea.size() == 0) {
			}
			else {
				std::vector<double> phi_c_x;
				std::vector<double> phi_c_y;
				std::vector<double> weight_c2;
				for (int fea_local = 0; fea_local < local_fea.size(); fea_local++) {
					int index = local_fea[fea_local];
					double s = S_locx[index] / double(grid_size) - floor(S_locx[index] / double(grid_size));
					double t = S_locy[index] / double(grid_size) - floor(S_locy[index] / double(grid_size));
					int k = (i + 1) - floor(S_locx[index] / double(grid_size)); // 0-3
					int l = (j + 1) - floor(S_locy[index] / double(grid_size));
					double wc = B_spline_weight(s, k)*B_spline_weight(t, l);
					weight_c2.push_back(pow(wc, 2.0));
					double wab = 0.0;
					for (int a = 0; a < 4; a++) {
						for (int b = 0; b < 4; b++) {
							wab += pow(B_spline_weight(s, a)*B_spline_weight(t, b), 2.0);
						}
					}
					phi_c_x.push_back(wc*D_locx[index] / wab);
					phi_c_y.push_back(wc*D_locy[index] / wab);
				}
				double list_x = 0.0;
				double list_y = 0.0;
				double wc2_sum = 0.0;
				for (int fea_local = 0; fea_local < local_fea.size(); fea_local++) {
					wc2_sum += weight_c2[fea_local];
					list_x += weight_c2[fea_local] * phi_c_x[fea_local];
					list_y += weight_c2[fea_local] * phi_c_y[fea_local];
				}
				list_x /= wc2_sum;
				list_y /= wc2_sum;
				Graph_tran_grid_x.at<double>(i, j) += list_x;
				Graph_tran_grid_y.at<double>(i, j) += list_y;
			}
		}
	}
}

