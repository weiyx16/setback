#pragma once
#ifndef Utils
#define Utils

#include <math.h>
#include <vector>
#include "image_inout.h"
using namespace std;
using namespace cv;

//weight for bicubic intepolate
double bicubic_weight(double distance);
//U funtion
double U_calcu(double r2);
//Use for B_spline
double B_spline_weight(double distance, int order);
void grid_init(std::vector<double> D_locx, std::vector<double> S_locx, std::vector<double> D_locy, std::vector<double> S_locy, 
	cv::Mat &Graph_tran_grid_x, cv::Mat &Graph_tran_grid_y, int grid_size, double max_s_locx, double min_s_locx, double max_s_locy, double min_s_locy);
#endif // !Utils.h