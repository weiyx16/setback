#pragma once
#ifndef MATRIX
#define MATRIX

#include <math.h>
#include "image_inout.h"

using namespace std;
using namespace cv;

class Matrix {
public:

	void cvMat2stdMat(cv::Mat img_mat);
	void stdFea2stdMat_v(std::vector<Location_fea> fea_vector, int with_one = 1);
	void stdFea2stdMat_h(std::vector<Location_fea> fea_vector);
	void Mat_set(const std::vector<vector<double>> Mat_input);
	void stdFea2U(std::vector<Location_fea> fea_vector);

	void Mat_row_concat(std::vector<vector<double>> Mat_concat);
	void Mat_col_concat(std::vector<vector<double>> Mat_concat);
	void Mat_show();

	std::vector<double> Mat_solve(std::vector<double> Mat_b);
	
	std::vector<vector<double>> Mat_mul(std::vector<vector<double>> Mat_right);
	std::vector<vector<double>> Mat_return();

	Matrix() {};

private:
	std::vector<vector<double>> Mat_data;
	std::vector<vector<double>> Mat_data_general;
	double U_calcu(double r2);
	void change_rows(int x);
};

#endif // !MATRIX.h