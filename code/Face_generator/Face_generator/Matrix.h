#pragma once
#ifndef MATRIX
#define MATRIX

#include <math.h>
#include "image_inout.h"
#include "Utils.h"

using namespace std;
using namespace cv;

class Matrix {
public:
	// feature structure to myMat vertical
	void stdFea2stdMat_v(std::vector<Location_fea> fea_vector, int with_one = 1);
	// feature structure to myMat horizonal
	void stdFea2stdMat_h(std::vector<Location_fea> fea_vector);
	void Mat_set(const std::vector<vector<double>> Mat_input);
	// Create a U matrix from feature structure
	void stdFea2U(std::vector<Location_fea> fea_vector);

	// concat the mat in row
	void Mat_row_concat(std::vector<vector<double>> Mat_concat);
	// concat the mat in column
	void Mat_col_concat(std::vector<vector<double>> Mat_concat);
	void Mat_show();

	// solve the AX=B
	std::vector<double> Mat_solve(std::vector<double> Mat_b);
	// Mat multiply
	std::vector<vector<double>> Mat_mul(std::vector<vector<double>> Mat_right);
	std::vector<vector<double>> Mat_return();

	Matrix() {};

private:
	std::vector<vector<double>> Mat_data;
	std::vector<vector<double>> Mat_data_general;
	void change_rows(int x);
};

#endif // !MATRIX.h