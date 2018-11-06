#include "Matrix.h"

using namespace std;
using namespace cv;

//TODO
void Matrix::cvMat2stdMat(cv::Mat img_mat)
{
	
}

// 68*3 or 68*2
void Matrix::stdFea2stdMat_v(std::vector<Location_fea> fea_vector, int with_one = 1)
{
	int point_num = fea_vector.size();
	if (with_one == 1){
		for (int i=0; i<point_num; i++){
			std::vector<double> temp;
			Location_fea temp_fea = fea_vector[i];
			temp.push_back(1.0);
			temp.push_back(temp_fea.loc_x);
			temp.push_back(temp_fea.loc_y);
			Mat_data.push_back(temp);
		}
	}
	else{
		for (int i=0; i<point_num; i++){
			std::vector<double> temp;
			Location_fea temp_fea = fea_vector[i];
			temp.push_back(temp_fea.loc_x);
			temp.push_back(temp_fea.loc_y);
			Mat_data.push_back(temp);
		}
	}
}

// 3*68
void Matrix::stdFea2stdMat_h(std::vector<Location_fea> fea_vector)
{
	int point_num = fea_vector.size();
	std::vector<double> vector_x;
	std::vector<double> vector_y;
	std::vector<double> vector_one;
	for (int i=0; i<point_num; i++){
		Location_fea temp_fea = fea_vector[i];
		vector_one.push_back(1.0);
		vector_x.push_back(temp_fea.loc_x);
		vector_y.push_back(temp_fea.loc_y);
	}
	Mat_data.push_back(vector_one);
	Mat_data.push_back(vector_x);
	Mat_data.push_back(vector_y);
}

void Matrix::Mat_set(const std::vector<vector<double>> Mat_input)
{
	Mat_data = Mat_input;
}

//TODO
void Matrix::Mat_inv()
{ }

// concatenation in the row
void Matrix::Mat_row_concat(std::vector<vector<double>> Mat_concat)
{
	int row_num = Mat_concat.size();
	for (int i=0; i<row_num; i++){
		Mat_data.push_back(Mat_concat[i]);
	}
}

// concatenation in the column
void Matrix::Mat_col_concat(std::vector<vector<double>> Mat_concat)
{
	int row_num = Mat_data.size();
	int col_num = Mat_concat[0].size();
	for (int i=0; i<row_num; i++){
		for (int j=0; j<col_num; j++){
			Mat_data[i].push_back(Mat_concat[i][j]);
		}
	}
}

//TODO
std::vector<double> Matrix::Mat_solve(std::vector<double> Mat_b)
{
	double m = 0;
	double flag = 0;
	int rows = Mat_data.size();
	int cols = rows + 1;

	Mat_data_general = Mat_data;

	// Create the Generalized matrix
	for (int i = 0; i < rows; i++) {
		Mat_data_general[i].push_back(Mat_b[i]);
	}

	for (int i = 0; i < rows; i++) {
		change_rows(i);
		if (Mat_data_general[i][i] != 0) {
			for (int j = i + 1; j<rows; j++) {
				//m = 1;
				if (Mat_data_general[j][i] == 0)
					continue;
				m = abs(1.0 * Mat_data_general[i][i] / Mat_data_general[j][i]);

				if ((Mat_data_general[j][i] > 0 && Mat_data_general[i][i] > 0) || (Mat_data_general[j][i] < 0 && Mat_data_general[i][i] < 0)) {
					flag = 1;
				}
				else{
					flag = -1;
				}
				for (int k = i; k<cols; k++) {
					Mat_data_general[j][k] = Mat_data_general[j][k] * m - flag*Mat_data_general[i][k];
				}
			}
		}
	}

	std::vector<double> ans;
	double sum = 0;

	for (int i = 1; i < cols; i++) {
		sum = Mat_data_general[cols - 1 - i][cols - 1];
		for (int j = 1; j < i; j++) {
			sum -= ans[cols - j] * Mat_data_general[rows - i][rows - j];
		}
		ans[cols - i] = sum / Mat_data_general[cols - 1 - i][cols - 1 - i];
	}
	return ans;
}

void Matrix::change_rows(int x) {
	double temp = 0;
	int rows = Mat_data_general.size();
	int cols = rows + 1;
	double max = abs(Mat_data_general[x][x]);
	int loc = x;
	for (int i = x; i < rows; i++) {
		if (abs(Mat_data_general[i][x]) > max) {
			max = abs(Mat_data_general[i][x]);
			loc = i;
		}
	}
	if (x != loc) {
		for (int i = x; i < cols; i++) {
			temp = Mat_data_general[x][i];
			Mat_data_general[x][i] = Mat_data_general[loc][i];
			Mat_data_general[loc][i] = temp;
		}
	}
}

//TODO
std::vector<vector<double>> Matrix::Mat_mul(std::vector<vector<double>> Mat_right)
{
	return std::vector<vector<double>>();
}

std::vector<vector<double>> Matrix::Mat_return()
{
	return Mat_data;
}

double Matrix::U_calcu(double r2)
{
	if (abs(r2) < 0.00001) {
		return 0.0;
	}
	else {
		return (r2*log(r2));
	}
}
	
void Matrix::stdFea2U(std::vector<Location_fea> fea_vector)
{
	int point_num = fea_vector.size();
	std::vector<vector<double>> fea_temp;
	for (int i=0; i<point_num; i++){
		std::vector<double> temp;
		Location_fea fea = fea_vector[i];
		temp.push_back(fea.loc_x);
		temp.push_back(fea.loc_y);
		fea_temp.push_back(temp);
	}
	
	for (int i=0; i<point_num; i++){
		std::vector<double> temp;
		double x_i = fea_temp[i][0];
		double y_i = fea_temp[i][1];
		for (int j=0; j<point_num; j++){
			double x_j = fea_temp[j][0];
			double y_j = fea_temp[j][1];
			double U = U_calcu(pow(x_i- x_j,2)+pow(y_i - y_j,2));
			temp.push_back(U);
		}
		Mat_data.push_back(temp);
	}		
	
}
