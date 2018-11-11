#include "image_inout.h"
#include "Matrix.h"
#include "Utils.h"
#include <math.h>
#include <time.h>

using namespace std;
using namespace cv;

int main() {

	time_t timeBegin, timeEnd;
	timeBegin = time(NULL);
	
	// Load source image with its features;
	cv::String img_path_s = "D:\\value_analysis\\picture\\1.jpg";
	// New a object to load the image;
	image_inout img_io_s(img_path_s); 

	// Load the img to Mat type
	cv::Mat img_s = img_io_s.image_load();//type = CV_8U_C3
	// Load the correspond image feature
	std::vector<Location_fea> img_fea_s = img_io_s.image_features_load();

	// We can show the source_image here;
	/*
	cv::namedWindow("The source image");
	cv::imshow("The source image", img_s);
	waitKey(3000);
	*/

	// Load target image with its features;
	cv::String img_path_t = "D:\\value_analysis\\picture\\2.jpg";
	// New a object to load the image;
	image_inout img_io_t(img_path_t);

	// Load the img to Mat type
	cv::Mat img_t = img_io_t.image_load();
	// Load the correspond image feature
	std::vector<Location_fea> img_fea_t = img_io_t.image_features_load();

	// We can show the target_image here;
	/*
	cv::namedWindow("The target image");
	cv::imshow("The target image", img_t);
	waitKey(3000);
	*/

	// LPS
	// Create the Matrix of L
	Matrix Mat_L_t;
	Mat_L_t.stdFea2U(img_fea_t);
	Matrix Mat_Fea_t_v;
	Mat_Fea_t_v.stdFea2stdMat_v(img_fea_t, 1);
	Matrix Mat_Fea_t_h;
	Mat_Fea_t_h.stdFea2stdMat_h(img_fea_t);
	std::vector<vector<double>> Mat_Fea_t = Mat_Fea_t_h.Mat_return(); //2*68
	vector<vector<double>> one_matrix;
	for (int i = 0; i < 3; i++) {
		vector<double> temp;
		for (int j = 0; j < 3; j++) {
			temp.push_back(0.0);
		}
		one_matrix.push_back(temp);
	}
	// Concat these matrices for L 71*71
	Mat_L_t.Mat_col_concat(Mat_Fea_t_v.Mat_return());
	Mat_Fea_t_h.Mat_col_concat(one_matrix);
	Mat_L_t.Mat_row_concat(Mat_Fea_t_h.Mat_return());

	// Create the Matrix of Y 71*2
	Matrix Mat_Fea_s_v;
	Mat_Fea_s_v.stdFea2stdMat_v(img_fea_s, 0);
	vector<vector<double>> one_matrix_32;
	for (int i = 0; i < 3; i++) {
		vector<double> temp;
		for (int j = 0; j < 2; j++) {
			temp.push_back(0.0);
		}
		one_matrix_32.push_back(temp);
	}
	Mat_Fea_s_v.Mat_row_concat(one_matrix_32);

	//-----------2018.11.11 22:37 before here: right


	// Calculate the inverse of L to get paras

	// Choose to directly solve the formula
	std::vector<double> paras_1; // 1*71
	std::vector<double> paras_2; // 1*71
	std::vector<double> feas_1; // 1*71
	std::vector<double> feas_2; // 1*71
	int paras_num = Mat_Fea_s_v.Mat_return().size();
	for (int i = 0; i < paras_num; i++) {
		feas_1.push_back(Mat_Fea_s_v.Mat_return()[i][0]);
	}
	for (int i = 0; i < paras_num; i++) {
		feas_2.push_back(Mat_Fea_s_v.Mat_return()[i][1]);
	}

	// Calculate the paras
	paras_1 = Mat_L_t.Mat_solve(feas_1);
	paras_2 = Mat_L_t.Mat_solve(feas_2);

	std::vector<vector<double>> paras;
	paras.push_back(paras_1);
	paras.push_back(paras_2);

	Matrix Mat_paras;
	Mat_paras.Mat_set(paras); //2*71

	//-----------2018.11.11 22:13 before here: right

	int height_t = img_t.size().height;
	int width_t = img_t.size().width;
	int height_s = img_s.size().height;
	int width_s = img_s.size().width;

	enum inter_kind { nearest, bilinear, bicubic }inter_input = bicubic;
	if (inter_input != nearest&&inter_input != bilinear&&inter_input != bicubic) {
		std::cout << "Please input one kind of valid way for interpolate" << endl;
		std::cout << "Possible choice: nearest, bilinear, bicubic" << endl;
		return -1;
	}

	// Multiply the paras and locs to find the correspond location in source images
	int fea_num = Mat_Fea_t[0].size();
	for (int loc_i = 0; loc_i < height_t; loc_i++) {

		for (int loc_j = 0; loc_j < width_t; loc_j++) {
			if (loc_i % 100 == 0 && loc_j % 100 == 0) {
				cout << "At " << loc_i << " and " << loc_j << endl;
			}

			// k<68 for every U(x,y)
			std::vector<vector<double>> loc;
			for (int k = 0; k < fea_num; k++) {
				double loc_x = Mat_Fea_t[1][k];
				double loc_y = Mat_Fea_t[2][k];
				double r2 = pow(loc_x - loc_i, 2) + pow(loc_y - loc_j, 2);
				double U_r2 = r2*log(r2);
				std::vector<double> temp;
				temp.push_back(U_r2);
				loc.push_back(temp);
			}
			std::vector<double> tem_1, tem_i, tem_j;
			tem_1.push_back(1.0);
			loc.push_back(tem_1);
			tem_i.push_back(double(loc_i));
			loc.push_back(tem_i);
			tem_j.push_back(double(loc_j));
			loc.push_back(tem_j);

			std::vector<vector<double>> loc_double = Mat_paras.Mat_mul(loc);
			double loc_x_double = loc_double[0][0];
			double loc_y_double = loc_double[1][0];
			// Interpolate here with three ways
			switch (inter_input) {
				case nearest: {
					int loc_x = 0;
					int loc_y = 0;
					if (loc_x_double >= 0 && loc_x_double <= height_t + 1) (loc_x_double - floor(loc_x_double)) > 0.5 ? loc_x = floor(loc_x_double) : loc_x = floor(loc_x_double) + 1;
					else loc_x = -1;
					if (loc_y_double >= 0 && loc_y_double <= width_t + 1) (loc_y_double - floor(loc_y_double)) > 0.5 ? loc_y = floor(loc_y_double) : loc_y = floor(loc_y_double) + 1;
					else loc_y = -1;

					if (loc_x > height_s - 1 || loc_x < 0 || loc_y > width_s - 1 || loc_y < 0) {
						cv::Vec3b black = cv::Vec3b(0,0,0);
						//char black[3] = {0, 0, 0};
						img_t.at<Vec3b>(loc_i, loc_j) = black;
						// value = [0,0,0]

					}
					else {
						//cv::Vec3d value = img_s.at<Vec3d>(loc_x, loc_y);
						img_t.at<Vec3b>(loc_i, loc_j) = img_s.at<Vec3b>(loc_x, loc_y);// value;
						// value = [loc_x][loc_y]

					}

					break;
				}
				case bilinear: {
					if (loc_x_double >= 0 && loc_x_double < height_s - 1 && loc_y_double >= 0 && loc_y_double < width_s - 1) {
						double loc_x_1 = floor(loc_x_double);
						double delta_x_1 = loc_x_double - loc_x_1;
						double loc_y_1 = floor(loc_y_double);
						double delta_y_1 = loc_y_double - loc_y_1;
						double loc_x_2 = floor(loc_x_double) + 1;
						double delta_x_2 = 1 - delta_x_1;
						double loc_y_2 = floor(loc_y_double) + 1;
						double delta_y_2 = 1 - delta_y_1;

						cv::Vec3b value_1 = img_s.at<Vec3b>(loc_x_1, loc_y_1);
						cv::Vec3b value_2 = img_s.at<Vec3b>(loc_x_2, loc_y_1);
						cv::Vec3b value_3 = img_s.at<Vec3b>(loc_x_1, loc_y_2);
						cv::Vec3b value_4 = img_s.at<Vec3b>(loc_x_2, loc_y_2);
						img_t.at<Vec3b>(loc_i, loc_j) = delta_x_2*delta_y_2*value_1 + delta_x_1*delta_y_2*value_2 + delta_x_2*delta_y_1*value_3 + delta_x_1*delta_y_1*value_4;

						// value = delta_x_2*delta_y_2*[loc_x_1][loc_y_1] + delta_x_1*delta_y_2*[loc_x_2][loc_y_1] + 
						// delta_x_2*delta_y_1*[loc_x_1][loc_y_2] + delta_x_1*delta_y_1*[loc_x_2][loc_y_2];
					}
					else {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						img_t.at<Vec3b>(loc_i, loc_j) = black;
						// value = [0,0,0]
					}
					break;
				}
				case bicubic: {
					// Do not interpolate at the corner of the image for lacking of confident points
					if (loc_x_double >= 1 && loc_x_double < height_s - 2 && loc_y_double >= 1 && loc_y_double < width_s -2) {
						std::vector<int> loc_x;
						for (int loc_tem = -1; loc_tem < 3; loc_tem++) {
							loc_x.push_back(int(floor(loc_x_double)) + loc_tem);
						}

						std::vector<int> loc_y;
						for (int loc_tem = -1; loc_tem < 3; loc_tem++) {
							loc_y.push_back(int(floor(loc_y_double)) + loc_tem);
						}

						std::vector<double> delta_x;
						double temp_x = loc_x_double - double(floor(loc_x_double));
						for (int loc_tem = -1; loc_tem < 3; loc_tem++) {
							delta_x.push_back(bicubic_weight(double(temp_x - loc_tem)));
						}
						
						std::vector<double> delta_y;
						double temp_y = loc_y_double - double(floor(loc_y_double));
						for (int loc_tem = -1; loc_tem < 3; loc_tem++) {
							delta_y.push_back(bicubic_weight(double(temp_y - loc_tem)));
						}

						cv::Vec3b value = Vec3b(0,0,0);
						for (int si = 0; si < 4; si++) {
							for (int sj = 0; sj < 4; sj++) {
								cv::Vec3b value_temp = img_s.at<Vec3b>(loc_x[si],loc_y[sj]);
								value += delta_x[si] * delta_y[sj] * value_temp;
								// value[i,j] = delta_x[si]*delta_y[sj]*value[loc_x[si],loc_y[sj]];
							}
						}
						img_t.at<Vec3b>(loc_i, loc_j) = value;
					}
					else {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						img_t.at<Vec3b>(loc_i, loc_j) = black;
						// value = [0,0,0]
					}
					break;
				}
				default: {
					// default use nearest
					int loc_x = 0;
					int loc_y = 0;
					if (loc_x_double >= 0) (loc_x_double - floor(loc_x_double)) > 0.5 ? loc_x = floor(loc_x_double) : loc_x = floor(loc_x_double) + 1;
					else loc_x = -1;
					if (loc_y_double >= 0) (loc_y_double - floor(loc_y_double)) > 0.5 ? loc_y = floor(loc_y_double) : loc_y = floor(loc_y_double) + 1;
					else loc_y = -1;

					if (loc_x > height_s - 1 || loc_x < 0 || loc_y > width_s - 1 || loc_y < 0) {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						img_t.at<Vec3b>(loc_i, loc_j) = black;
						// value = [0,0,0]

					}
					else {
						cv::Vec3b value = img_s.at<Vec3b>(loc_x, loc_y);
						img_t.at<Vec3b>(loc_i, loc_j) = value;
						// value = [loc_x][loc_y]

					}
					break;
				}
			}
		}
	}

	timeEnd = time(NULL);
	cout << "The total process time: " << timeEnd - timeBegin << endl;

	// We can show the altered_image here;

	cv::namedWindow("The new image");
	cv::imshow("The new image", img_t);
	waitKey(3000);

	// Save the altered image;
	img_io_s.image_save(img_t);

}
