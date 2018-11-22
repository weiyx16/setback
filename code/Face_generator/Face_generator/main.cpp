#include "image_inout.h"
#include "Matrix.h"
#include "Utils.h"
#include <math.h>
#include <time.h>

using namespace std;
using namespace cv;

int main(int argc, char* argv[]) {
	/*Output the help info*/
	
	if (argc == 2) {
		std::string help = argv[1];
		std::string help_ = "-h";
		if (!help.compare(help_)) {
			std::cout << ">> Help infomation belows: " << endl;
			std::cout << ">> Please input four parameters." << endl;
			std::cout << ">> Including source image path, target image path, interpolate methods and warp methods" << endl;
			std::cout << ">> Use absolute path please" << endl;
			std::cout << ">> Possible interpolate choice: nearest, bilinear, bicubic" << endl;
			std::cout << ">> Possible warp choice: TPS, BS, BSNaive" << endl;
			std::cout << ">> NOTICE: If you want to use BSNaive, use nearest only" << endl;
			std::cout << ">> E.g. D:\\1.jpg D:\\2.jpg bicubic TPS" << endl;
			return -1;
		}
		else {
			std::cout << ">> Unvalid input! Try -h please" << endl;
			return -1;
		}
	}
	
	/*Read parameters from command*/
	
	if (argc != 5) {
		std::cout << ">> Please input enough parameters" << endl;
		return -1;
	}
	
	
	std::string img_source_path = argv[1];
	std::string img_target_path = argv[2];
	std::string interpolate = argv[3];
	std::string warping = argv[4];
	
	/*Interpolate method and wrap method choice*/
	enum inter_kind { nearest, bilinear, bicubic }inter_input = bicubic;
	
	std::string n = "nearest";
	std::string l = "bilinear";
	std::string c = "bicubic";
	if (!interpolate.compare(n)) {
		inter_input = nearest;
	}
	else {
		if (!interpolate.compare(l)) {
			inter_input = bilinear;
		}
		else {
			if (!interpolate.compare(c)) {
				inter_input = bicubic;
			}
			else {
				std::cout << ">> Please input one kind of valid way for interpolate" << endl;
				std::cout << ">> Possible choice: nearest, bilinear, bicubic" << endl;
				return -1;
			}
		}
	}
	
	enum warp_kind { TPS, BS, BSNaive }warp_input = TPS;
	
	std::string t = "TPS";
	std::string b = "BS";
	std::string bn = "BSNaive";
	if (!warping.compare(t)) {
		warp_input = TPS;
	}
	else {
		if (!warping.compare(b)) {
			warp_input = BS;
		}
		else {
			if (!warping.compare(bn)) {
				warp_input = BSNaive;
			}
			else {
				std::cout << ">> Please input one kind of valid way for warping!" << endl;
				std::cout << ">> Possible choice: TPS, BS, BSNaive" << endl;
				return -1;
			}
		}
	}
	
	/*Image and features input*/
	time_t timeBegin, timeEnd;
	timeBegin = time(NULL);
	
	// Load source image with its features;
	cv::String img_path_s = img_source_path;//"D:\\ATsinghua\\value_analysis\\picture\\1.jpg"; //
	// New a object to load the image;
	image_inout img_io_s(img_path_s); 

	// Load the img to Mat type
	cv::Mat img_s = img_io_s.image_load();//type = CV_8U_C3
	// Load the correspond image feature
	std::vector<Location_fea> img_fea_s = img_io_s.image_features_load();

	// We can show the source_image here;
	
	cv::namedWindow("The source image");
	cv::imshow("The source image", img_s);
	cv::waitKey(3000);
	cv::destroyAllWindows();

	// Load target image with its features;
	cv::String img_path_t = img_target_path; // "D:\\ATsinghua\\value_analysis\\picture\\2.jpg"; //
	// New a object to load the image;
	image_inout img_io_t(img_path_t);

	// Load the img to Mat type
	cv::Mat img_t = img_io_t.image_load();
	// Load the correspond image feature
	std::vector<Location_fea> img_fea_t = img_io_t.image_features_load();

	// We can show the target_image here;
	
	cv::namedWindow("The target image");
	cv::imshow("The target image", img_t);
	cv::waitKey(3000);
	cv::destroyAllWindows();

	switch (warp_input)
	{
	case TPS: {
		// LPS
		// Create the Matrix of L
		Matrix Mat_L_t;
		Mat_L_t.stdFea2U(img_fea_t);
		Matrix Mat_Fea_t_v;
		Mat_Fea_t_v.stdFea2stdMat_v(img_fea_t, 1);
		Matrix Mat_Fea_t_h;
		Mat_Fea_t_h.stdFea2stdMat_h(img_fea_t);
		std::vector<vector<double>> Mat_Fea_t = Mat_Fea_t_h.Mat_return(); //3*68 take [1][2]row
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
		std::vector<vector<double>> Mat_Fea_s = Mat_Fea_s_v.Mat_return(); // 68*2 
		vector<vector<double>> one_matrix_32;
		for (int i = 0; i < 3; i++) {
			vector<double> temp;
			for (int j = 0; j < 2; j++) {
				temp.push_back(0.0);
			}
			one_matrix_32.push_back(temp);
		}
		Mat_Fea_s_v.Mat_row_concat(one_matrix_32);

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

		int height_t = img_t.size().height;
		int width_t = img_t.size().width;
		int height_s = img_s.size().height;
		int width_s = img_s.size().width;

		// Multiply the paras and locs to find the correspond location in source images
		int fea_num = Mat_Fea_t[0].size();
		for (int loc_i = 0; loc_i < height_t; loc_i++) {
			for (int loc_j = 0; loc_j < width_t; loc_j++) {
				// k<68 for every U(x,y)
				std::vector<vector<double>> loc;
				for (int k = 0; k < fea_num; k++) {
					double loc_x = Mat_Fea_t[1][k];
					double loc_y = Mat_Fea_t[2][k];
					double r2 = pow(loc_x - loc_i, 2) + pow(loc_y - loc_j, 2);
					double U_r2 = U_calcu(r2);
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
					if (loc_x_double >= 0 && loc_x_double <= height_s + 1) (loc_x_double - floor(loc_x_double)) > 0.5 ? loc_x = floor(loc_x_double) + 1 : loc_x = floor(loc_x_double);
					else loc_x = -1;
					if (loc_y_double >= 0 && loc_y_double <= width_s + 1) (loc_y_double - floor(loc_y_double)) > 0.5 ? loc_y = floor(loc_y_double) + 1 : loc_y = floor(loc_y_double);
					else loc_y = -1;

					if (loc_x > height_s - 1 || loc_x < 0 || loc_y > width_s - 1 || loc_y < 0) {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						img_t.at<Vec3b>(loc_i, loc_j) = black;
					}
					else {
						img_t.at<Vec3b>(loc_i, loc_j) = img_s.at<Vec3b>(loc_x, loc_y);
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
						cv::Vec3b value = delta_x_2*delta_y_2*value_1 + delta_x_1*delta_y_2*value_2 + delta_x_2*delta_y_1*value_3 + delta_x_1*delta_y_1*value_4;
						if (value[0] < 0.0 || value[0] > 255.0) value[0] = 0;
						else value[0] = int(value[0]);
						if (value[1] < 0.0 || value[1] > 255.0) value[1] = 0;
						else value[1] = int(value[1]);
						if (value[2] < 0.0 || value[2] > 255.0) value[2] = 0;
						else value[2] = int(value[2]);
						img_t.at<Vec3b>(loc_i, loc_j) = value;
					}
					else {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						img_t.at<Vec3b>(loc_i, loc_j) = black;
					}
					break;
				}
				case bicubic: {
					// Do not interpolate at the corner of the image for lacking of confident points
					if (loc_x_double >= 1 && loc_x_double < height_s - 2 && loc_y_double >= 1 && loc_y_double < width_s - 2) {
						cv::Vec3b value;
						int i = floor(loc_x_double);
						int j = floor(loc_y_double);
						double u = loc_x_double - i;
						double v = loc_y_double - j;
						double SU[4] = { bicubic_weight(u + 1),bicubic_weight(u),bicubic_weight(u - 1),bicubic_weight(u - 2) };
						cv::Mat A(1, 4, CV_64FC1, SU);
						double SV[4] = { bicubic_weight(v + 1),bicubic_weight(v),bicubic_weight(v - 1),bicubic_weight(v - 2) };
						cv::Mat C(4, 1, CV_64FC1, SV);

						// For mat multiply convience, deal with single color one time
						Mat B_single = Mat::zeros(4, 4, CV_64FC1);
						for (int channel = 0; channel < 3; channel++) {
							for (int m = -1; m < 3; m++) {
								for (int n = -1; n < 3; n++) {
									B_single.at<double>(m + 1, n + 1) = img_s.at<Vec3b>(i + m, j + n)[channel]; //  < 255 ? img_s.at<Vec3b>(i + m, j + n)[channel] : 255;
								}
							}
							Mat single_color = (A*B_single*C);
							double color = single_color.at<double>(0);
							if (color < 0.0 || color > 255.0) color = 0;
							else (color - floor(color)) > 0.5 ? color = floor(color) + 1 : color = floor(color);
							value[channel] = color;
						}
						img_t.at<Vec3b>(loc_i, loc_j) = value;
					}
					else {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						img_t.at<Vec3b>(loc_i, loc_j) = black;
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
					}
					break;
				}
				}
			}
		}

		timeEnd = time(NULL);
		cout << "The total process time: " << timeEnd - timeBegin << endl;

		// We can show the altered_image here;

		cv::namedWindow("The altered image");
		cv::imshow("The altered image", img_t);
		cv::waitKey(3000);
		cv::destroyAllWindows();

		// Save the altered image;
		img_io_s.image_save(img_t, img_path_t);
		break;
	}
	case BS: {	
		// B spline here
		
		// normalize the feature points
		// Find the boundary of the face area
		// struct feature to matrix first
		Matrix Mat_Fea_t_h;
		Mat_Fea_t_h.stdFea2stdMat_h(img_fea_t);
		std::vector<vector<double>> Mat_Fea_t = Mat_Fea_t_h.Mat_return(); //3*68 take [1][2]row
		Matrix Mat_Fea_s_v;
		Mat_Fea_s_v.stdFea2stdMat_v(img_fea_s, 0);
		std::vector<vector<double>> Mat_Fea_s = Mat_Fea_s_v.Mat_return(); // 68*2

		std::vector<double> T_locx = Mat_Fea_t[1]; //1*68
		std::vector<double> T_locy = Mat_Fea_t[2];
		double max_t_locx = *std::max_element(std::begin(T_locx), std::end(T_locx));
		double max_t_locy = *std::max_element(std::begin(T_locy), std::end(T_locy));
		double min_t_locx = *std::min_element(std::begin(T_locx), std::end(T_locx));
		double min_t_locy = *std::min_element(std::begin(T_locy), std::end(T_locy));

		std::vector<double> S_locx;
		std::vector<double> S_locy;
		for (int i = 0; i < Mat_Fea_s.size(); i++) {
			S_locx.push_back(Mat_Fea_s[i][0]);
			S_locy.push_back(Mat_Fea_s[i][1]);
		}
		double max_s_locx = *std::max_element(std::begin(S_locx), std::end(S_locx));
		double max_s_locy = *std::max_element(std::begin(S_locy), std::end(S_locy));
		double min_s_locx = *std::min_element(std::begin(S_locx), std::end(S_locx));
		double min_s_locy = *std::min_element(std::begin(S_locy), std::end(S_locy));

		// Normalize the target area to source area.
		// Only modify the relative location
		// I think it is better
		double scale_x = (max_s_locx - min_s_locx) / (max_t_locx - min_t_locx);
		double tran_x = (max_t_locx*min_s_locx - max_s_locx*min_t_locx) / (max_t_locx - min_t_locx);
		double scale_y = (max_s_locy - min_s_locy) / (max_t_locy - min_t_locy);
		double tran_y = (max_t_locy*min_s_locy - max_s_locy*min_t_locy) / (max_t_locy - min_t_locy);
		for (int i = 0; i < T_locx.size(); i++) {
			double temp_x = T_locx[i];
			double temp_y = T_locy[i];
			temp_x = temp_x * scale_x + tran_x;
			temp_y = temp_y * scale_y + tran_y;
			T_locx[i] = temp_x;
			T_locy[i] = temp_y;
		}

		/*
			lets show the modified feature point on the image to compare please.
		*/
		// Notice the difference in the coordinates!!!
		cv::Mat img_feature = img_s.clone();
		IplImage img_feature_arr = img_feature;
		for (int i = 0; i < S_locx.size(); i++){
			CvPoint pointplot_T;
			pointplot_T.y=int(T_locx[i]); // x and y is just different here...
			pointplot_T.x=int(T_locy[i]);
			cvCircle(&img_feature_arr, pointplot_T ,3 , cv::Scalar(0,255,0));//, 1, 8, 3 );
			CvPoint pointplot_S;
			pointplot_S.y=int(S_locx[i]);
			pointplot_S.x=int(S_locy[i]);
			cvCircle(&img_feature_arr, pointplot_S ,3 , cv::Scalar(255,0,0));//, 1, 8, 3 );
		}
		cv::Mat feature_sav = cvarrToMat(&img_feature_arr);
		cv::namedWindow("The feature image");
		cv::imshow("The feature image", feature_sav);
		cv::imwrite("D:\\ATsinghua\\value_analysis\\features.jpg", feature_sav);
		cvWaitKey(1000);
		
		Matrix Mat_Fea_t_scale;
		std::vector<vector<double>> Fea_t_scale; // 2*68
		Fea_t_scale.push_back(T_locx);
		Fea_t_scale.push_back(T_locy);
		Mat_Fea_t_scale.Mat_set(Fea_t_scale);

		// Create two grid translation graph (X Y direction divided)
		int height_s = img_s.size().height;
		int width_s = img_s.size().width;
		int height_t = img_t.size().height;
		int width_t = img_t.size().width;
		std::vector<double> D_locx;
		std::vector<double> D_locy;
		double max_derative = 0; 
		for (int i = 0; i < S_locx.size(); i++) {
			D_locx.push_back(T_locx[i] - S_locx[i]);
			D_locy.push_back(T_locy[i] - S_locy[i]);
			double derative = sqrt(pow((T_locx[i] - S_locx[i]), 2) + pow((T_locy[i] - S_locy[i]), 2));
			if (derative > max_derative) {
				max_derative = derative;
			}
		}
		
		int grid_size = floor(max_derative);
		cout << "Grid size: " << grid_size << endl;
		cv::Mat Graph_tran_grid_x = cv::Mat::zeros(floor(height_s / grid_size), floor(width_s / grid_size), CV_64FC1);
		cv::Mat Graph_tran_grid_y = cv::Mat::zeros(floor(height_s / grid_size), floor(width_s / grid_size), CV_64FC1);
		int grid_height = Graph_tran_grid_x.size().height;
		int grid_width = Graph_tran_grid_x.size().width;

		// Only initialize the grid to solve the grid and plot the calculated feature points to prove the right of grid

		grid_init(D_locx, S_locx, D_locy, S_locy, Graph_tran_grid_x, Graph_tran_grid_y, grid_size, max_s_locx, min_s_locx, max_s_locy, min_s_locy);

		// Create two pixel translation graph (X Y direction divided) from grid

		cv::Mat Graph_tran_x = cv::Mat::zeros(height_s, width_s, CV_64FC1);
		cv::Mat Graph_tran_y = cv::Mat::zeros(height_s, width_s, CV_64FC1);
		for (int i = 0; i < height_s; i++) {
			for (int j = 0; j < width_s; j++) {
				// grid: x-1 ~ x+2 y-1 ~ y+2
				// if grid = 1, then the u = x/grid - floor(x/grid) = 0 then the weight is const number
				double u = double(i) / grid_size - floor(double(i) / grid_size);
				double v = double(j) / grid_size - floor(double(j) / grid_size);
				double grid_i = floor(double(i) / grid_size);
				double grid_j = floor(double(j) / grid_size);
				double x_list = 0.0;
				double y_list = 0.0;
				for (int near_i = - 1; near_i < 3; near_i++) {
					for (int near_j = - 1; near_j < 3; near_j++) {
						if (grid_i + near_i < 0 || grid_i + near_i > grid_height - 1 || grid_j + near_j < 0 || grid_j + near_j > grid_width - 1) {
							x_list += 0.0;
							y_list += 0.0;
						}
						else {
							x_list += B_spline_weight(u, near_i + 1)*B_spline_weight(v, near_j + 1)*Graph_tran_grid_x.at<double>(grid_i + near_i, grid_j + near_j);
							y_list += B_spline_weight(u, near_i + 1)*B_spline_weight(v, near_j + 1)*Graph_tran_grid_y.at<double>(grid_i + near_i, grid_j + near_j);
						}					
					}
				}
				//x_list = fabs(x_list) > double(2.0*grid_size) ? x_list / fabs(x_list) * grid_size : x_list;
				//y_list = fabs(y_list) > double(2.0*grid_size) ? y_list / fabs(y_list) * grid_size : y_list;
				Graph_tran_x.at<double>(i, j) = x_list;
				Graph_tran_y.at<double>(i, j) = y_list;
			}
		}

		// Alter the image
		cv::Mat img_output = cv::Mat::zeros(height_s, width_s, CV_8UC3); //cv::Mat::zeros(height_s, width_s, CV_8UC3);
		for (int loc_i = 0; loc_i < height_s; loc_i++) {
			for (int loc_j = 0; loc_j < width_s; loc_j++) {
				// i j is x' y' -> x,y

				double x_list = Graph_tran_x.at<double>(loc_i, loc_j);
				double y_list = Graph_tran_y.at<double>(loc_i, loc_j);
				/*
				double loc_x_double = loc_i - x_list;
				double loc_y_double = loc_j - y_list;
				*/

				/*
				Another way
				*/

				double x_ = loc_i - x_list; // double x_ = loc_i - x_list;
				double y_ = loc_j - y_list; // double y_ = loc_j - y_list; maybe this make more sense
				double i_ = floor(x_);
				double u_ = x_ - i_;
				double j_ = floor(y_);
				double v_ = y_ - j_;

				double x_list_ = 0.0;
				double y_list_ = 0.0;
				for (int near_i = -1; near_i < 3; near_i++) {
					for (int near_j = -1; near_j < 3; near_j++) {
						if (i_ + near_i < 0 || i_ + near_i > height_s - 1 || j_ + near_j < 0 || j_ + near_j > width_s - 1) {
							x_list_ += 0.0;
							y_list_ += 0.0;
						}
						else {
							x_list_ += B_spline_weight(u_, near_i + 1)*B_spline_weight(v_, near_j + 1)
								*Graph_tran_x.at<double>(i_ + near_i, j_ + near_j);
							y_list_ += B_spline_weight(u_, near_i + 1)*B_spline_weight(v_, near_j + 1)
								*Graph_tran_y.at<double>(i_ + near_i, j_ + near_j);
						}
					}
				}

				double loc_x_double = loc_i - x_list_;
				double loc_y_double = loc_j - y_list_;

				// Interpolate here with three ways
				switch (inter_input) {
				case nearest: {
					int loc_x = 0;
					int loc_y = 0;
					if (loc_x_double >= 0 && loc_x_double <= height_s + 1) (loc_x_double - floor(loc_x_double)) > 0.5 ? loc_x = floor(loc_x_double) + 1 : loc_x = floor(loc_x_double);
					else loc_x = -1;
					if (loc_y_double >= 0 && loc_y_double <= width_s + 1) (loc_y_double - floor(loc_y_double)) > 0.5 ? loc_y = floor(loc_y_double) + 1 : loc_y = floor(loc_y_double);
					else loc_y = -1;

					if (loc_x > height_s - 1 || loc_x < 0 || loc_y > width_s - 1 || loc_y < 0) {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						//char black[3] = {0, 0, 0};
						img_output.at<Vec3b>(loc_i, loc_j) = black;
						// value = [0,0,0]
					}
					else {
						//cv::Vec3d value = img_s.at<Vec3d>(loc_x, loc_y);
						img_output.at<Vec3b>(loc_i, loc_j) = img_s.at<Vec3b>(loc_x, loc_y);
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
						cv::Vec3b value = delta_x_2*delta_y_2*value_1 + delta_x_1*delta_y_2*value_2 + delta_x_2*delta_y_1*value_3 + delta_x_1*delta_y_1*value_4;
						if (value[0] < 0.0 || value[0] > 255.0) value[0] = 0;
						else value[0] = int(value[0]);
						if (value[1] < 0.0 || value[1] > 255.0) value[1] = 0;
						else value[1] = int(value[1]);
						if (value[2] < 0.0 || value[2] > 255.0) value[2] = 0;
						else value[2] = int(value[2]);
						img_output.at<Vec3b>(loc_i, loc_j) = value;
						// value = delta_x_2*delta_y_2*[loc_x_1][loc_y_1] + delta_x_1*delta_y_2*[loc_x_2][loc_y_1] + 
						// delta_x_2*delta_y_1*[loc_x_1][loc_y_2] + delta_x_1*delta_y_1*[loc_x_2][loc_y_2];
					}
					else {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						img_output.at<Vec3b>(loc_i, loc_j) = black;
						// value = [0,0,0]
					}
					break;
				}
				case bicubic: {
					// Do not interpolate at the corner of the image for lacking of confident points
					if (loc_x_double >= 1 && loc_x_double < height_s - 2 && loc_y_double >= 1 && loc_y_double < width_s - 2) {
						cv::Vec3b value;
						int i = floor(loc_x_double);
						int j = floor(loc_y_double);
						double u = loc_x_double - i;
						double v = loc_y_double - j;
						double SU[4] = { bicubic_weight(u + 1),bicubic_weight(u),bicubic_weight(u - 1),bicubic_weight(u - 2) };
						cv::Mat A(1, 4, CV_64FC1, SU);
						double SV[4] = { bicubic_weight(v + 1),bicubic_weight(v),bicubic_weight(v - 1),bicubic_weight(v - 2) };
						cv::Mat C(4, 1, CV_64FC1, SV);

						// For mat multiply convience, deal with single color one time
						Mat B_single = Mat::zeros(4, 4, CV_64FC1);
						for (int channel = 0; channel < 3; channel++) {
							for (int m = -1; m < 3; m++) {
								for (int n = -1; n < 3; n++) {
									B_single.at<double>(m + 1, n + 1) = img_s.at<Vec3b>(i + m, j + n)[channel]; //  < 255 ? img_s.at<Vec3b>(i + m, j + n)[channel] : 255;
								}
							}
							Mat single_color = (A*B_single*C);
							double color = single_color.at<double>(0);
							if (color < 0.0 || color > 255.0) color = 0;
							else (color - floor(color)) > 0.5 ? color = floor(color) + 1 : color = floor(color);
							value[channel] = color;
						}
						img_output.at<Vec3b>(loc_i, loc_j) = value;

					}
					else {
						cv::Vec3b black = cv::Vec3b(0, 0, 0);
						img_output.at<Vec3b>(loc_i, loc_j) = black;
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
						img_output.at<Vec3b>(loc_i, loc_j) = black;
						// value = [0,0,0]

					}
					else {
						cv::Vec3b value = img_s.at<Vec3b>(loc_x, loc_y);
						img_output.at<Vec3b>(loc_i, loc_j) = value;
						// value = [loc_x][loc_y]

					}
					break;
				}
				}
			}
		}

		
		// We can show the altered_image here;
		cv::namedWindow("The new image");
		cv::imshow("The new image", img_output);
		cv::waitKey(3000);

		// Save the altered image;
		img_io_s.image_save(img_output, img_path_t);
		break;
	}

	case BSNaive: {
		// B spline here Naive version

		// normalize the feature points
		// Find the boundary of the face area
		// struct feature to matrix first
		Matrix Mat_Fea_t_h;
		Mat_Fea_t_h.stdFea2stdMat_h(img_fea_t);
		std::vector<vector<double>> Mat_Fea_t = Mat_Fea_t_h.Mat_return(); //3*68 take [1][2]row
		Matrix Mat_Fea_s_v;
		Mat_Fea_s_v.stdFea2stdMat_v(img_fea_s, 0);
		std::vector<vector<double>> Mat_Fea_s = Mat_Fea_s_v.Mat_return(); // 68*2

		std::vector<double> T_locx = Mat_Fea_t[1]; //1*68
		std::vector<double> T_locy = Mat_Fea_t[2];
		double max_t_locx = *std::max_element(std::begin(T_locx), std::end(T_locx));
		double max_t_locy = *std::max_element(std::begin(T_locy), std::end(T_locy));
		double min_t_locx = *std::min_element(std::begin(T_locx), std::end(T_locx));
		double min_t_locy = *std::min_element(std::begin(T_locy), std::end(T_locy));

		std::vector<double> S_locx;
		std::vector<double> S_locy;
		for (int i = 0; i < Mat_Fea_s.size(); i++) {
			S_locx.push_back(Mat_Fea_s[i][0]);
			S_locy.push_back(Mat_Fea_s[i][1]);
		}
		double max_s_locx = *std::max_element(std::begin(S_locx), std::end(S_locx));
		double max_s_locy = *std::max_element(std::begin(S_locy), std::end(S_locy));
		double min_s_locx = *std::min_element(std::begin(S_locx), std::end(S_locx));
		double min_s_locy = *std::min_element(std::begin(S_locy), std::end(S_locy));

		// Normalize the target area to source area.
		// Only modify the relative location
		// I think it is better
		double scale_x = (max_s_locx - min_s_locx) / (max_t_locx - min_t_locx);
		double tran_x = (max_t_locx*min_s_locx - max_s_locx*min_t_locx) / (max_t_locx - min_t_locx);
		double scale_y = (max_s_locy - min_s_locy) / (max_t_locy - min_t_locy);
		double tran_y = (max_t_locy*min_s_locy - max_s_locy*min_t_locy) / (max_t_locy - min_t_locy);
		for (int i = 0; i < T_locx.size(); i++) {
			double temp_x = T_locx[i];
			double temp_y = T_locy[i];
			temp_x = temp_x * scale_x + tran_x;
			temp_y = temp_y * scale_y + tran_y;
			T_locx[i] = temp_x;
			T_locy[i] = temp_y;
		}

		/*
		lets show the modified feature point on the image to compare please.
		*/
		// Notice the difference in the coordinates!!!
		cv::Mat img_feature = img_s.clone();
		IplImage img_feature_arr = img_feature;
		for (int i = 0; i < S_locx.size(); i++) {
			CvPoint pointplot_T;
			pointplot_T.y = int(T_locx[i]); // x and y is just different here...
			pointplot_T.x = int(T_locy[i]);
			cvCircle(&img_feature_arr, pointplot_T, 3, cv::Scalar(0, 255, 0));//, 1, 8, 3 );
			CvPoint pointplot_S;
			pointplot_S.y = int(S_locx[i]);
			pointplot_S.x = int(S_locy[i]);
			cvCircle(&img_feature_arr, pointplot_S, 3, cv::Scalar(255, 0, 0));//, 1, 8, 3 );
		}
		cv::Mat feature_sav = cvarrToMat(&img_feature_arr);
		cv::namedWindow("The feature image");
		cv::imshow("The feature image", feature_sav);
		cv::imwrite("D:\\ATsinghua\\value_analysis\\features.jpg", feature_sav);
		cvWaitKey(1000);
		
		Matrix Mat_Fea_t_scale;
		std::vector<vector<double>> Fea_t_scale; // 2*68
		Fea_t_scale.push_back(T_locx);
		Fea_t_scale.push_back(T_locy);
		Mat_Fea_t_scale.Mat_set(Fea_t_scale);

		// Create two grid translation graph (X Y direction divided)
		int height_s = img_s.size().height;
		int width_s = img_s.size().width;
		int height_t = img_t.size().height;
		int width_t = img_t.size().width;
		std::vector<double> D_locx;
		std::vector<double> D_locy;
		double max_derative = 0;
		for (int i = 0; i < S_locx.size(); i++) {
			D_locx.push_back(T_locx[i] - S_locx[i]);
			D_locy.push_back(T_locy[i] - S_locy[i]);
			double derative = sqrt(pow((T_locx[i] - S_locx[i]), 2) + pow((T_locy[i] - S_locy[i]), 2));
			if (derative > max_derative) {
				max_derative = derative;
			}
		}

		int grid_size = 30; // floor(max_derative);
		cout << "Grid size: " << grid_size << endl;
		cv::Mat Graph_tran_grid_x = cv::Mat::zeros(floor(height_s / grid_size), floor(width_s / grid_size), CV_64FC1);
		cv::Mat Graph_tran_grid_y = cv::Mat::zeros(floor(height_s / grid_size), floor(width_s / grid_size), CV_64FC1);
		int grid_height = Graph_tran_grid_x.size().height;
		int grid_width = Graph_tran_grid_x.size().width;

		// Use the interative way to solve the grid and plot the calculated feature points to prove the right of grid

		grid_init(D_locx, S_locx, D_locy, S_locy, Graph_tran_grid_x, Graph_tran_grid_y, grid_size, max_s_locx, min_s_locx, max_s_locy, min_s_locy);
		double de_sum_x = 1000.0;
		double de_sum_y = 1000.0;
		int max_inter = 1000;
		int inter_time = 0;
		while (de_sum_x > 10.0 || de_sum_y > 10.0) {
			inter_time += 1;
			if (inter_time > max_inter) {
				break;
			}
			// warp the New S_loc first;
			std::vector<double> T_locx_new;
			std::vector<double> T_locy_new;
			de_sum_x = 0.0;
			de_sum_y = 0.0;
			for (int p = 0; p < S_locx.size(); p++) {
				// i,j is the cooresponding pixel in the img
				int i = S_locx[p];
				int j = S_locy[p];
				double u = double(i) / grid_size - floor(double(i) / grid_size);
				double v = double(j) / grid_size - floor(double(j) / grid_size);
				double grid_i = floor(double(i) / grid_size);
				double grid_j = floor(double(j) / grid_size);
				double x_list = 0.0;
				double y_list = 0.0;
				for (int near_i = -1; near_i < 3; near_i++) {
					for (int near_j = -1; near_j < 3; near_j++) {
						if (grid_i + near_i < 0 || grid_i + near_i > grid_height - 1 || grid_j + near_j < 0 || grid_j + near_j > grid_width - 1) {
							x_list += 0.0;
							y_list += 0.0;
						}
						else {
							x_list += B_spline_weight(u, near_i + 1)*B_spline_weight(v, near_j + 1)*Graph_tran_grid_x.at<double>(grid_i + near_i, grid_j + near_j);
							y_list += B_spline_weight(u, near_i + 1)*B_spline_weight(v, near_j + 1)*Graph_tran_grid_y.at<double>(grid_i + near_i, grid_j + near_j);
						}
					}
				}
				T_locx_new.push_back(S_locx[p] + x_list);
				T_locy_new.push_back(S_locy[p] + y_list);
				// To check if the recalculate the list of feature point is close to what we want
				de_sum_x += pow((T_locx[p] - T_locx_new[p]), 2.0);
				de_sum_y += pow((T_locy[p] - T_locy_new[p]), 2.0);
			}

			std::cout << ">> For the " << inter_time << " the derative is: " << de_sum_x << " " << de_sum_y << endl;
			std::vector<double> D_locx_tem;
			std::vector<double> D_locy_tem;
			for (int i = 0; i < S_locx.size(); i++) {
				D_locx_tem.push_back((T_locx[i] - T_locx_new[i]));
				D_locy_tem.push_back((T_locy[i] - T_locy_new[i]));
			}

			grid_init(D_locx_tem, S_locx, D_locy_tem, S_locy, Graph_tran_grid_x, Graph_tran_grid_y, grid_size, max_s_locx, min_s_locx, max_s_locy, min_s_locy);
			// I need to plot the possible target feature points.
			if (de_sum_x < 10.0 && de_sum_y < 10.0) {
				cv::Mat img_feature_calc = img_s.clone();
				IplImage img_feature_arr = img_feature_calc;
				for (int i = 0; i < S_locx.size(); i++) {
					CvPoint pointplot_T;
					pointplot_T.y = int(T_locx[i]); // x and y is just different here...
					pointplot_T.x = int(T_locy[i]);
					cvCircle(&img_feature_arr, pointplot_T, 3, cv::Scalar(0, 255, 0));
					CvPoint pointplot_T_calc;
					pointplot_T_calc.y = int(T_locx_new[i]);
					pointplot_T_calc.x = int(T_locy_new[i]);
					cvCircle(&img_feature_arr, pointplot_T_calc, 2, cv::Scalar(0, 0, 255));//, 1, 8, 3 );
				}
				cv::Mat feature_calc_sav = cvarrToMat(&img_feature_arr);
				cv::namedWindow("The feature calculate image");
				cv::imshow("The feature calculate image", feature_calc_sav);
				cv::imwrite("D:\\ATsinghua\\value_analysis\\features_calc.jpg", feature_calc_sav);
				cvWaitKey(3000);
			}
		}

		cv::Mat Graph_tran_x = cv::Mat::zeros(height_s, width_s, CV_64FC1);
		cv::Mat Graph_tran_y = cv::Mat::zeros(height_s, width_s, CV_64FC1);
		for (int i = 0; i < height_s; i++) {
			for (int j = 0; j < width_s; j++) {
				// grid: x-1 ~ x+2 y-1 ~ y+2
				// if grid = 1, then the u = x/grid - floor(x/grid) = 0 then the weight is const number
				double u = double(i) / grid_size - floor(double(i) / grid_size);
				double v = double(j) / grid_size - floor(double(j) / grid_size);
				double grid_i = floor(double(i) / grid_size);
				double grid_j = floor(double(j) / grid_size);
				double x_list = 0.0;
				double y_list = 0.0;
				for (int near_i = -1; near_i < 3; near_i++) {
					for (int near_j = -1; near_j < 3; near_j++) {
						if (grid_i + near_i < 0 || grid_i + near_i > grid_height - 1 || grid_j + near_j < 0 || grid_j + near_j > grid_width - 1) {
							x_list += 0.0;
							y_list += 0.0;
						}
						else {
							x_list += B_spline_weight(u, near_i + 1)*B_spline_weight(v, near_j + 1)*Graph_tran_grid_x.at<double>(grid_i + near_i, grid_j + near_j);
							y_list += B_spline_weight(u, near_i + 1)*B_spline_weight(v, near_j + 1)*Graph_tran_grid_y.at<double>(grid_i + near_i, grid_j + near_j);
						}
					}
				}
				//x_list = fabs(x_list) > double(2.0*grid_size) ? x_list / fabs(x_list) * grid_size : x_list;
				//y_list = fabs(y_list) > double(2.0*grid_size) ? y_list / fabs(y_list) * grid_size : y_list;
				Graph_tran_x.at<double>(i, j) = x_list;
				Graph_tran_y.at<double>(i, j) = y_list;
			}
		}

		cv::Mat img_output_ini = cv::Mat::zeros(height_s, width_s, CV_8UC3);
		for (int loc_i = 0; loc_i < height_s; loc_i++) {
			for (int loc_j = 0; loc_j < width_s; loc_j++) {
				double x_list = Graph_tran_x.at<double>(loc_i, loc_j);
				double y_list = Graph_tran_y.at<double>(loc_i, loc_j);
				double loc_x_double = loc_i + x_list;
				double loc_y_double = loc_j + y_list;

				double loc_x_new = (loc_x_double - floor(loc_x_double)) > 0.5 ? floor(loc_x_double) + 1 : floor(loc_x_double);
				double loc_y_new = (loc_y_double - floor(loc_y_double)) > 0.5 ? floor(loc_y_double) + 1 : floor(loc_y_double);
				if (loc_x_new > -1 && loc_x_new < height_s && loc_y_new > -1 && loc_y_new < width_s) {
					img_output_ini.at<Vec3b>(loc_x_new, loc_y_new) = img_s.at<Vec3b>(loc_i, loc_j);
				}
			}
		}

		// Alter the image

		cv::Mat img_output = cv::Mat::zeros(height_s, width_s, CV_8UC3); //cv::Mat::zeros(height_s, width_s, CV_8UC3);
		for (int loc_i = 0; loc_i < height_s; loc_i++) {
			for (int loc_j = 0; loc_j < width_s; loc_j++) {
				double loc_x_double = 0.0;
				double loc_y_double = 0.0;
				if (img_output_ini.at<Vec3b>(loc_i, loc_j) == cv::Vec3b(0, 0, 0))
				{
					// Interpolate here with three ways
					switch (inter_input) {
					case nearest: {
						int loc_x = loc_i;
						int loc_y = loc_j;
						int near_n = 10;
						for (int near_i = -2; near_i < 3; near_i++) {
							for (int near_j = -2; near_j < 3; near_j++) {
								int temp_x = loc_x + near_i;
								int temp_y = loc_y + near_j;
								temp_x = temp_x > -1 ? temp_x : 0;
								temp_y = temp_y > -1 ? temp_y : 0;
								temp_x = temp_x < height_s ? temp_x : height_s - 1;
								temp_y = temp_y < width_s ? temp_y : width_s - 1;
								if (img_output_ini.at<Vec3b>(temp_x, temp_y) != cv::Vec3b(0, 0, 0) && (abs(near_j) + abs(near_i)) < near_n) {
									loc_x = temp_x;
									loc_y = temp_y;
									near_n = abs(near_j) + abs(near_i);
								}
							}
						}

						if (loc_x > height_s - 1 || loc_x < 0 || loc_y > width_s - 1 || loc_y < 0) {
							cv::Vec3b black = cv::Vec3b(0, 0, 0);
							//char black[3] = {0, 0, 0};
							img_output.at<Vec3b>(loc_i, loc_j) = black;
							// value = [0,0,0]
						}
						else {
							//cv::Vec3d value = img_s.at<Vec3d>(loc_x, loc_y);
							img_output.at<Vec3b>(loc_i, loc_j) = img_output_ini.at<Vec3b>(loc_x, loc_y);
						}
						break;
					}
					}
				}
				else {
					img_output.at<Vec3b>(loc_i, loc_j) = img_output_ini.at<Vec3b>(loc_i, loc_j);
				}
			}
		}

		// We can show the altered_image here;
		cv::namedWindow("The new image");
		cv::imshow("The new image", img_output);
		cv::waitKey(3000);

		// Save the altered image;
		img_io_s.image_save(img_output, img_path_t);
		break;
	}

	default:
		break;
	}
	
	return 0;
}
