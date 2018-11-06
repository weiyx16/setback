#include "image_inout.h"
#include "Matrix.h"

using namespace std;
using namespace cv;

void main() {

	// Load source image with its features;
	cv::String img_path_s = "D:\\value_analysis\\picture\\1.jpg";
	// New a object to load the image;
	image_inout img_io_s(img_path_s);

	// Load the img to Mat type
	cv::Mat img_s = img_io_s.image_load();
	// Load the correspond image
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
	// Load the correspond image
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
	vector<vector<double>> one_matrix;
	for (int i = 0; i < 3; i++) {
		vector<double> temp;
		for (int j = 0; j < 2; j++) {
			temp.push_back(0.0);
		}
		one_matrix.push_back(temp);
	}
	Mat_Fea_s_v.Mat_row_concat(one_matrix);


	// Calculate the inverse of L to get paras
	// Mat_L_t.Mat_inv();
	// vector<vector<double>> paras;
	// paras = Mat_L_t.Mat_mul(Mat_Fea_s_v.Mat_return());

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



	Matrix Mat_paras;
	Mat_paras.Mat_set(paras);


	// Save the altered image;
	// img_io_s.image_save(img_ori);

}
