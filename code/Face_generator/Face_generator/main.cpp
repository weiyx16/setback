#include "image_inout.h"

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

	// Save the altered image;
	// img_io_s.image_save(img_ori);

}
