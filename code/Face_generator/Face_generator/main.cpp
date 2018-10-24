#include "image_inout.h"
#include <iostream>

using namespace std;
using namespace cv;

void main() {

	// Load image with its features;
	cv::String img_path = "D:\\ATsinghua\\value_analysis\\picture\\1.jpg";
	image_inout img_io(img_path);

	cv::Mat img_ori = img_io.image_load();
	std::vector<Location_fea> img_fea = img_io.image_features_load();
	cv::namedWindow("The origin image");
	cv::imshow("The origin image", img_ori);
	waitKey(3000);





	// Save the altered image;
	img_io.image_save(img_ori);
	
	

}
