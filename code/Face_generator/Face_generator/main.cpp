#include "image_inout.h"
#include <iostream>
#include <Eigen/Dense> 

using namespace std;
using namespace cv;
using namespace Eigen;

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

	/* ---- how to use Eigen
	Eigen::Matrix2d a;
	a << 1, 2,
		3, 4;

	Eigen::MatrixXd b(2, 2);
	b << 2, 3,
		1, 4;

	std::cout << "a + b =\n" << a + b << std::endl;
	std::cout << "a - b =\n" << a - b << std::endl;
	std::cout << "Doing a += b;" << std::endl;
	a += b;
	std::cout << "Now a =\n" << a << std::endl;

	Eigen::Vector3d v(1, 2, 3);
	Eigen::Vector3d w(1, 0, 0);

	std::cout << "-v + w - v =\n" << -v + w - v << std::endl;
	*/
}
