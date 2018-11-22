#pragma once
#ifndef IMAGE_INOUT
#define IMAGE_INOUT

#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace cv;

// In the feature file given: y(column) first then the x(row)
struct Location_fea
{
	float loc_x;
	float loc_y;
};

class image_inout
{
public:
	// load img
	cv::Mat image_load();
	// load feature txt
	std::vector<Location_fea> image_features_load();
	// save img
	void image_save(cv::Mat img_save, std::string img_path_t);
	image_inout(cv::String path) : img_path(path) {
		type_loc = img_path.find_first_of('.');
		img_name = img_path.substr(0, type_loc);
	};
private:
	cv::String img_path;
	std::size_t type_loc = 0;
	std::string img_name = "no_img";
};

#endif // !IMAGE_INOUT.h
