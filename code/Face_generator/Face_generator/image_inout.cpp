#include "image_inout.h"

using namespace std;
using namespace cv;

// Load the image
cv::Mat image_inout::image_load()
{
	std::cout << ">> Loading image from " << img_path << endl;
	Mat img_in = cv::imread(img_path, IMREAD_COLOR);
	std::cout << ">> Finished loading" << endl;
	return img_in;
}

// Load the correspond feature.txt
std::vector<Location_fea> image_inout::image_features_load()
{
	std::vector<Location_fea> fea_data;
	std::string file_path;
	file_path = img_name + ".txt";

	int num = 0;
	std::cout << ">> Loading image features from " << file_path << endl;

	std::ifstream fea_load;
	fea_load.open(file_path);
	if (fea_load.is_open())
	{
		while (!fea_load.eof())
		{
			Location_fea fea_point;
			fea_load >> fea_point.loc_y;
			fea_load >> fea_point.loc_x;
			fea_data.push_back(fea_point);
			num++;
			// cout << "Input:: " << num << fea_point.loc_x << ' ' << fea_point.loc_y << endl;
		}
		fea_load.close();
		fea_data.pop_back(); // Seems to have a problem of reloading the final line.
		std::cout << ">> Successfully load " << num-1 << " points" << endl;
	}
	else std::cout << ">> Fail to load image features";
	return fea_data;
}

// Save the img_altered with a new file name
void image_inout::image_save(cv::Mat img_save, std::string img_path_t)
{
	std::size_t loc_l = img_path_t.find_last_of('\\');
	if (loc_l == -1){
		loc_l = img_path_t.find_last_of('/');
	}
	std::size_t loc_h = img_path_t.find_last_of('.');
	std::string img_t_name = img_path_t.substr(loc_l + 1, loc_h - loc_l - 1);

	cv::String img_save_path = img_name + "_alteredby_" + img_t_name +".jpg";
	std::cout << ">> Saving image to " << img_save_path << endl;
	cv::imwrite(img_save_path, img_save);
	std::cout << ">> Finished saving" << endl;
}

