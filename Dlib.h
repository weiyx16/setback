#include <dlib/opencv.h>
#include <opencv2/opencv.hpp>
#include <dlib/image_processing/frontal_face_detector.h>
#include <dlib/image_processing/render_face_detections.h>
#include <dlib/image_processing.h>
#include <dlib/gui_widgets.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
//#include <opencv2/highgui.h> 

using namespace dlib;
using namespace std;
using namespace cv;

int main()
{
		//image_window win;
 
		// Load face detection and pose estimation models.
		dlib::frontal_face_detector detector = get_frontal_face_detector();
		dlib::shape_predictor pose_model;
		string modelpath = "./shape_predictor_68_face_landmarks.dat";
		dlib::deserialize(modelpath) >> pose_model;
 
		// load an image using opencv and convert it to bgr_pixel;
		
		string imgpath = "";
		cv::Mat img = cv::imread(imgpath);
		// dlib::cv_image<rgb_pixel> dlib_img(img); // only stores pointer, no deep copy  
		dlib::array2d<bgr_pixel> = dlib_img(img.rows, img.cols);
		// 2-for-loop
		
		
		// Detect faces 
		std::vector<rectangle> faces = detector(dlib_img);
		// Find the pose of each face.
		std::vector<full_object_detection> shapes;
		for (unsigned long face_i = 0; face_i < faces.size(); ++face_i)
		{
			shapes.push_back(pose_model(dlib_img, faces[face_i]));
			for (int fea_i = 0; fea_i < 68; fea_i++) {
				circle(img, cvPoint(shapes[face_i].part(fea_i).x(), shapes[face_i].part(fea_i).y()), 3, cv::Scalar(0, 0, 255), -1);
				//	shapes[0].part(i).x();//68个
				cout << shapes[face_i].part(fea_i).x() << ' ' << shapes[face_i].part(fea_i).y() << endl;
			}
		}

		//Display it all on the screen
		imshow("Dlib特征点", img);
}