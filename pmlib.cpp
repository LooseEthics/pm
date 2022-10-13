
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
//#include <opencv2/opencv.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace cv;
using namespace std;


// Geometry conversion functions

// Frame parser functions

void fp_brute_force(Mat* img, int* max_row_arr){
	// Process all columns sequentially
}

// Column search functions

int cs_full_col_max(Mat* img, int col) {
	// Brute force search entire column for brightest pixel
	int maxrow = 0;
	int maxbr = 0;
	for (int row = 0; row < (* img).rows; row++) {
		Vec3b& px = (* img).at<Vec3b>(row, col);
		if (px[0] + px[1] + px[2] > maxbr) {
			maxbr = px[0] + px[1] + px[2];
			maxrow = row;
		}
	}
	return maxrow;
}

int cs_neigh_thres_max(Mat* img, int col, int last_max_row, int thres, int range) {
	// Thresholded search within {range} pixels of {last_max_row}, if not found do brute force
	return 0;
}

int cs_thres_center(Mat* img, int col, int last_max_row, int thres) {
	// Find range around {last_max_row} above {thres}, return center pixel of range
	return 0;
}

int cs_thres_median(Mat* img, int col, int last_max_row, int thres) {
	// Find range around {last_max_row} above {thres}, return median pixel of range
	return 0;
}

// Curvature correction functions