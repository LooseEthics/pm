// opencv test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>

using namespace cv;
using namespace std;

int main_img()
{
    string path = "Resources/test.png";
    Mat img = imread(path);
    imshow("Image", img);
    waitKey(0);
    return 0;
}

int main_vid() {
    string path = "Resources/vid.mp4";
    VideoCapture cap(path);
    Mat img;

    while (true){
        cap.read(img);
        imshow("Image", img);
        waitKey(10);
    }
    return 0;
}

int main_cam() {
    VideoCapture cap(0);
    Mat img;

    while (true) {
        cap.read(img);
        imshow("Image", img);
        waitKey(10);
    }
    return 0;
}

int main() {
    string path = "Resources/PM_test_1.png";
    Mat img = imread(path);


    int optimize = 0;
    int tag = 1;
    // Brightest array setup
    int *maxrow_arr = (int *)malloc(sizeof(int) * img.cols);
    for (int i = 0; i < img.cols; i++) {
        maxrow_arr[i] = -1;
    }

    if (optimize) {
        //float thrcols[1] = { 0.5 };
        float thrcols[3] = { 0.25, 0.5, 0.75 };
        float thrscale = 0.7;
        int thrsum = 0;
        int thrcnt = 0;

        // Bruteforce threshold columns
        // Get brightness threshold
        for (float ratio : thrcols) {
            int maxrow = 0;
            int maxbr = 0;
            int i = (int)(ratio * img.cols);
            for (int j = 0; j < img.rows; j++) {
                Vec3b& px = img.at<Vec3b>(j, i);
                if (px[0] + px[1] + px[2] > maxbr) {
                    maxbr = px[0] + px[1] + px[2];
                    maxrow = j;
                }
            }
            maxrow_arr[i] = maxrow;
            thrsum += maxbr;
            thrcnt += 1;
            
            if (tag) {
                Vec3b& px = img.at<Vec3b>(maxrow, i);
                px.val[0] = 0;
                px.val[1] = 0;
                px.val[2] = 255;
            }
        }

        // Optimization constants
        int thr = (int)(thrscale * thrsum / thrcnt);
        int voffset = 8;
        int hoffset = 8;


        // Left pass
        for (int i = (int)(0.5*img.cols - 1); i >= 0; i--) {
            int maxrow = 0;
            int maxbr = 0;
            //// TODO
        }

        // Right pass
        //// TODO

    }
    else {
        for (int i = 0; i < img.cols; i++) {
            int maxrow = 0;
            int maxbr = 0;
            for (int j = 0; j < img.rows; j++) {
                Vec3b& px = img.at<Vec3b>(j, i);
                if (px[0] + px[1] + px[2] > maxbr) {
                    maxbr = px[0] + px[1] + px[2];
                    maxrow = j;
                }
            }
            maxrow_arr[i] = maxrow;
            
            if (tag) {
                Vec3b& px = img.at<Vec3b>(maxrow, i);
                px.val[0] = 0;
                px.val[1] = 0;
                px.val[2] = 255;
            }
        }
    }

    for (int i = 0; i < img.cols; i++) {
        cout << maxrow_arr[i] << ", ";
    }
    imshow("Image", img);
   
    waitKey(0);
    return 0;
}
