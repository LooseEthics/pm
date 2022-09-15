// Profilometer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//Program flow:
// Initialize PMeter
// Set geometry constants - sigma, dist, radius
// Set camera constants - pixel dimensions, fov angles
// Precalculate other constants - trig expressions and such
// ---Precalculate curvature correction---
//   Nevermind, this varies with horizontal position and therefore both x and y, would need full screen 2d array
//   Maybe there's an approximation that's fixed in one of those and varies in the other?
// Per frame:
//   Preprocessing - Calculate threshold for identifying intensity spikes
//   Build incidence pixel (IP) array
//     Laser likely centered, may not extend across whole fov -> start in middle of x range, not edge
//     Prioritize neighborhood search to cut down number of operations
//     Option 1: maximal intensity pixel is IP
//     Option 2: find groups of pixels above threshold, use weighted average of the segment as IP
//   Calculate real space depth of IPs
//   Apply curvature correction
//   Return depth array
//

#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <opencv2/opencv.hpp>

class PMeter {
public:
    int height, width;
    // camera pixel dimensions
    double alpha, beta;
    // camera angle dimensions
    // full fov is 2* this angle, vertical and horizontal respectively
    // 0 < alpha < sigma
    // alpha >= sigma causes phi singularity
    // 0 < beta < pi/2
    double sigma;
    // angle between laser plane and camera center line
    // 0 < sigma < pi/2
    double dist;
    // distance between camera and laser

    double tan_alpha;
    double tan_beta;
    double* tan_phi_arr;
    double* tan_theta_arr;
    double sin_sigma;
    double cos_sigma;
    double p_h;


    void fastinit(int h, int w, double a, double b, double s, double d) {
        // class variable initializer
        height = h;
        width = w;
        alpha = a;
        beta = b;
        sigma = s;
        dist = d;
        trig_precalc();
    }

    void trig_precalc() {
        tan_alpha = tan(alpha);
        tan_beta = tan(beta);
        tan_phi_arr = (double*)malloc(sizeof(double) * height);
        for (int i = 0; i < height; i++) {
            //tan_phi_arr[i] = tan_alpha / height * (2 * i + 1 - height);
            tan_phi_arr[i] = tan_alpha / height * (2 * i - height);
        }
        tan_theta_arr = (double*)malloc(sizeof(double) * width);
        for (int i = 0; i < width; i++) {
            //tan_theta_arr[i] = -tan_beta / width * (2 * i - 1 - width);
            tan_theta_arr[i] = -tan_beta / width * (2 * i - width);
        }
        sin_sigma = sin(sigma);
        cos_sigma = cos(sigma);
        p_h = 2 * dist / height * tan_alpha / sin_sigma;

    }

    void print_conf() {
        std::cout << "Height: " << height << "\n";
        std::cout << "Width: " << width << "\n";
        std::cout << "Alpha: " << alpha << "\n";
        std::cout << "Beta: " << beta << "\n";
        std::cout << "Sigma: " << sigma << "\n";
        std::cout << "Dist: " << dist << "\n";
        //std::cout << "Phi: " << phi << "\n";
        //std::cout << "Theta: " << theta << "\n";
    }
    
    double get_phi(int y) {
        // phi is angle below camera center line of given y coordinate - above center is negative
        // returns phi at pixel center, removing the +1 shifts this to low edge
        // alpha >= phi >= -alpha
        double phi = atan(tan_phi_arr[y]);
        return phi;
    }

    double get_theta(int x) {
        // theta is horizontal angle between center line and given x coordinate
        // returns theta at pixel center, removing the -1 shifts this to low edge
        // beta >= theta >= -beta
        // if viewed from above, positive theta means pixel vector is to the left of center, i.e. x < width/2
        double theta = atan(tan_theta_arr[x]);
        return theta;
    }

    double sin_sp(int y) {
        // sin( sigma + phi)
        return sin_sigma + cos_sigma * tan_phi_arr[y] / sqrt(1 + tan_phi_arr[y]);
    }

    double cos_sp(int y) {
        // cos( sigma + phi)
        return cos_sigma + sin_sigma * tan_phi_arr[y] / sqrt(1 + tan_phi_arr[y]);
    }

    double pixel_delta_approx(int y) {
        // returns approximation of change in elevation from y to y+1
        // assumes that epsilon=0
        //double delta = 2 * dist * tan(alpha) * cos(phi) / (height * sin(sigma + phi) * sin(sigma + phi));
        double delta = p_h * sin_sigma / (sin_sigma + cos_sigma * tan_phi_arr[y]);
        return delta;
    }

    double pixel_delta_exact(int y) {
        // returns exact change in elevation from y to y+1
        //double tan_eps = 2 * tan(alpha) * cos(phi) * cos(phi) / (height + 2 * tan(alpha) * sin(phi) * cos(phi));
        double tan_phi = tan_phi_arr[y];
        double tan_eps = 2 * tan_alpha / (height * (1 + tan_phi * tan_phi) + 2 * tan_alpha * tan_phi);
        std::cout << "Tan(eps) at phi=" << get_phi(y) << ": " << tan_eps << "\n";
        //double delta = 2 * dist * tan(alpha) * cos(phi) * (cos(phi) - sin(phi) * tan_eps) / (height * sin(sigma + phi) * (sin(sigma + phi) + cos(sigma + phi) * tan_eps));
        //double delta = dist * tan_eps / (sin_sp(y) * (sin_sp(y) + tan_eps * cos_sp(y)));
        double delta = 2 * dist * tan_alpha * (1 - tan_phi_arr[y] * tan_eps) / (height * sin_sp(y) * (sin_sp(y) + cos_sp(y) * tan_eps));
        return delta;
    }
    double get_elevation(int y) {
        // elevation above horizontal camera plane
        // singularity if phi = -sigma
        // flat projection surface -> x coordinate irrelevant
        //double e = dist / tan(phi + sigma);
        double e = dist * cos_sp(y) / sin_sp(y);
        return e;
    }
    
    double get_horizontal_dist(int x, int y) {
        // get horizontal distance from center line, same sign convention as theta
        // generally depends on both x and y
        //double f = dist * tan(theta) / sin(sigma + phi);
        double f = dist * tan_theta_arr[x] / sin_sp(y);
        return f;
    }

    /*
    void approx_test_phi(double new_phi) {
        // print exact and approximate delta and the error at given phi
        // changes phi
        phi = new_phi;
        double da = pixel_delta_approx();
        std::cout << "Delta approx.: " << da << "\n";
        double de = pixel_delta_exact();
        std::cout << "Delta exact: " << de << "\n";
        std::cout << "Error absolute: " << (de - da) << "\n";
        std::cout << "Error relative: " << (de - da)/de << "\n";
    }
    */

    void approx_test(int y) {
        double da = pixel_delta_approx(y);
        std::cout << "Delta approx.: " << da << "\n";
        double de = pixel_delta_exact(y);
        std::cout << "Delta exact: " << de << "\n";
        std::cout << "Error absolute: " << (de - da) << "\n";
        std::cout << "Error relative: " << (de - da) / de << "\n";
    }

    /*
    void edge_test() {
        // run approx_test_phi at bottom edge, phi=0, and top edge
        // preserves phi
        double old_phi = phi;
        std::cout << "\nMax phi\n";
        approx_test_phi(alpha);
        std::cout << "\nZero phi\n";
        approx_test_phi(0.0);
        std::cout << "\nMin phi\n";
        approx_test_phi(-1*alpha);
        phi = old_phi;
    }
    */

    void edge_test() {
        // run approx_test_phi at y = {0, height/2, height-1}
        std::cout << "\ny = h - 1\n";
        approx_test(height - 1);
        std::cout << "\ny = h/2\n";
        approx_test(height/2);
        std::cout << "\ny = 0\n";
        approx_test(0);
    }
};


int main()
{
    PMeter PM;
    // fastinit(int h, int w, double a, double b, double s, double d)
    PM.fastinit(3120, 4208, M_PI / 180 * 21.5, M_PI / 180 * 43, M_PI / 180 * 51, 0.05);
    PM.print_conf();


    std::cout << "phi at y=50 is " << PM.get_phi(50) << "\n";
    PM.edge_test();


    ////cv::Mat image = cv::imread("./PM_test_0.bmp");
    //cv::Mat image = cv::imread("./PM_test_1.png");
    ////cv::Mat image = cv::imread("./PM_test_2.png");
    //
    //cv::String windowName = "imtest";
    //cv::namedWindow(windowName);
    //cv::imshow(windowName, image);
    //cv::waitKey(0);


    std::cout << "\n===\nEnd of main()\n===\n";
}

