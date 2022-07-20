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
    double phi, theta;
    // coordinate angles

    void fastinit(int h, int w, double a, double b, double s, double d, double p, double t) {
        //TODO: overload this
        // class variable initializer
        height = h;
        width = w;
        alpha = a;
        beta = b;
        sigma = s;
        dist = d;
        phi = p;
        theta = t;
    }

    void print_conf() {
        std::cout << "Height: " << height << "\n";
        std::cout << "Width: " << width << "\n";
        std::cout << "Alpha: " << alpha << "\n";
        std::cout << "Beta: " << beta << "\n";
        std::cout << "Sigma: " << sigma << "\n";
        std::cout << "Dist: " << dist << "\n";
        std::cout << "Phi: " << phi << "\n";
        std::cout << "Theta: " << theta << "\n";
    }

    double get_phi(int y) {
        // phi is angle below camera center line of given y coordinate - above center is negative
        // returns phi at pixel center, removing the +1 shifts this to low edge
        // alpha >= phi >= -alpha
        phi = atan(tan(alpha)/height * (2 * y + 1 - height));
        return phi;
    }

    double get_theta(int x) {
        // theta is horizontal angle between center line and given x coordinate
        // returns theta at pixel center, removing the -1 shifts this to low edge
        // beta >= theta >= -beta
        // if viewed from above, positive theta means pixel vector is to the left of center, i.e. x < width/2
        theta = atan(tan(beta) / width * (width - 2 * x - 1));
        return theta;
    }

    double pixel_delta_approx() {
        // returns approximation of change in elevation from y to y+1
        // assumes that epsilon=0
        double delta = 2 * dist * tan(alpha) * cos(phi) / (height * sin(sigma + phi) * sin(sigma + phi));
        return delta;
    }

    double pixel_delta_exact() {
        // returns exact change in elevation from y to y+1
        double tan_eps = 2 * tan(alpha) * cos(phi) * cos(phi) / (height + 2 * tan(alpha) * sin(phi) * cos(phi));
        std::cout << "Tan(eps) at phi=" << phi << ": " << tan_eps << "\n";
        double delta = 2 * dist * tan(alpha) * cos(phi) * (cos(phi) - sin(phi) * tan_eps) / (height * sin(sigma + phi) * (sin(sigma + phi) + cos(sigma + phi) * tan_eps));
        return delta;
    }

    double get_elevation() {
        // elevation above horizontal camera plane
        // singularity if phi = -sigma
        // flat projection surface -> x coordinate irrelevant
        double e = dist / tan(phi + sigma);
        return e;
    }

    double get_horizontal_dist() {
        // get horizontal distance from center line, same sign convention as theta
        // generally depends on both x and y
        double f = dist * tan(theta) / sin(sigma + phi);
        return f;
    }

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

};


int main()
{
    PMeter PM;
    PM.fastinit(1000, 1000, M_PI / 6, M_PI / 6, M_PI / 3, 5.0, 0.0, 0.0);
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

