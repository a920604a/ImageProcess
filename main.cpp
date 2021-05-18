/*
 * @Author: yuan
 * @Date: 2021-05-13 02:27:58
 * @LastEditTime: 2021-05-18 09:27:14
 * @FilePath: /面試題目_20210121/main.cpp
 */
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <iostream>
using namespace cv;
using namespace std;
#include "histogram.cpp"
#include "rotate.cpp"
#include "scale.cpp"
#include "noise.cpp"
#include "inpaint.cpp"

int main()
{

    std::string image_path = "./lena.bmp";
    Mat img = imread(image_path, IMREAD_GRAYSCALE);

    histogram_fun(img); // https://followtutorials.com/2013/01/intensity-histogram-using-c-and-opencv-image-processing.html

    Mat ret_rotate;
    ret_rotate = rotate(img, 45);
    imwrite("rotate.bmp", ret_rotate);
    cout << "-------------rotate finished--------------" << endl;
    Mat ret_scale;
    ret_scale = scale(rotate(img, 45), 800, 600); //https://blog.csdn.net/pentiumCM/article/details/104720100
    imwrite("r_scale.bmp", ret_scale);
    cout << "-------------scale finished--------------" << endl;
    // ret = rotate(scale(img, 800, 600), 45); // not shift

    noise_type r = noise(img);
    imwrite("noise.bmp", r.src);
    imwrite("filter.bmp", r.dst);
    cout << "-------------noise finished--------------" << endl;
    Mat ret_repair;
    ret_repair = mask_inpaint(img);
    imwrite("repair.bmp", ret_repair);
    cout << "-------------inpaint finished--------------" << endl;

    return 0;
}
