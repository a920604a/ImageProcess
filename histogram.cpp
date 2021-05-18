/*
 * @Author: yuan
 * @Date: 2021-05-18 04:59:46
 * @LastEditTime: 2021-05-18 09:25:55
 * @FilePath: /面試題目_20210121/histogram.cpp
 */

#include <opencv2/imgproc/imgproc.hpp> //line Point

void histogram_fun(Mat img)
{
    int histogram[256];

    for (int i = 0; i < 256; i++)
    {
        histogram[i] = 0;
    }
    for (int y = 0; y < img.rows; y++)
    {
        for (int x = 0; x < img.cols; x++)
        {
            histogram[(int)img.at<uchar>(y, x)]++;
        }
    }

    int hist_w = 256;
    int hist_h = 512;

    Mat histImage(hist_h, hist_w, CV_8UC1, Scalar(255, 255, 255));

    int max = histogram[0];
    for (int i = 1; i < 256; i++)
    {
        if (max < histogram[i])
        {
            max = histogram[i];
        }
    }

    for (int i = 0; i < 255; i++)
    {
        histogram[i] = ((double)histogram[i] / max) * histImage.rows;
    }

    for (int i = 0; i < 255; i++)
    {
        line(histImage,
             Point(i, hist_h),                //起點
             Point(i, hist_h - histogram[i]), //終點
             Scalar(0, 255, 0), 1, 8, 0);
    }
    imwrite("hist.bmp", histImage);
}
