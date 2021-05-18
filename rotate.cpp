/*
 * @Author: yuan
 * @Date: 2021-05-18 04:59:21
 * @LastEditTime: 2021-05-18 09:29:51
 * @FilePath: /面試題目_20210121/rotate.cpp
 */

#define PI 3.1415926


vector<int> affine(int x, int y, Mat m)
{
    // 𝚍𝚜𝚝(x,y)=𝚜𝚛𝚌(𝙼11x+𝙼12y+𝙼13,𝙼21x+𝙼22y+𝙼23)

    double new_x = m.at<double>(0, 0) * x + m.at<double>(0, 1) * y + m.at<double>(0, 2);
    double new_y = m.at<double>(1, 0) * x + m.at<double>(1, 1) * y + m.at<double>(1, 2);

    return {(int)round(new_x), (int)round(new_y)};
}




Mat rotate(Mat src, double angle) 
{

    int src_h = src.rows, src_w = src.cols;
    double alpha = cos(angle * PI / 180.0);
    double beta = sin(angle * PI / 180.0);


    Mat r = (cv::Mat1d(2, 3) << alpha, beta, (1 - alpha) * src.cols / 2 - beta * src.rows / 2,
             -beta, alpha, beta * src.cols / 2 + (1 - alpha) * src.rows / 2);

    Mat_<uchar> img = src;
    int c = 0;
    Mat det = Mat(cv::Size(src_w, src_h), CV_8UC1, Scalar(0));
    for (int y = 0; y < src.rows; ++y)
    { // src_h

        for (int x = 0; x < src.cols; ++x)
        { // src_w
            vector<int> yuan = affine(y, x, r);
            if (yuan[1] >= src.cols || yuan[0] >= src.rows || yuan[0] < 0 || yuan[1] < 0)
                continue;
            else
            {
                uchar p = img(y, x);
                det.at<uchar>(yuan[0], yuan[1]) = p;
            }
        }
    }
    return det;
}