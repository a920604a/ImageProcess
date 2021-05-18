

Mat scale(Mat src, int dst_h, int dst_w)
{

    Mat det = Mat(cv::Size(dst_h, dst_w), CV_8UC1, Scalar(0));

    int src_h = src.rows, src_w = src.cols;
    float scale_w = (float)dst_w / src_w, scale_h = (float)dst_h / src_h;
    double inf_w = 1.0 / scale_w, inf_h = 1.0 / scale_h;

    for (int y = 0; y < dst_w; ++y)
    {
        uchar *d = src.ptr<uchar>(y / scale_w);
        uchar *e = src.ptr<uchar>((y + 1) / scale_w);
        for (int x = 0; x < dst_h; ++x)
        {

            float src_x = x * (float)src_h / dst_h;
            float src_y = y * (float)src_w / dst_w;

            int src_x_int = floor(src_x);
            int src_y_int = floor(src_y);

            float src_x_float = src_x - src_x_int;
            float src_y_float = src_y - src_y_int;
            if (src_x_int + 1 == src_w || src_y_int + 1 == src_h)
            {
                det.at<uchar>(y, x) = d[(int)(x / scale_h)];
                continue;
            }
            det.at<uchar>(y, x) =
                (1. - src_y_float) * (1. - src_x_float) * d[(int)(x / scale_h)] +
                (1. - src_y_float) * src_x_float * d[(int)((x + 1) / scale_h)] +
                src_y_float * (1. - src_x_float) * e[(int)(x / scale_h)] +
                src_y_float * src_x_float * e[(int)((x + 1) / scale_h)];
        }
    }


    return det;
}