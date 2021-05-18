

#include <iterator>
#include <random>
#include <tuple>

#define KERNEL_SIZE 3

Mat applyFilter(Mat src, Mat kernel)
{
    Mat dst(src.size(), src.type());

    for (int i = 0; i < src.cols; i++)
    {
        for (int j = 0; j < src.rows; j++)
        {
            for (int k = i; k < i + kernel.cols; k++)
            {
                for (int l = j; l < j + kernel.rows; l++)
                {

                    dst.at<uchar>(i, j) += kernel.at<double>(k - i, l - j) * (int)src.at<uchar>(k, l);
                }
            }
        }
    }

    return dst;
}

void generateGaussianTemplate(double window[3][3], int ksize, double sigma)
{
    int center = ksize / 2; // kernel的中心位置，也就是坐标的原点
    double x2, y2;
    double s = 2 * sigma * sigma;
    double sum = 0;
    for (int i = 0; i < ksize; i++)
    {
        x2 = pow(i - center, 2);
        for (int j = 0; j < ksize; j++)
        {
            y2 = pow(j - center, 2);
            double g = exp(-(x2 + y2) / s);
            g /= s * PI;
            window[i][j] = g;
            sum += g;
        }
    }
    for (int i = 0; i < ksize; i++)
    {
        for (int j = 0; j < ksize; j++)
        {
            window[i][j] /= sum;
        }
    }
}

void fill_guassian(Mat &src, Mat &dst, const double mean, const double stddev)
{
    src.copyTo(dst);
    // double g  = 1/(sigma * sqrt( 2* PI)) * exp( exp, ( - pow((x-sigma), 2) / 2*pow(sigma,2)) );

    // Define random generator with Gaussian distribution

    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);

    // Add Gaussian noise
    for (int y = 0; y < dst.rows; ++y)
    {
        uchar *d = dst.ptr<uchar>(y);
        for (int x = 0; x < dst.cols; ++x)
        {
            if ((int)d[x] + dist(generator) < 0 || (int)d[x] + dist(generator) > 255)
                continue;
            else
                dst.at<uchar>(y, x) = d[x] + (uchar)dist(generator);
        }
    }
}
struct noise_type
{
    Mat src;
    Mat dst;
};

noise_type noise(Mat src)
{

    // https://stackoverflow.com/questions/32889309/adding-gaussian-noise
    Mat dst;
    fill_guassian(src, dst, 2, 20);

    double window[KERNEL_SIZE][KERNEL_SIZE] = {0};
    // https: //www.cnblogs.com/wangguchangqing/p/6407717.html
    generateGaussianTemplate(window, KERNEL_SIZE, 0.8);

    Mat windows = Mat(KERNEL_SIZE, KERNEL_SIZE, cv::DataType<double>::type, window);

    // https://gist.github.com/OmarAflak/aca9d0dc8d583ff5a5dc16ca5cdda86a
    Mat det = applyFilter(dst, windows);

    return noise_type{dst, det};
}
