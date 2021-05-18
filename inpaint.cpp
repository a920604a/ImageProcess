/*
 * @Author: yuan
 * @Date: 2021-05-13 10:18:01
 * @LastEditTime: 2021-05-18 09:27:42
 * @FilePath: /面試題目_20210121/inpaint.cpp
 */

#include <math.h>
#include <vector>
#include <queue>
#include <float.h>

struct point
{
    float x;
    float y;
};

struct Item
{
    int i = 0;
    int j = 0;
    float p = 0;

    Item(int i, int j, int p) : i(i), j(j), p(p){};
};

class Compare
{
public:
    bool operator()(const Item &lhs, const Item &rhs)
    {
        if (lhs.p == rhs.p)
        {
            return lhs.i > rhs.i;
        }
        return lhs.p < rhs.p;
    }
};

template <typename T>
bool sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

void mask(Mat &src, Mat &dst, int msize, int center_i, int center_j)
{

    src.copyTo(dst);
    int left = center_j - msize / 2 - 1, top = center_i - msize / 2 - 1;

    for (int i = left; i <= left + msize; ++i)
    {
        for (int j = top; j <= top + msize; ++j)
        {
            dst.at<uchar>(i, j) = 0;
        }
    }
}

void Mask(Mat &src, Mat &dst, int msize, int center_i, int center_j)
{
    src.copyTo(dst);

    int left = center_j - msize / 2 - 1, top = center_i - msize / 2 - 1;

    for (int i = 0; i < src.rows; ++i)
    {
        for (int j = 0; j < src.cols; ++j)
        {
            if (i >= left && i <= left + msize && j >= top && j <= top + msize)
                dst.at<uchar>(i, j) = 255;
            else
                dst.at<uchar>(i, j) = 0;
        }
    }
}

static float VectorScalMult(point &v1, point &v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

static float VectorLength(point &v1)
{
    return sqrt(v1.x * v1.x + v1.y * v1.y);
}

float solver(int i1, int j1, int i2, int j2, vector<vector<char>> &visit, vector<vector<double>> &t)
{

    float sol = exp(6);
    if (visit[i1][j1] == 'N')
    {
        if (visit[i2][j2] == 'N')
        {
            float r = sqrt(2 - ((t[i1][j1] - t[i2][j2]) * (t[i1][j1] - t[i2][j2])));
            float s = (t[i1][j1] + t[i2][j2] - r) / 2;
            if (s >= t[i1][j1] && s >= t[i2][j2])
                sol = s;
            else
            {
                s += r;
                if (s >= t[i1][j1] && s >= t[i2][j2])
                    sol = s;
            }
        }
        else
            sol = 1 + t[i1][j1];
    }
    else if (visit[i2][j2] == 'N')
        sol = 1 + t[i2][j2];

    return sol;
}

void paint(int i, int j, Mat &out, vector<vector<char>> &visit, vector<vector<double>> &t)
{
    // search neighborhood(l,k) of (i,j) such that not in outside

    point gradT, gradI;
    if (visit[i][j + 1] != 'I')
    {
        if (visit[i][j - 1] != 'I')
        {
            gradT.x = t[i][j + 1] - t[i][j - 1];
        }
        else
        {
            gradT.x = t[i][j + 1] - t[i][j];
        }
    }
    else
    {
        if (visit[i][j - 1] != 'I')
        {
            gradT.x = t[i][j] - t[i][j - 1];
        }
        else
        {
            gradT.x = 0;
        }
    }

    if (visit[i + 1][j] != 'I')
    {
        if (visit[i - 1][j] != 'I')
        {
            gradT.y = t[i + 1][j] - t[i - 1][j];
        }
        else
        {
            gradT.y = t[i + 1][j] - t[i][j];
        }
    }
    else
    {
        if (visit[i - 1][j] != 'I')
        {
            gradT.y = t[i][j] - t[i - 1][j];
        }
        else
        {
            gradT.y = 0;
        }
    }
    int range = 1;
    point r;
    float Ia = 0, Jx = 0, Jy = 0, s = 1.0e-20f, w, dst, lev, dir, sat;
    for (int k = i - range; k <= i + range; k++)
    {
        int km = k - 1 + (k == 1), kp = k - 1 - (k == t.size() - 2);
        for (int l = j - range; l <= j + range; l++)
        {
            int lm = l - 1 + (l == 1), lp = l - 1 - (l == t[0].size() - 2);
            if (k > 0 && l > 0 && k < t.size() - 1 && l < t[0].size() - 1)
            {
                if ((visit[k][l] != 'I') && ((l - j) * (l - j) + (k - i) * (k - i) <= range * range))
                {
                    r.x = (float)(i - k);
                    r.y = (float)(j - l);

                    dir = VectorScalMult(r, gradT);
                    if (fabs(dir) <= 0.01)
                        dir = 0.000001f;
                    dst = 1. / (VectorLength(r) * VectorLength(r));
                    // dst = (float)(1. / (VectorLength(r) * sqrt(VectorLength(r))));
                    lev = (float)(1. / (1 + fabs(t[k][l]) - t[i][j]));

                    w = (float)fabs(dst * lev * dir);

                    if (visit[k][l + 1] != 'I')
                    {
                        if (visit[k][l - 1] != 'I')
                        {
                            gradI.x = (float)(out.at<uchar>(km, lp + 1) - out.at<uchar>(km, lm - 1)) * 2.0f;
                        }
                        else
                        {
                            gradI.x = (float)(out.at<uchar>(km, lp + 1) - out.at<uchar>(km, lm));
                        }
                    }
                    else
                    {
                        if (visit[k][l - 1] != 'I')
                        {
                            gradI.x = (float)(out.at<uchar>(km, lp) - out.at<uchar>(km, lm - 1));
                        }
                        else
                        {
                            gradI.x = 0;
                        }
                    }
                    if (visit[k + 1][l] != 'I')
                    {
                        if (visit[k - 1][l] != 'I')
                        {
                            gradI.y = (float)(out.at<uchar>(kp + 1, lm) - out.at<uchar>(km - 1, lm)) * 2.0f;
                        }
                        else
                        {
                            gradI.y = (float)(out.at<uchar>(kp + 1, lm) - out.at<uchar>(km, lm));
                        }
                    }
                    else
                    {
                        if (visit[k - 1][l] != 'I')
                        {
                            gradI.y = (float)(out.at<uchar>(kp, lm) - out.at<uchar>(km - 1, lm));
                        }
                        else
                        {
                            gradI.y = 0;
                        }
                    }

                    Ia += (float)w * (float)(out.at<uchar>(km, lm));
                    Jx -= (float)w * (float)(gradI.x * r.x);
                    Jy -= (float)w * (float)(gradI.y * r.y);
                    s += w;
                }
            }
        }
    }
    sat = (float)((Ia / s + (Jx + Jy) / (sqrt(Jx * Jx + Jy * Jy) + 1.0e-20f) + 0.5f));
    {
        out.at<uchar>(i, j) = uchar(sat);
    }
}

Mat _dilate(Mat mask)
{
    Mat ret = Mat(mask.size(), mask.type(), Scalar(0));
    vector<vector<int>> dilation_matrix = {
        {0, 1, 0},
        {1, 1, 1},
        {0, 1, 0}};

    vector<pair<int, int>> coor;
    int sum, x, y, tmp;
    for (int i = 0; i < dilation_matrix.size(); ++i)
    {
        for (int j = 0; j < dilation_matrix[0].size(); ++j)
        {
            if (dilation_matrix[i][j] == 1)
                coor.push_back(make_pair(i, j));
        }
    }
    for (int i = 1; i < mask.rows - 1; ++i)
    {
        uchar *d = mask.ptr<uchar>(i);
        for (int j = 1; j < mask.cols - 1; ++j)
        {
            sum = 0;
            // kernel element multiply
            for (int k = 0; k < coor.size(); ++k)
            {
                auto p = coor[k];
                x = p.first;
                y = p.second;
                tmp = dilation_matrix[x][y] * (int)mask.at<uchar>(i + x, j + k);
                sum = max(tmp, sum);
            }

            ret.at<uchar>(i, j) = (uchar)sum;
        }
    }

    return ret;
}

Mat _subtract(Mat mask_edge, Mat mask)
{

    Mat ret = Mat(mask.size(), mask.type(), Scalar(0));

    int x = 0;
    for (int i = 0; i < mask.rows; ++i)
    {
        uchar *d = mask.ptr<uchar>(i);
        uchar *e = mask_edge.ptr<uchar>(i);
        for (int j = 0; j < mask.cols; ++j)
        {
            x = max((int)e[j] - (int)d[j], 0);
            ret.at<uchar>(i, j) = (uchar)x;
        }
    }

    return ret;
}

void inpaint(Mat &src, Mat &mask, Mat &out)
{
    src.copyTo(out);

    priority_queue<Item, vector<Item>, Compare> omega;

    int r = mask.rows, c = mask.cols;

    // initalize  matrix t, visit
    vector<vector<char>> visit(r, vector<char>(c, 'N')); //  B:band  N:known  I:inside
    vector<vector<double>> t(r, vector<double>(c, 1.0e6f));
    for (int i = 0; i < mask.rows; ++i)
    {
        for (int j = 0; j < mask.cols; ++j)
        {
            visit[i][j] = 'N';
            t[i][j] = 1.0e6f;
        }
    }

    vector<pair<int, int>> mask_coor;
    vector<pair<int, int>> msk;

    Mat mask_dilation = _dilate(mask);

    Mat edge = _subtract(mask_dilation, mask);
    for (int i = 0; i < mask.rows; ++i)
    {
        uchar *d = mask.ptr<uchar>(i);
        uchar *e = edge.ptr<uchar>(i);

        for (int j = 0; j < mask.cols; ++j)
        {
            if ((int)e[j] == 255)
            {
                omega.push(Item(i, j, FLT_MAX));
                visit[i][j] = 'B';
                t[i][j] = 0;
            }
            if ((int)d[j] == 255)
            {
                visit[i][j] = 'I';
                t[i][j] = 0;
            }
        }
    }

    int x = 0, y = 0, i = 0, j = 0;
    float dist = 0;
    while (!omega.empty())
    {
        auto p = omega.top();
        x = p.i;
        y = p.j;
        omega.pop();
        visit[x][y] = 'N';
        // traverse x,y surrounding
        for (int q = 0; q < 4; q++)
        {
            if (q == 0)
            {
                i = x - 1;
                j = y;
            }
            else if (q == 1)
            {
                i = x;
                j = y - 1;
            }
            else if (q == 2)
            {
                i = x + 1;
                j = y;
            }
            else if (q == 3)
            {
                i = x;
                j = y + 1;
            }
            if ((i <= 0) || (j <= 0) || (i > out.rows) || (j > out.cols))
                continue;

            if (visit[i][j] != 'N')
            {
                if (visit[i][j] == 'I')
                {

                    paint(i, j, out, visit, t);
                    visit[i][j] = 'B';
                }

                dist = min(
                    min(
                        solver(i - 1, j, i, j - 1, visit, t),
                        solver(i + 1, j, i, j - 1, visit, t)),
                    min(
                        solver(i - 1, j, i, j + 1, visit, t),
                        solver(i + 1, j, i, j + 1, visit, t)));
                t[i][j] = dist;
                omega.push(Item(i, j, dist));
            }
        }
    }
}

Mat mask_inpaint(Mat src)
{

    int center_i = src.rows / 2, center_j = src.cols / 2;

    Mat dst;
    mask(src, dst, 32, center_i, center_j);
    imwrite("mask.bmp", dst);
    Mat mask;
    Mask(src, mask, 32, center_i, center_j);
    imwrite("Mask.bmp", mask);

    Mat ret;
    inpaint(dst, mask, ret);
    return ret;
}
