#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <algorithm> // max
#include<cmath> // pow
#include"meshFast.h"
#include"ZNCC.h"
#include"CImg.h"
using namespace cimg_library;

std::vector<double> abs(std::vector<double> a) {
    std::vector<double> absA;
    for (int i = 0; i < a.size(); i++)
    {
        if (a[i] >= 0)
        {
            absA.push_back(a[i]);
        }
        else {
            absA.push_back(-a[i]);
        }
    }
    return absA;
}

int max(std::vector<int> a) {
    double maximum = a[0];
    for (int i = 1; i < a.size(); i++)
    {
        if (a[i] > maximum)
        {
            maximum = a[i];
        }
    }
    return maximum;
}


void opticalMethod(cimg_library::CImg<unsigned int> image1, cimg_library::CImg<unsigned int> image2, double* meshY, double* meshX, double* y_area_top,
    double* x_area_rig, double* y_area_bot, double* x_area_lef, double* Ny1, double* Nx1, int count) {

    double scale = 1.0;
    double e = 0.000000001;

    int it = 1;

    std::ofstream velocityFile;
    velocityFile.open("C:/Users/1/source/repos/Course work/PIV/velocity.txt"); // open file for writing velocity

    // std::vector<double> meshX; std::vector<double> meshY - x and y - coordinate of mesh points
    // std::vector<int> Nx1; std::vector<int> Ny1 - x and y - size of IW
    // std::vector<double> y_area_top - top-size of serching area
    // std::vector<double> x_area_rig - right-size of serching area
    // std::vector<double> y_area_bot - bottom-size of serching area
    // std::vector<double> x_area_lef - left-size of serching are

    //int count = sizeof(*Nx1) / sizeof(Nx1[0]);

    // int Nx_max = max(Nx1);
     //int Ny_max = max(Ny1);
    int* meshXint = meshFast(count, Nx1, meshX); // acceleration without interpolation
    int* meshYint = meshFast(count, Ny1, meshY); // acceleration without interpolation

    // END - number of points in x - and y - directions for median test

    //double* x{ new double[2]{ 0.0, 0.0 } };  double* y{ new double[2] { 0.0, 0.0 } };  double* dx{ new double[2] { 0.0, 0.0 } };  double* dy{ new double[2] { 0.0, 0.0 } };
    //y = x; dx = x; dy = x;
    //double* leftx = x; double* topy = x; // при интерпол€ции интенсивности
    double* uPredict{ new double[count] {} }; // prediction of x - displacement
    double* vPredict{ new double[count] {} }; // prediction of y - displacement
    double* u1{ new double[count] {} };
    double* v1{ new double[count] {} };
    double* u2{ new double[count] {} };
    double* v2{ new double[count] {} };
    double** SNR{ new double*[count] {} }; // ratio of 2 local maximus of ZNCC
    for (int i = 0; i < count; i++)
    {
        SNR[i] = new double[1] {};
    }

    double* Nx2 = Nx1; // x - size of IW
    double* Ny2 = Ny1; // y - size of IW

    double** im1{ new double* [image1.height()] {} }; double** im2{ new double* [image2.height()] {} }; // arrays of the whole images

    for (int i = 0; i < image1.height(); i++)
    {
        im1[i] = new double[image1.width()] {};
    }

    for (int i = 0; i < image2.height(); i++)
    {
        im2[i] = new double[image2.width()] {};
    }

    for (int i = 0; i < image1.height(); i++)
    {
        for (int j = 0; j < image1.width(); j++)
        {
            im1[i][j] = double(image1(j, i));
        }
    }

    for (int i = 0; i < image2.height(); i++)
    {
        for (int j = 0; j < image2.width(); j++)
        {
            im2[i][j] = double(image2(j, i));
        }
    }

    int Iby = image2.height();
    int Ibx = image2.width();

    double*** GSTable1 = GST(im1, Iby, Ibx);
    double*** GSTable2 = GST(im2, Iby, Ibx);

    // GlobalSums11 - GST for image 1
    // GlobalSums12 - GST for image 2
    // GlobalSums21 - GST of squares for image 1
    // GlobalSums22 - GST of squares image 2

    double** GlobalSums11{ new double* [Iby] {} }; // memory allocation
    double** GlobalSums12{ new double* [Iby] {} };
    double** GlobalSums21{ new double* [Iby] {} };
    double** GlobalSums22{ new double* [Iby] {} };

    for (int i = 0; i < Iby; i++)
    {
        GlobalSums11[i] = new double[Ibx] {};
        GlobalSums12[i] = new double[Ibx] {};
        GlobalSums21[i] = new double[Ibx] {};
        GlobalSums22[i] = new double[Ibx] {};
    }

    for (int i = 0; i < Iby; i++)
    {
        for (int j = 0; j < Ibx; j++)
        {
            GlobalSums11[i][j] = GSTable1[i][j][0];
            GlobalSums12[i][j] = GSTable2[i][j][0];
            GlobalSums21[i][j] = GSTable1[i][j][1];
            GlobalSums22[i][j] = GSTable2[i][j][1];
        }
    }

    double aver1 = 0.0; double D1 = 0.0; // the average and dispersion of the 1th image
    double aver2 = 0.0; double D2 = 0.0; // the average and dispersion of the 2th image

    aver1 = GlobalSums11[image1.height() - 1][image1.width() - 1] / (image1.height() * image1.width()) + e;
    aver2 = GlobalSums12[image2.height() - 2][image2.width() - 1] / (image2.height() * image2.width()) + e;
    D1 = sqrt(GlobalSums21[image1.height() - 1][image1.width() - 1] + image1.height() * image1.width() * aver1 * aver1) + e;
    D2 = sqrt(GlobalSums22[image2.height() - 1][image2.width() - 1] + image2.height() * image2.width() * aver2 * aver2) + e;


    // -BEGIN - first iteration of displacement calculation
    if (it == 1) {

        for (int q = 0; q < 1; q++)
        {
            std::cout << "number of iteration is " << it << "   number of frame is " << q << std::endl;

            int ii = q; // current image




            for (int p = 0; p < count; p++)
            {
                double x_lefAndRig[4]; double y_topAndBot[4];
                // x_lefAndRig[0] - left x - coordinate of 1st image
                // x_lefAndRig[1] - left x - coordinate of 2st image
                // x_lefAndRig[2] - right x - coordinate of 1st image
                // x_lefAndRig[3] - right x - coordinate of 2nd image

                // y_topAndBot[0] - top y - coordinate of 1st image
                // y_topAndBot[1] - top y - coordinate of 2st image
                // y_topAndBot[2] - bottom y - coordinate of 1st image
                // y_topAndBot[3] - bottom y - coordinate of 2nd image


                x_lefAndRig[0] = meshXint[p] - Nx2[p] / 2.0 + 0.5; // left x - coordinate of 1st image
                x_lefAndRig[1] = x_lefAndRig[0] - x_area_lef[p]; // left x - coordinate of 2nd image
                x_lefAndRig[2] = x_lefAndRig[0] + Nx2[p] - 1.0; // right x - coordinate of 1st image
                x_lefAndRig[3] = x_lefAndRig[2] + x_area_rig[p]; // right x - coordinate of 2nd image

                y_topAndBot[0] = meshYint[p] - Ny2[p] / 2.0 + 0.5; // top y - coordinate of 1st image
                y_topAndBot[1] = y_topAndBot[0] - y_area_top[p]; // top y - coordinate of 2st image
                y_topAndBot[2] = y_topAndBot[0] + Ny2[p] - 1.0; // bottom y - coordinate of 1st image
                y_topAndBot[3] = y_topAndBot[2] + y_area_bot[p]; // bottom y - coordinate of 2st image

                int Nxf1 = x_lefAndRig[2] - x_lefAndRig[0] + 1; int Nxf2 = x_lefAndRig[3] - x_lefAndRig[1] + 1;
                int Nyf1 = y_topAndBot[2] - y_topAndBot[0] + 1; int Nyf2 = y_topAndBot[3] - y_topAndBot[1] + 1;

                int l_border1 = (x_lefAndRig[0])*(copysign(1, x_lefAndRig[0]) - 1)*0.5;
                int r_border1 = (Nxf1 * 0.5 * (copysign(1, Ibx - x_lefAndRig[2]) +
                    copysign(1, Ibx - x_lefAndRig[2]) * copysign(1, Ibx - x_lefAndRig[2]))) + 0.5 * (Ibx - (int)abs(x_lefAndRig[0]) - 1) *
                    (copysign(1, x_lefAndRig[2] - Ibx) + copysign(1, x_lefAndRig[2] - Ibx) * copysign(1, x_lefAndRig[2] - Ibx));
                int t_border1 = y_topAndBot[0] * (copysign(1, y_topAndBot[0]) - 1) * 0.5;
                int b_border1 = (Nyf1 * 0.5 * (copysign(1, Iby - y_topAndBot[2]) + 
                    copysign(1, Iby - y_topAndBot[2]) * copysign(1, Iby - y_topAndBot[2]))) + 0.5 * (Iby - (int)abs(y_topAndBot[0]) - 1) * 
                    (copysign(1, y_topAndBot[2] - Iby) + copysign(1, y_topAndBot[2] - Iby) * copysign(1, y_topAndBot[2] - Iby));
                
                int l_border2 = (x_lefAndRig[1]) * (copysign(1, x_lefAndRig[1]) - 1) * 0.5;
                int r_border2 = (Nxf2 * 0.5 * (copysign(1, Ibx - x_lefAndRig[3]) +
                    copysign(1, Ibx - x_lefAndRig[3]) * copysign(1, Ibx - x_lefAndRig[3]))) + 0.5 * (Ibx - (int)abs(x_lefAndRig[1]) - 1) *
                    (copysign(1, x_lefAndRig[3] - Ibx) + copysign(1, x_lefAndRig[3] - Ibx) * copysign(1, x_lefAndRig[3] - Ibx));
                int t_border2 = y_topAndBot[1] * (copysign(1, y_topAndBot[1]) - 1) * 0.5;
                int b_border2 = (Nyf2 * 0.5 * (copysign(1, Iby - y_topAndBot[3]) +
                    copysign(1, Iby - y_topAndBot[3]) * copysign(1, Iby - y_topAndBot[3]))) + 0.5 * (Iby - (int)abs(y_topAndBot[1]) - 1) *
                    (copysign(1, y_topAndBot[3] - Iby) + copysign(1, y_topAndBot[3] - Iby) * copysign(1, y_topAndBot[3] - Iby));

                double** f1{ new double* [Nyf1] {} };
                double** f2{ new double* [Nyf2] {} };
                for (int i = 0; i < Nyf1; i++)
                {
                    f1[i] = new double[Nxf1] {};
                }
                for (int i = 0; i < Nyf2; i++)
                {
                    f2[i] = new double[Nxf2] {};
                }

                for (int i = 0; i < Nyf1; i++)
                {
                    for (int j = 0; j < Nxf1; j++)
                    {
                        f1[i][j] = 0.0;
                    }
                }
                for (int i = 0; i < Nyf2; i++)
                {
                    for (int j = 0; j < Nxf2; j++)
                    {
                        f2[i][j] = 0.0;
                    }
                }

                for (int i = t_border1; i < b_border1; i++)
                {
                    for (int j = l_border1; j < r_border1; j++)
                    {
                        f1[i][j] = im1[i + (int)y_topAndBot[0]][j + (int)x_lefAndRig[0]];
                    }
                }
                for (int i = t_border2; i < b_border2; i++)
                {
                    for (int j = l_border2; j < r_border2; j++)
                    {
                        f2[i][j] = im2[i + (int)y_topAndBot[1]][j + (int)x_lefAndRig[1]];
                    }
                }

                /*for (int i = (int)y_topAndBot[0]; i < (int)y_topAndBot[2] + 1; i++)
                //{
                    for (int j = (int)x_lefAndRig[0]; j < (int)x_lefAndRig[2] + 1; j++)
                    {
                         f1[i - (int)y_topAndBot[0]][j - (int)x_lefAndRig[0]] = im1[i][j];
                    }
                }
                for (int i = (int)y_topAndBot[1]; i < (int)y_topAndBot[3] + 1; i++)
                {
                    for (int j = (int)x_lefAndRig[1]; j < (int)x_lefAndRig[3] + 1; j++)
                    {
                        f2[i - (int)y_topAndBot[1]][j - (int)x_lefAndRig[1]] = im2[i][j];
                    }
                }*/

                /*for (int i = 0; i < Nyf1; i++)
                {
                    for (int j = 0; j < Nxf1; j++)
                    {
                        std::cout << f1[i][j] << "  ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl; std::cout << std::endl;
                for (int i = 0; i < Nyf2; i++)
                {
                    for (int j = 0; j < Nxf2; j++)
                    {
                        std::cout << f2[i][j] << "  ";
                    }
                    std::cout << std::endl;
                }*/

                /*aver1 = (GlobalSums11[(int)y_topAndBot[2] - 1][(int)x_lefAndRig[2] - 1] - GlobalSums11[(int)y_topAndBot[2] - 1][(int)x_lefAndRig[0] - 1] - GlobalSums11[(int)y_topAndBot[0] - 1][(int)x_lefAndRig[2] - 1] + GlobalSums11[(int)y_topAndBot[0] - 1][(int)x_lefAndRig[0] - 1]) / (Nyf1 * Nxf1) + e;
                aver2 = (GlobalSums12[(int)y_topAndBot[2] - 1][(int)x_lefAndRig[2] - 1] - GlobalSums12[(int)y_topAndBot[2] - 1][(int)x_lefAndRig[0] - 1] - GlobalSums12[(int)y_topAndBot[0] - 1][(int)x_lefAndRig[2] - 1] + GlobalSums12[(int)y_topAndBot[0] - 1][(int)x_lefAndRig[0] - 1]) / (Nyf1 * Nxf1) + e;
                D1 = sqrt((GlobalSums21[(int)y_topAndBot[2] - 1][(int)x_lefAndRig[2] - 1] - GlobalSums21[(int)y_topAndBot[2] - 1][(int)x_lefAndRig[0] - 1] - GlobalSums21[(int)y_topAndBot[0] - 1][(int)x_lefAndRig[2] - 1] + GlobalSums21[(int)y_topAndBot[0] - 1][(int)x_lefAndRig[0] - 1]) + Nyf1 * Nxf1 * aver1 * aver1) + e;
                D2 = sqrt((GlobalSums22[(int)y_topAndBot[2] - 1][(int)x_lefAndRig[2] - 1] - GlobalSums22[(int)y_topAndBot[2] - 1][(int)x_lefAndRig[0] - 1] - GlobalSums22[(int)y_topAndBot[0] - 1][(int)x_lefAndRig[2] - 1] + GlobalSums22[(int)y_topAndBot[0] - 1][(int)x_lefAndRig[0] - 1]) + Nyf1 * Nxf1 * aver2 * aver2) + e;*/

                /*aver1 = (GlobalSums11[b_border1 - 1 + (int)abs(y_topAndBot[0])][r_border1 - 1 + (int)abs(x_lefAndRig[0])] - GlobalSums11[b_border1 - 1 + (int)abs(y_topAndBot[0])][l_border1 - 1 + (int)abs(x_lefAndRig[0])] - GlobalSums11[t_border1 - 1 + (int)abs(y_topAndBot[0])][r_border1 - 1 + (int)abs(x_lefAndRig[0])] + GlobalSums11[t_border1 - 1 + (int)abs(y_topAndBot[0])][l_border1 - 1 + (int)abs(x_lefAndRig[0])]) / (Nyf1 * Nxf1) + e;
                aver2 = (GlobalSums12[b_border2 - 1 + (int)abs(y_topAndBot[1])][r_border2 - 1 + (int)abs(x_lefAndRig[1])] - GlobalSums12[b_border2 - 1 + (int)abs(y_topAndBot[1])][l_border2 - 1 + (int)abs(x_lefAndRig[1])] - GlobalSums12[t_border2 - 1 + (int)abs(y_topAndBot[1])][r_border2 - 1 + (int)abs(x_lefAndRig[1])] + GlobalSums12[t_border2 - 1 + (int)abs(y_topAndBot[1])][l_border2 - 1 + (int)abs(x_lefAndRig[1])]) / (Nyf1 * Nxf1) + e;
                D1 = sqrt((GlobalSums21[b_border1 - 1 + (int)abs(y_topAndBot[0])][r_border1 - 1 + (int)abs(x_lefAndRig[0])] - GlobalSums21[b_border1 - 1 + (int)abs(y_topAndBot[0])][l_border1 - 1 + (int)abs(x_lefAndRig[0])] - GlobalSums21[t_border1 - 1 + (int)abs(y_topAndBot[0])][r_border1 - 1 + (int)abs(x_lefAndRig[0])] + GlobalSums21[t_border1 - 1 + (int)abs(y_topAndBot[0])][l_border1 - 1 + (int)abs(x_lefAndRig[0])]) + Nyf1 * Nxf1 * aver1 * aver1) + e;
                D2 = sqrt((GlobalSums22[b_border2 - 1 + (int)abs(y_topAndBot[1])][r_border2 - 1 + (int)abs(x_lefAndRig[1])] - GlobalSums22[b_border2 - 1 + (int)abs(y_topAndBot[1])][l_border2 - 1 + (int)abs(x_lefAndRig[1])] - GlobalSums22[t_border2 - 1 + (int)abs(y_topAndBot[1])][r_border2 - 1 + (int)abs(x_lefAndRig[1])] + GlobalSums22[t_border2 - 1 + (int)abs(y_topAndBot[1])][l_border2 - 1 + (int)abs(x_lefAndRig[1])]) + Nyf1 * Nxf1 * aver2 * aver2) + e;*/

                double* ZNCC = pivZNCC(f1, f2, ((int)x_lefAndRig[0] - (int)x_lefAndRig[1]), ((int)y_topAndBot[0] - (int)y_topAndBot[1]), Nxf1, Nxf2, Nyf1, Nyf2, aver1, aver2, D1, D2);
                u1[p] = ZNCC[0]; v1[p] = ZNCC[1]; SNR[p][ii] = ZNCC[2];

                for (int i = 0; i < Nyf1; ++i) 
                {
                    delete[] f1[i];
                }                  
                delete[] f1;
                for (int i = 0; i < Nyf2; ++i)
                {
                    delete[] f2[i];
                }
                delete[] f2;
            }


            for (int i = 0; i < count; i++)
            {
                u2[i] = uPredict[i] + u1[i];
                v2[i] = vPredict[i] + v1[i];
            }

            //[u2, v2, error(:, 1)] = median_test([s_out name], u2, v2, error(:, 1), 2.0, 1 * bx, 1 * by);

        }


        for (int i = 0; i < count; i++)
        {
            velocityFile << u2[i] << "  " << v2[i] << std::endl; // writing x and y shifts to a file
        }
        velocityFile << std::endl;
        for (int i = 0; i < count; i++)
        {
            for (int j = 0; j < 1; j++)
            {
                velocityFile << SNR[i][j] << "  ";
            }
            velocityFile << std::endl;
        }
        velocityFile << std::endl;

        it = it + 1;

    }

    // -END - first iteration of displacement calculation

    delete [] meshXint; delete [] meshYint;
   /* delete [] x; delete [] y; delete [] dx; delete [] dy;
    delete [] leftX; delete [] topY;*/
    delete [] uPredict; delete [] vPredict;
    delete [] u1; delete [] v1; delete [] u2; delete [] v2;
    for (int i = 0; i < count; ++i)
    {
        delete[] SNR[i];
    }
    delete[] SNR;
    delete[] Nx2; delete[] Ny2;
    for (int i = 0; i < image1.height(); ++i)
    {
        delete[] im1[i];
    }
    delete[] im1;
    for (int i = 0; i < image2.height(); ++i)
    {
        delete[] im2[i];
    }
    delete[] im2;
    for (int i = 0; i < Iby; ++i)
    {
        delete[] GlobalSums11[i];
        delete[] GlobalSums12[i];
        delete[] GlobalSums21[i];
        delete[] GlobalSums22[i];
    }
    delete[] GlobalSums11;
    delete[] GlobalSums12;
    delete[] GlobalSums21;
    delete[] GlobalSums22;
    for (int i = 0; i < Iby; ++i) {
        for (int j = 0; j < Ibx; ++j) 
        {
            delete[] GSTable1[i][j];
            delete[] GSTable2[i][j];
        }
        delete[] GSTable1[i];
        delete[] GSTable2[i];
    }
    delete[] GSTable1;
    delete[] GSTable2;

    velocityFile.close();

}