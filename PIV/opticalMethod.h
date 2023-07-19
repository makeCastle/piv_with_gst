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


void opticalMethod(cimg_library::CImg<unsigned int> image1, cimg_library::CImg<unsigned int> image2, std::vector<double> meshY, std::vector<double> meshX, std::vector<double> y_area_top,
    std::vector<double> x_area_rig, std::vector<double> y_area_bot, std::vector<double> x_area_lef, std::vector<double> Ny1, std::vector<double> Nx1) {

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

    int count = Nx1.size();

    // int Nx_max = max(Nx1);
     //int Ny_max = max(Ny1);
    std::vector<int> meshXint = meshFast(count, Nx1, meshX); // acceleration without interpolation
    std::vector<int> meshYint = meshFast(count, Ny1, meshY); // acceleration without interpolation

    // END - number of points in x - and y - directions for median test

    //int rad1[2] = { (meshXint[1] - meshXint[0]), (meshYint[Nx] - meshYint[0]) }; // дл€ увеличени€ окрестности поиска
    //int rad2[2] = { 1, 1 }; // узла дл€ интерпол€ции geedfree

    std::vector<double> x = { 0.0, 0.0 }; std::vector<double> y = x; std::vector<double> dx = x; std::vector<double> dy = x;
    std::vector<double> leftX = x; std::vector<double> topY = x; // при интерпол€ции интенсивности
    std::vector<double> uPredict(count, 0.0); // prediction of x - displacement
    std::vector<double> vPredict(count, 0.0); // prediction of y - displacement
    std::vector<double> u1;
    std::vector<double> v1;
    std::vector<double> u2;
    std::vector<double> v2;
    std::vector<std::vector<double>> SNR; // ratio of 2 local maximus of ZNCC
    for (int i = 0; i < count; i++) {
        std::vector<double> row(1, 0.0);
        SNR.push_back(row);
    }

    std::vector<double> Nx2 = Nx1; // x - size of IW
    std::vector<double> Ny2 = Ny1; // y - size of IW

    std::vector<std::vector<double>> im1; std::vector<std::vector<double>> im2; // arrays of the whole images
    for (int i = 0; i < image1.height(); i++)
    {
        std::vector<double> row1;
        for (int j = 0; j < image1.width(); j++)
        {
            row1.push_back(double(image1(j, i)));
        }
        im1.push_back(row1);
    }

    for (int i = 0; i < image2.height(); i++)
    {
        std::vector<double> row2;
        for (int j = 0; j < image2.width(); j++)
        {
            row2.push_back(double(image2(j, i)));
        }
        im2.push_back(row2);
    }

    int Iby = image2.height();
    int Ibx = image2.width();

    std::vector<std::vector<std::vector<double>>> GSTable1 = GST(im1, Iby, Ibx);
    std::vector<std::vector<std::vector<double>>> GSTable2 = GST(im2, Iby, Ibx);

    // GlobalSums11 - GST for image 1
    // GlobalSums12 - GST for image 2
    // GlobalSums21 - GST of squares for image 1
    // GlobalSums22 - GST of squares image 2

    std::vector<std::vector<double>> GlobalSums11;
    std::vector<std::vector<double>> GlobalSums12;
    std::vector<std::vector<double>> GlobalSums21;
    std::vector<std::vector<double>> GlobalSums22;

    GlobalSums11.resize(Iby); // memory allocation
    GlobalSums12.resize(Iby);
    GlobalSums21.resize(Iby);
    GlobalSums22.resize(Iby);
    for (int i = 0; i < Iby; ++i) {
        GlobalSums11[i].resize(Ibx);
    }
    for (int i = 0; i < Iby; ++i) {
        GlobalSums12[i].resize(Ibx);
    }
    for (int i = 0; i < Iby; ++i) {
        GlobalSums21[i].resize(Ibx);
    }
    for (int i = 0; i < Iby; ++i) {
        GlobalSums22[i].resize(Ibx);
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

    // -BEGIN - first iteration of displacement calculation
    if (it == 1) {

        for (int i = 0; i < 1; i++)
        {
            std::cout << "number of iteration is " << it << "   number of frame is " << i << std::endl;

            int ii = i; // current image




            for (int p = 0; p < count; p++)
            {
                std::vector<double> x_lefAndRig; std::vector<double> y_topAndBot;
                // x_lefAndRig[0] - left x - coordinate of 1st image
                // x_lefAndRig[1] - left x - coordinate of 2st image
                // x_lefAndRig[2] - right x - coordinate of 1st image
                // x_lefAndRig[3] - right x - coordinate of 2nd image

                // y_topAndBot[0] - top y - coordinate of 1st image
                // y_topAndBot[1] - top y - coordinate of 2st image
                // y_topAndBot[2] - bottom y - coordinate of 1st image
                // y_topAndBot[3] - bottom y - coordinate of 2nd image

                x_lefAndRig.push_back(meshXint[p] - Nx2[p] / 2.0 + 0.5); // left x - coordinate of 1st image
                x_lefAndRig.push_back(x_lefAndRig[0] - x_area_lef[p]); // left x - coordinate of 2nd image
                x_lefAndRig.push_back(x_lefAndRig[0] + Nx2[p] - 1.0); // right x - coordinate of 1st image
                x_lefAndRig.push_back(x_lefAndRig[2] + x_area_rig[p]); // right x - coordinate of 2nd image

                y_topAndBot.push_back(meshYint[p] - Ny2[p] / 2.0 + 0.5); // top y - coordinate of 1st image
                y_topAndBot.push_back(y_topAndBot[0] - y_area_top[p]); // top y - coordinate of 2st image
                y_topAndBot.push_back(y_topAndBot[0] + Ny2[p] - 1.0); // bottom y - coordinate of 1st image
                y_topAndBot.push_back(y_topAndBot[2] + y_area_bot[p]); // bottom y - coordinate of 2st image

                std::vector<std::vector<double>> f1;
                std::vector<std::vector<double>> f2;

                int Nxf1 = x_lefAndRig[2] - x_lefAndRig[0] + 1; int Nxf2 = x_lefAndRig[3] - x_lefAndRig[1] + 1;
                int Nyf1 = y_topAndBot[2] - y_topAndBot[0] + 1; int Nyf2 = y_topAndBot[3] - y_topAndBot[1] + 1;



                for (int i = (int)y_topAndBot[0]; i < (int)y_topAndBot[2] + 1; i++)
                {
                    std::vector<double> rowf1;
                    for (int j = (int)x_lefAndRig[0]; j < (int)x_lefAndRig[2] + 1; j++)
                    {
			if ((i > image1.height() || i < 0) || (j > image1.width() || j < 0)) {
				rowf1.push_back(0.0);
			} else { 
                        	rowf1.push_back(im1[i][j]);
			}
                    }
                    f1.push_back(rowf1);
                }
                for (int i = (int)y_topAndBot[1]; i < (int)y_topAndBot[3] + 1; i++)
                {
                    std::vector<double> rowf2;
                    for (int j = (int)x_lefAndRig[1]; j < (int)x_lefAndRig[3] + 1; j++)
                    {
			if ((i > image2.height() || i < 0) || (j > image2.width() || j < 0)) {
				rowf2.push_back(0.0);
			} else { 
                        	rowf2.push_back(im2[i][j]);
			}
                    }
                    f2.push_back(rowf2);
                }

                aver1 = (GlobalSums11[y_topAndBot[2] - 1][x_lefAndRig[2] - 1] - GlobalSums11[y_topAndBot[2] - 1][x_lefAndRig[0] - 1] - GlobalSums11[y_topAndBot[0] - 1][x_lefAndRig[2] - 1] + GlobalSums11[y_topAndBot[0] - 1][x_lefAndRig[0] - 1]) / (Nyf1 * Nxf1) + e;
                aver2 = (GlobalSums12[y_topAndBot[2] - 1][x_lefAndRig[2] - 1] - GlobalSums12[y_topAndBot[2] - 1][x_lefAndRig[0] - 1] - GlobalSums12[y_topAndBot[0] - 1][x_lefAndRig[2] - 1] + GlobalSums12[y_topAndBot[0] - 1][x_lefAndRig[0] - 1]) / (Nyf1 * Nxf1) + e;
                D1 = sqrt((GlobalSums21[y_topAndBot[2] - 1][x_lefAndRig[2] - 1] - GlobalSums21[y_topAndBot[2] - 1][x_lefAndRig[0] - 1] - GlobalSums21[y_topAndBot[0] - 1][x_lefAndRig[2] - 1] + GlobalSums21[y_topAndBot[0] - 1][x_lefAndRig[0] - 1]) + Nyf1 * Nxf1 * aver1 * aver1) + e;
                D2 = sqrt((GlobalSums22[y_topAndBot[2] - 1][x_lefAndRig[2] - 1] - GlobalSums22[y_topAndBot[2] - 1][x_lefAndRig[0] - 1] - GlobalSums22[y_topAndBot[0] - 1][x_lefAndRig[2] - 1] + GlobalSums22[y_topAndBot[0] - 1][x_lefAndRig[0] - 1]) + Nyf1 * Nxf1 * aver2 * aver2) + e;


                std::vector<double> ZNCC = pivZNCC(f1, f2, (x_lefAndRig[0] - x_lefAndRig[1]), (y_topAndBot[0] - y_topAndBot[1]), Nxf1, Nxf2, Nyf1, Nyf2, aver1, aver2, D1, D2);
                u1.push_back(ZNCC[0]); v1.push_back(ZNCC[1]); SNR[p][ii] = ZNCC[2];
            }


            for (int i = 0; i < count; i++)
            {
                u2.push_back(uPredict[i] + u1[i]);
                v2.push_back(vPredict[i] + v1[i]);
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

    velocityFile.close();

}