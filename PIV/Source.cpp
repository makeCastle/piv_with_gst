#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <ctime>      
#include"CImg.h"
#include"C:\Users\1\source\repos\Course work\PIV\meshFast.h"
#include"C:\Users\1\source\repos\Course work\PIV\opticalMethod.h"

int main() {
   /* cimg_library::CImg<unsigned char> image1("C:/Users/1/OneDrive/Документы/информатика, математика и физика/Курсовая/1.bmp");
    //std::cout << (int)image1(6, 81);
    //image1.display();
	std::ofstream fout;
	fout.open("Текст.txt");
	std::vector<std::vector<double>> a = { {5, 3, 1}, {6, 7, 8}, {9, 0, 2} };
	std::vector<std::vector<double>> b = { {6, 7, 1}, {3, 5, 9}, {4, 1, 0} };

	std::vector<std::vector<double>> c = { {3, 5, 9}, {-6, -7, -8}};
	std::vector<std::vector<double>> d = { {2, 0, 6}, {-3, -5, -9}};

    std::vector<double> z1(5, 0.0);  std::vector< std::vector<double>> a1;
    for (int i = 0; i < 1; i++)
    {
        a1.push_back(z1);
    }

    for (int i = 1; i < 5 - 1; ++i)
    {
        std::vector<double> row1;
        for (int j = 0; j < 1; j++)
        {
            row1.push_back(0);
        }
        for (int j = 1; j < 5 - 1; j++)
        {
            row1.push_back(a[i - 1][j - 1]);
        }
        for (int j = 5 - 1; j < 5; j++)
        {
            row1.push_back(0);
        }
        a1.push_back(row1);
    }

    for (int i = 0; i < 1; i++)
    {
        a1.push_back(z1);
    }

    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            std::cout << a1[i][j];
        }
        std::cout << std::endl;
    }

	fout.close(); */
    clock_t start, finish;
    start = clock();

    cimg_library::CImg<unsigned int> image1("D:/Tests/test files/photos/1.bmp");
    cimg_library::CImg<unsigned int> image2("D:/Tests/test files/photos/2.bmp");

    std::ifstream meshFile;
    meshFile.open("D:/Tests/test files/100000 vectors/mesh_test 100000 v 32x32 mesh 16x16 roi.txt"); // open file with mesh
    std::ifstream roiFile;
    roiFile.open("D:/Tests/test files/100000 vectors/ROI_test 100000 v 16x16.txt"); // open file with ROI-sizes

    double z = 0.0;

    std::vector<double> meshX; std::vector<double> meshY; // x and y - coordinate of mesh points
    std::vector<double> Nx1; std::vector<double> Ny1; // x and y - size of IW
    std::vector<double> y_area_top; // top-size of serching area
    std::vector<double> x_area_rig; // right-size of serching area
    std::vector<double> y_area_bot; // bottom-size of serching area
    std::vector<double> x_area_lef; // left-size of serching are

    while (!meshFile.eof()) // filling arrays of mesh
    {
        for (int i = 0; i < 5; i++)
        {
            meshFile >> z;
            if (i == 1)
            {
                meshX.push_back(z);
            }
            else if (i == 2)
            {
                meshY.push_back(z);
            }
            else if (i == 3)
            {
                Nx1.push_back(z);
            }
            else if (i == 4)
            {
                Ny1.push_back(z);
            }
        }
    }
    while (!roiFile.eof()) // filling arrays of ROI-sizes
    {
        for (int i = 0; i < 4; i++)
        {
            roiFile >> z;
            if (i == 0)
            {
                y_area_top.push_back(z);
            }
            else if (i == 1)
            {
                x_area_rig.push_back(z);
            }
            else if (i == 2)
            {
                y_area_bot.push_back(z);
            }
            else if (i == 3)
            {
                x_area_lef.push_back(z);
            }
        }
    }

    opticalMethod(image1, image2, meshY, meshX, y_area_top,x_area_rig, y_area_bot, x_area_lef, Ny1, Nx1);

    meshFile.close();
    roiFile.close();

    finish = clock();
    std::cout << "runtime = " << (double)(finish - start) / CLOCKS_PER_SEC << std::endl; // время работы программы  

    return 0;
}

