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

    cimg_library::CImg<unsigned int> image1("D:/Tests/Область поиска выходит за пределы/image51.bmp");
    cimg_library::CImg<unsigned int> image2("D:/Tests/Область поиска выходит за пределы/image52.bmp");

    std::ifstream meshFile;
    meshFile.open("D:/Tests/Область поиска выходит за пределы/mesh_test5.txt"); // open file with mesh
    std::ifstream roiFile;
    roiFile.open("D:/Tests/Область поиска выходит за пределы/ROI_test.txt"); // open file with ROI-sizes

    double z = 0.0; int k = 0; int c = 0;

    while (!meshFile.eof()) // filling arrays of mesh
    {
        for (int i = 0; i < 5; i++)
        {
            meshFile >> z;
        }
        ++k;
    }

    meshFile.close();
    meshFile.open("D:/Tests/Область поиска выходит за пределы/mesh_test5.txt");

    double* meshX{ new double[k] {} }; double* meshY{ new double[k] {} }; // x and y - coordinate of mesh points
    double* Nx1{ new double[k] {} }; double* Ny1{ new double[k] {} }; // x and y - size of IW
    double* y_area_top{ new double[k] {} }; // top-size of serching area
    double* x_area_rig{ new double[k] {} }; // right-size of serching area
    double* y_area_bot{ new double[k] {} }; // bottom-size of serching area
    double* x_area_lef{ new double[k] {} }; // left-size of serching are

    while (!meshFile.eof()) // filling arrays of mesh
    {
        for (int i = 0; i < 5; i++)
        {
            meshFile >> z;
            if (i == 1)
            {
                meshX[c] = z;
            }
            else if (i == 2)
            {
                meshY[c] = z;
            }
            else if (i == 3)
            {
                Nx1[c] = z;
            }
            else if (i == 4)
            {
                Ny1[c] = z;
            }
        }
        ++c;
    }
    c = 0;
    while (!roiFile.eof()) // filling arrays of ROI-sizes
    {
        for (int i = 0; i < 4; i++)
        {
            roiFile >> z;
            if (i == 0)
            {
                y_area_top[c] = z;
            }
            else if (i == 1)
            {
                x_area_rig[c] = z;
            }
            else if (i == 2)
            {
                y_area_bot[c] = z;
            }
            else if (i == 3)
            {
                x_area_lef[c] = z;
            }
        }
        ++c;
    }

    opticalMethod(image1, image2, meshY, meshX, y_area_top,x_area_rig, y_area_bot, x_area_lef, Ny1, Nx1, k);

    meshFile.close();
    roiFile.close();

    finish = clock();
    std::cout << "runtime = " << (double)(finish - start) / CLOCKS_PER_SEC << std::endl; // время работы программы  

    //delete[] meshX; delete[] meshY;
    //delete[] Nx1; delete[] Ny1;
    //delete[] y_area_top;
    //delete[] x_area_rig;
    //delete[] y_area_bot;
    //delete[] x_area_lef;

    return 0;
}

