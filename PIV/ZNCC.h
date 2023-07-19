#pragma once

#include<iostream>
#include<vector>
#include"CImg.h"



std::vector<std::vector<std::vector<double>>> GST(std::vector<std::vector<double>> a, int Ny, int Nx) {
	// filling GST and GST of squares for fragment of the image
	std::vector<std::vector<std::vector<double>>> GSTable;
	GSTable = std::vector<std::vector<std::vector<double>>>(Ny + 1, std::vector<std::vector<double>>(Nx + 1, std::vector<double>(2)));

	double result1 = 0.0;
	double result2 = 0.0;
	int y1;
	int x1;

	// BEGIN - memory allocation
	for (int i = 0; i <= Nx; i++)
	{
		GSTable[0][i][0] = 0.0;
		GSTable[0][i][1] = 0.0;
	}

	for (int i = 0; i <= Ny; i++)
	{
		GSTable[i][0][0] = 0.0;
		GSTable[i][0][1] = 0.0;
	}
	// END - memory allocation

	for (int y = 1; y <= Ny; y++)
	{
		for (int x = 1; x <= Nx; x++)
		{
			y1 = y - 1;
			x1 = x - 1;
			result1 = a[y1][x1];
			result2 = a[y1][x1] * a[y1][x1];

			GSTable[y][x][0] = GSTable[y1][x][0] + GSTable[y][x1][0] - GSTable[y1][x1][0] + result1; // GST
			GSTable[y][x][1] = GSTable[y1][x][1] + GSTable[y][x1][1] - GSTable[y1][x1][1] + result2; // GST of squares
		}
	}

	std::vector<std::vector<std::vector<double>>> Result;
	Result = std::vector<std::vector<std::vector<double>>>(Ny, std::vector<std::vector<double>>(Nx, std::vector<double>(2)));

	for (int i = 1; i <= Ny; i++)
	{
		for (int j = 1; j <= Nx; j++)
		{
			Result[i - 1][j - 1][0] = GSTable[i][j][0];
			Result[i - 1][j - 1][1] = GSTable[i][j][1];
		}
	}

	return Result;
}

std::vector<double> pivZNCC(std::vector<std::vector<double>> f1, std::vector<std::vector<double>> f2, int w0, int h0, int Nx1, int Nx2, int Ny1, int Ny2,
	double aver1, double aver2, double D1, double D2) {
	// the function defines the similarity function of images f1 and f2 and preliminary shifts (w and h)
	// f1 - fragment 1th image
	// f2 - fragment 2th image

	double e = 0.000000001; // to avoid division by 0, log(-e), log(0)

	std::vector<std::vector<double>> ff;
	for (int i = 0; i < Ny2 - Ny1 + 1; i++) {
		std::vector<double> rowf(Nx2 - Nx1 + 1, 0.0);
		ff.push_back(rowf);
	}

	int NxROI = Nx2 - Nx1 + 1; // x-size of searching area
	int NyROI = Ny2 - Ny1 + 1; // y-size of searching area
	int Nx11 = Nx1; // x-size of computational domain
	int Ny11 = Ny1; // y-size of computational domain

	// BEGIN - calculation of ZNCC
	for (int j = 0; j < NyROI; ++j) {
		for (int i = 0; i < NxROI; ++i) {
			for (int y = 0; y < Ny11; ++y) {
				for (int x = 0; x < Nx11; ++x) {			
					ff[j][i] = ff[j][i] + (f1[y][x] - aver1) * (f2[j + y][i + x] - aver2);
				}
			}
			ff[j][i] = ff[j][i] / (D1 * D2) + 2.0 * e; // to avoid division by log(-e)
		}
	}
	// END - calculation of ZNCC

	for (int i = 0; i < NxROI; i++) // correction of the CC in x at the top border
	{
		ff[0][i] = ff[1][i];
	}

	for (int i = 0; i < NxROI; i++) // correction of the CC in x at the bottom border
	{
		ff[NyROI - 1][i] = ff[NyROI - 2][i];
	}

	for (int i = 0; i < NyROI; i++) // correction of the CC in x at the left border
	{
		ff[i][0] = ff[i][1];
	}

	for (int i = 0; i < NyROI; i++) // correction of the CC in x at the right border
	{
		ff[i][NxROI - 1] = ff[i][NxROI - 2];
	}

	//std::vector<double> displace; std::vector<int> hVector; // searching of 1th minimum FS in x
	//hVector.resize(Nx2 - Nx1 - 3); displace.resize(Nx2 - Nx1 - 3);
	//for (int i = 1; i < Nx2 - Nx1 - 1; i++)
	//{
	//	displace.push_back(ff[1][i]);
	//	hVector.push_back(1);
	//}

	//for (int i = 1; i < Nx2 - Nx1 - 1; i++)
	//{
	//	for (int j = 2; j < Ny2 - Ny1 - 1; j++)
	//	{
	//		if (ff[j][i] > ff[j - 1][i])
	//		{
	//			displace[i - 1] = ff[j][i];
	//			hVector[i - 1] = j;
	//		}
	//	}
	//}

	//int w = 0;
	//for (int i = 1; i < Nx2 - Nx1 - 3; i++)
	//{
	//	if (displace[i] > displace[i - 1])
	//	{
	//		w = i;
	//	}
	//}

	//int h = hVector[w];

	//w = w + 1;
	//h = h + 1;

	double max1ff = -1.0; double max2ff = -1.0;
	double w = 0.0; double h = 0.0;

	for (int j = 1; j < Ny2 - Ny1; ++j) {
		for (int i = 1; i < Nx2 - Nx1; ++i) {
			if (((ff[j][i] >= ff[j][i + 1]) && (ff[j][i] >= ff[j][i - 1]) && (ff[j][i] >= ff[j + 1][i]) && (ff[j][i] >= ff[j - 1][i]) && (ff[j][i] > max1ff))) {
				max1ff = ff[j][i];
				w = i;
				h = j;
			}
		}

	}
	for (int j = 1; j < Ny2 - Ny1; ++j) {
		for (int i = 1; i < Nx2 - Nx1; ++i) {
			if (((ff[j][i] >= ff[j][i + 1]) && (ff[j][i] >= ff[j][i - 1]) && (ff[j][i] >= ff[j + 1][i]) && (ff[j][i] >= ff[j - 1][i]) && (ff[j][i] > max2ff))) {
				if (!((i == w) && (j == h))) {
					max2ff = ff[j][i];
				}
			}

		}

	}
	
	w = w + 1;
	h = h + 1;

	if (max2ff == -1) { max2ff = max1ff / 2.0; }
	double SNR = max1ff / max2ff;

	// BEGIN - sub-pixel interpolation of correct shifts

	double R2 = ff[h - 1][w - 1]; // FS in x in minimum (e is not subtracted to exclude the flat FS)
	double R1_w = ff[h - 1][w - 2] - e; // FS in x on the left
	double R3_w = ff[h - 1][w] - e; // FS in x on the right
	double R1_h = ff[h - 2][w - 1] - e; // FS in y on the left
	double R3_h = ff[h][w - 1] - e; // FS in y on the right
	w = w + (R1_w - R3_w) / (2.0 * R1_w - 4.0 * R2 + 2.0 * R3_w); // sub-pixel interpolation
	h = h + (R1_h - R3_h) / (2.0 * R1_h - 4.0 * R2 + 2.0 * R3_h); // sub-pixel interpolation


	double w_exact = -(w0 - w + 1); // correction of the column number with the minimum value
	double h_exact = h0 - h + 1; // correction of the row number with the minimum value
	// END - sub-pixel interpolation of correct shifts
	std::vector<double> res = { w_exact, h_exact, SNR };
	return res;
}



