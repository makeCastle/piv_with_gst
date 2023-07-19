#pragma once
#include<fstream>
#include<vector>


std::vector<std::vector<double>> medianTest(std::vector<double> w, std::vector<double> h, std::vector<double> error,
	double Thr, double bx, double by, int Nx, int Ny) {
	double eps = 0.1; // noise threshold level

	std::ifstream field;
	field.open("C:/Users/1/OneDrive/Документы/информатика, математика и физика/Курсовая/field.txt"); // pathname of mesh
	// Ny, Nx - sizes of mesh

	std::vector<std::vector<double>> w_field;
	std::vector<std::vector<double>> h_field; // extended flow field displacements
	
	for (int i = 0; i < Ny + 2 * by; i++)
	{
		std::vector<double> row(Nx + 2 * bx, 0.0);
		w_field.push_back(row);
	}
	h_field = w_field;
	std::vector<std::vector<double>> w_median = w_field; std::vector<std::vector<double>> h_median = w_field; // filtered flow field displacements

	std::vector<std::vector<double>> error_field; // number of outliers
	for (int i = 0; i < Ny; i++)
	{
		std::vector<double> row(Nx, 0.0);
		error_field.push_back(row);
	}

	for (int j = 0; j < Ny; ++j) {
		for (int i = 0; i < Nx; ++i) {
			int k = j * Nx + i + 1; // number of mesh point
			error_field[j][i] = error[k];
			w_field[j + by][i + bx] = w[k]; // fill in field format
			h_field[j + by][i + bx] = h[k]; // fill in field format
		}
	}
	for (int g = 0; g < by; ++g) {
		for (int i = 0; i < Nx + 2*bx; i++)
		{
			w_field[g][i] = w_field[by][i]; h_field[g][i] = h_field[by][i]; // fill the top borders
			w_field[Ny + by + g][i] = w_field[by + Ny - 1][i]; h_field[Ny + by + g][i] = h_field[by + Ny - 1][i]; // fill the bottom borders
		}
	}
	for (int g = 0; g < bx; ++g) {
		for (int i = 0; i < Ny + 2 * by; i++)
		{
			w_field[i][g] = w_field[i][bx];   h_field[i][g] = h_field[i][bx]; // fill the left borders
			w_field[i][Nx + bx + g] = w_field[i][bx + Nx - 1];   h_field[i][Nx + bx + g] = h_field[i][bx + Nx - 1]; // fill the right borders
		}
	}

	// BEGIN - filtered flow field
	for (int j = by; j < Ny + by; ++j) {
		for (int i = bx; i < Nx + bx; ++i) {
			// FOR w - components
			std::vector<std::vector<double>> Neigh; // data of neighborhood with center - point
			for (int k = j - by; k < j + by; ++k)
			{
				std::vector<double> row;
				for (int p = i - bx; p < i + bx; ++p)
				{
					row.push_back(w_field[j][i]);
				}
				Neigh.push_back(row);
			}
			std::vector<std::vector<double>> NeighCol = Neigh(:); // in column format
				NeighCol2_w = [NeighCol(1:(2 * bx + 1) * by + bx); NeighCol((2 * bx + 1) * by + bx + 2:end)];
			// neighborhood exluding center - point
				Median = median(NeighCol2_w); // median of neighborhood
				Fluct = w_field(j, i) - Median; // fluctuation with respect to median
				Res = NeighCol2_w - Median; // residual: neighborhood fluctuation
				MedianRes = median(abs(Res)); // median(absolute) value of residual
				Normfluct_w = abs(Fluct / (MedianRes + eps)); // normalized fluctuation of median residual

			// FOR h - components
				Neigh = h_field(j - by:j + by, i - bx : i + bx); // data of neighborhood with center - point
				NeighCol = Neigh(:); // in column format
				NeighCol2_h = [NeighCol(1:(2 * bx + 1) * by + bx); NeighCol((2 * bx + 1) * by + bx + 2:end)];
			// neighborhood exluding center - point
				Median = median(NeighCol2_h); // median of neighborhood
				Fluct = h_field(j, i) - Median; // fluctuation with respect to median
				Res = NeighCol2_h - Median; // residual: neighborhood fluctuation
				MedianRes = median(abs(Res)); // median(absolute) value of residual
				Normfluct_h = abs(Fluct / (MedianRes + eps)); // normalized fluctuation of median residual

				if (sqrt(Normfluct_w ^ 2 + Normfluct_h ^ 2) > Thr) {
					// detection criterion
						error_field(j - by, i - bx) = error_field(j - by, i - bx) + 1;
					w_median(j, i) = NaN;% mean(NeighCol2_w); // mean filter
						h_median(j, i) = NaN;% mean(NeighCol2_h); // mean filter
				}
				else {
					w_median(j, i) = w_field(j, i); // the same
						h_median(j, i) = h_field(j, i); // the same
				}
		}
	}
	// END - filtered flow field

	for (j = 1:1 : Ny) {
		for (i = 1 : 1 : Nx) {
			k = (j - 1) * Nx + i; // number of mesh point
			error(k) = error_field(j, i);
			w(k) = w_median(j + by, i + bx); // fill in column format
			h(k) = h_median(j + by, i + bx); // fill in column format
		}
	}
}