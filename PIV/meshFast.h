#pragma once
#include<iostream>
#include<vector>

std::vector<int> meshFast(int count, std::vector<double> NxorY, std::vector<double> meshXorY) {
	std::vector<int> meshXorYNew;  meshXorYNew.reserve(count);
	for (int p = 0; p < count; p++) {
		if ((NxorY[p] / 2.0 - floor(NxorY[p] / 2.0)) == 0.0) {
			meshXorYNew.push_back(floor(meshXorY[p]) + 0.5);
		}
		else {
			meshXorYNew.push_back(round(meshXorY[p]) + 0.5);
		}
	}
	return meshXorYNew;
}
