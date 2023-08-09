#pragma once
#include<iostream>
#include<vector>

int* meshFast(int count, double* NxorY, double* meshXorY) {
	int* meshXorYNew{ new int[count] {} };
	for (int p = 0; p < count; p++) {
		if ((NxorY[p] / 2.0 - floor(NxorY[p] / 2.0)) == 0.0) {
			meshXorYNew[p] = floor(meshXorY[p]) + 0.5;
		}
		else {
			meshXorYNew[p] = round(meshXorY[p]) + 0.5;
		}
	}
	return meshXorYNew;
}
