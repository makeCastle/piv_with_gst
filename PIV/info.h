#pragma once
#include<iostream>
 #include<fstream>

function[it, u, v, error] = info(path, it2, count, N_frame, Nor)
// ФУНКЦИЯ ВЫПОЛНЯЕТ ПОИСК МИНИМУМОВ ФП ПОСЛЕ ФИЛЬТРАЦИИ ПРЕДВАРИТЕЛЬНО ОПРЕДЕЛЕННЫХ СМЕЩЕНИЙ
// fx - ФП по x
// fy - ФП по y
// w_filt - фильтрованные значения смещений по x
// h_filt - фильтрованные значения смещений по y

u = zeros(count, Nor, N_frame); // calculated displacements
v = zeros(count, Nor, N_frame); // calculated displacements
error = zeros(count, it2);

s = [path '\velocity\# info.txt'];
if (exist(s, 'file')) {
	it = dlmread(s, '\t'); // pathname of mesh
	for (p = 1 : 1 : count) {
		u(p, :, : ) = (dlmread([path '\velocity\' num2str(it1) 'it\u_'  num2str(p) '.txt'],'\t'))';
		v(p, :, : ) = (dlmread([path '\velocity\' num2str(it1) 'it\v_'  num2str(p) '.txt'],'\t'))';

		err = dlmread([path '\velocity\error_' num2str(p) '.txt'], '\n');
		error(p, 1:it1) = err(1:it1);
	}
	it = it + 1;
}
else {
	it = 1;
}

end
