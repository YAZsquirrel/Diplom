#pragma once
#include <fstream>
#include <iostream>
#include "GridMaker.h"
#include <direct.h>
#include <iomanip>
#include <string>

using namespace mesh_comps;

/*Переопределение операторов для вывода и ввода в бинарный файл*/
__forceinline std::ostream& operator < (std::ostream& file, const double& data)
{
	file.write((char*)&data, sizeof(data));
	return file;
}
__forceinline std::ostream& operator < (std::ostream& file, const int& data)
{
	file.write((char*)&data, sizeof(data));
	return file;
}

__forceinline std::istream& operator > (std::istream& file, double& data)
{
	file.read((char*)&data, sizeof(data));
	return file;
}

__forceinline std::istream& operator > (std::istream& file, int& data)
{
	file.read((char*)&data, sizeof(data));
	return file;
}

void Output2D(int it, double* q, mesh_comps::Mesh *mesh)
{
	const int SCALE_PRINT = 1; // Не знаю что это такое
	std::cout << "!!!!!!!!!!!!!Scale = " << SCALE_PRINT << std::endl;
	int kuslov = mesh->knots.size(); //  Количество узлов сетки
	int kolel = mesh->hexas.size(); // Количество конечных элементов
	int kolVerInel = 8; // Количество вершин в конечном элементе

	auto finitElems = mesh->hexas.data();
	auto coords = mesh->knots.data();
	int i, j;
	std::string pathInput = "output2D_temperature";
	if (it == 0)
	{
		/*Создаём директорию для выходных файлов, которые будут использоваться построителем сетки
		Необходима библиотека #include <direct.h>*/
		system(std::string("rmdir /s /q " + pathInput).c_str());
		_mkdir(pathInput.c_str());


		/*Эту часть можно оставить без изменений*/
		std::string path = pathInput + "/inftry.dat";
		std::ofstream ofp;
		ofp.open(path, std::ios::binary);
		ofp << "\tISLAU=	0 INDKU1=	 0 INDFPO=	1" << std::endl;
		ofp << "KUZLOV= " << kuslov << "  KPAR= " << kolel << "    KT1= 0   KTR2= 0   KTR3= 0" << std::endl;
		ofp << "KISRS1= 0 KISRS2= 0 KISRS3= 0   KBRS= 0" << std::endl;
		ofp << "\tKT7= 0   KT10= 0   KTR4= 0  KTSIM= 0" << std::endl;
		ofp << "\tKT6= 0" << std::endl;
		ofp.close();
		ofp.clear();

		path = pathInput + "/nver.dat";
		ofp.open(path, std::ios::binary);
		/*Заполняем информацию о конечных элементах*/
		for (int i = 0; i < kolel; i++)
		{
			/*Заплняем информацию об i-ом конечном элементе*/
			for (int j = 0; j < kolVerInel; j++)
				ofp < (finitElems[i]->knots_num[j] + 1);


			/*Непонятно зачем но без этого вроде как работать не будет, так что это необходиомо оставить*/
			for (int j = 0; j < 6; j++)
				ofp < 1;

		}
		ofp.close();
		ofp.clear();

		path = pathInput + "/xyz.dat";
		ofp.open(path, std::ios::binary);
		/*Запрлняем информацию о координатах сетки*/
		for (int i = 0; i < kuslov; i++)
		{
			ofp < coords[i]->x;
			ofp < coords[i]->y;
			ofp < coords[i]->z;
		}

		ofp.close();
		ofp.clear();

		path = pathInput + "/nvkat.dat";
		ofp.open(path, std::ios::binary);
		/*Непонятно зачем но без этого вроде как работать не будет, так что это необходиомо оставить*/
		for (int i = 0; i < kolel; i++)
		{
			ofp < 1;
		}
		ofp.close();
		ofp.clear();

		path = pathInput + "\\smtr";
		ofp.open(path, std::ios::binary);
		/*Непонятно зачем но без этого вроде как работать не будет, так что это необходиомо оставить*/
		for (int uz = 0; uz < kuslov; uz++)
		{
			ofp < 1;
			ofp < 1;
		}
		ofp.close();
		ofp.clear();


		path = pathInput + "\\fields.cnf";
		ofp.open(path, std::ios::binary);
		/*Непонятно зачем но без этого вроде как работать не будет, так что это необходиомо оставить*/
		ofp << "4" << std::endl;
		ofp << std::endl;
		ofp << "Temperature" << std::endl;
		ofp << "Displacement - X" << std::endl;
		ofp << "Displacement - Y" << std::endl;
		ofp << "Displacement - Z" << std::endl;
		ofp << std::endl;
		ofp << std::endl;

		std::ofstream ouf;
		path = pathInput + "\\times_main";
		/*Заполняем временные слои вроде как (пока что временной слой один)*/
		ouf.open(path);
		ouf << 1 << std::scientific << std::setprecision(15) << std::endl;
		ouf << 1.0 << std::endl;
		/*
		* На случай нескольких временных слоёв
		ouf << TimeGrid.size() << scientific << setprecision(15) << std::endl;
		for (int i = 0; i < TimeGrid.size(); i++)
			ouf << (double)TimeGrid[i] / (3600.0 * 24.0) << std::endl;
		*/
		ouf.clear();
		ouf.close();
	}

	std::string path = pathInput + "\\sx." + std::to_string(it);
	std::ofstream ofp3(path, std::ios::binary);
	/*Заполняем значения давления на it временном слое*/
	for (int i = 0; i < kuslov; i++) {
		ofp3 < q[i];
	}

	/*
	for (int i = 0; i < mesh.koluz; i++) {
		ofp3 < T[i];
	}
	*/
	ofp3.close();
	ofp3.clear();
}
