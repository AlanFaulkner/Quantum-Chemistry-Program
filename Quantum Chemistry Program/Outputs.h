#pragma once

#include <iostream>       // STD C++ Lib
#include <fstream>		  // read and wirte to/from file
#include <iomanip>        // formating of output stream
#include <vector>		  // allows use of vectors - varible size arrys

#include "Input.h"

class Outputs
{
public:

	void Print_molecular_info();
	void Print_Geometry();

	Outputs();
	~Outputs();

private:

	void Print_Coordinates(std::vector<std::vector<double>> Coordinates, int Destination);
	void Print_Moment_of_inertia(std::vector<std::vector<double>> Coordinates, int NoAtoms, int Destination);

	double bond(std::vector<std::vector<double>>Input_Coordinates, int i, int j);
	double unit(std::vector<std::vector<double>>Input_Coordinates, int cart, int a, int b);
	double angle(std::vector<std::vector<double>>Input_Coordinates, int a, int b, int c);

	void Print_Bonds(std::vector<std::vector<double>> Coordinates, int Number_ofatoms, int Destination);
	void Print_Angles(std::vector<std::vector<double>> Coordinates, int Number_of_atoms, int Destination);
	void Print_OOP(std::vector<std::vector<double>> Coordinates, int Number_of_atoms, int Destination);
	void Print_Torsions(std::vector<std::vector<double>> Coordinates, int Number_of_atoms, int Destination);
};

