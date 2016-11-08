#pragma once

#include <iostream>       // STD C++ Lib
#include <fstream>		  // read and wirte to/from file
#include <vector>		  // allows use of vectors - varible size arrys
#include <algorithm>	  // Sort function
#include <string>         // for use of character strings
#include <istream>		  // for getline

#include <Eigen/Dense>    // needed to diagnalise inertia tensors
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

class Input
{
public:

	int Get_number_of_Atoms(std::string Filename);
	int Get_print_level(std::string Filename);
	std::vector<std::vector<double>> Get_Coordinates(std::string Filename, int Number_of_Atoms);
	std::vector<double> Calculate_centre_of_mass(std::vector<std::vector<double>> Coordinates, int Number_of_atoms);
	int Return_Line_Number(std::string Filename, std::string Search_Phrase);
	
	std::vector<std::vector<double>> Get_Moments_of_inertia(std::vector<std::vector<double>> Coordinates, int Number_of_atoms);

	Input();
	~Input();

private:
	
	double Atomic_Masses[11] = { 0.0, 1.0079, 4.0026, 6.941, 9.0122, 10.811, 12.0107, 14.0067, 15.9994, 18.9984, 20.1797 };
	
	

};
