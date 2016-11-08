#pragma once

#define _USE_MATH_DEFINES // Allows the use of M_PI for C++ - must be placed above includes or wont work!
#include <cmath>          // Need for M_PI  and exp function
#include <iostream>       // STD C++ Lib
#include <iomanip>        // formating of output stream
#include <fstream>		  // read and wirte to/from file
#include <vector>		  // allows use of vectors - varible size arrys
#include <algorithm>	  // Sort function

#include <string>         // for use of character strings
#include <sstream>

class Orbital_Overlaps
{

public:
	
	void build_Basis_Set_Data(std::vector<std::vector<double>>, int Number_of_Atoms);
	std::vector<std::vector<double>> Overlap_Matrix();
	std::vector<std::vector<double>> Calculate_Kinetic_Overlap();
	std::vector<std::vector<double>> Calculate_OneElectron_Overlap();


private:

	//Things used to build matrixes.
	std::vector<int> Orbital_List;
	int Total_Number_Of_Orbitals;
	std::string Elements[10] = { "None", "Hydrogen", "Helium", "Lithium", "Berylium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine" };
	int Return_Line_Number(std::string Filename, std::string Search_Phrase);
	
	//Matrixes required for calculating overlaps.
	std::vector<std::vector<double>> Extended_Coordinates;
	std::vector<std::vector<double>> RR;
	std::vector<std::vector<double>> Orbital_Coefficents;
	std::vector<std::vector<double>> Prime_Cutoffs;
	std::vector<std::vector<int>> Orbital_Descriptions;
	std::vector<double> Z;

	//functions used in calculating overlap matrixes
	int Factorial(int n);
	double Norm(double alpha, int orbital);

	//for OEO
	double abscissa(int n, int i);
	double Omega(int n, int i);
	double Fx(double X, int Atom1, int Atom2, int rr, int Orbital1, int Orbital2);
	double CHP(double eps, double M, int Atom1, int Atom2, int RR, int Alpha, int Beta);
};