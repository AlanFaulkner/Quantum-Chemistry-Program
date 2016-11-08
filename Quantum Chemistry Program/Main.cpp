/*
#####################################################################################
##                                                                                 ##
## Quantum Chemistry Program Version 1.0                    Alan Faulkner          ##
##                                                                                 ##
## This is based in part on tutorials presented at :                               ##
## http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming                  ##
##                                                                                 ##
## The aim of this project it to generate all relevant data and compute the energy ##
## of a given molecular system using the Hartree-Fock method and structral data    ##
## presented in a text file in the form of cartisian coordinates. All caulation    ##
## preformed using the STO-3G basis set -- this is hard coded at this point.       ##
##                                                                                 ##
#####################################################################################
*/

#define _USE_MATH_DEFINES // Allows the use of M_PI for C++ - must be placed above includes or wont work!
#include <cmath>          // Need for M_PI  and exp function
#include <iostream>       // STD C++ Lib
#include <iomanip>        // formating of output stream
#include <fstream>		  // read and wirte to/from file
#include <vector>		  // allows use of vectors - varible size arrys
#include <algorithm>	  // Sort function

#include <string>         // for use of character strings
#include <istream>		  // for getline
#include <Eigen/Dense>
#include "Input.h"
#include "Outputs.h"
#include "Orbital_Overlaps.h"


using Eigen::MatrixXd;

/*
#####################################################################################
##                              ---Test Data---                                    ##
##                       ---Data is for a water Molecule---                        ##
##   This should be read in from a text file and stored correctly at later date    ##
#####################################################################################
*/

//Coordinate Info
double Cart[7][3] = { 0.00000000,  1.43233673, -0.96104039,
					  0.00000000, -1.43233673, -0.96104039,
					  0.00000000,  0.00000000,  0.24026010,
					  0.00000000,  0.00000000,  0.24026010,
					  0.00000000,  0.00000000,  0.24026010,
				 	  0.00000000,  0.00000000,  0.24026010,
					  0.00000000,  0.00000000,  0.24026010 };

//RR values
double RR[3][3] = { 0.00000000,  1.43233673, -0.96104039,
					0.00000000, -1.43233673, -0.96104039,
					0.00000000,  0.00000000,  0.24026010 };

//Orbital Coefficents for the STO-3G Basis set
double OrbCoeff[7][3] = { 3.425250914,   0.6239137298,  0.168855404,
						  3.425250914,   0.6239137298,  0.168855404,
						130.709321400,  23.80886050,    6.443608313,
						  5.033151319,   1.169596125,  0.38038896,
						  5.033151319,   1.169596125,  0.38038896,
	                      5.033151319,   1.169596125,  0.38038896,
						  5.033151319,   1.169596125,  0.38038896 };

//Primative cut-off values for orbitals in the STO-3G Basis set
double PrimCut[7][3] = {  0.15432896730,  0.5353281423,  0.4446345422,
	   					  0.15432896730,  0.5353281423,  0.4446345422,
						  0.15432896730,  0.5353281423,  0.4446345422,
						 -0.09996722919,  0.3995128261,  0.7001154689,
						  0.15591627500,  0.6076837186,  0.3919573931,
						  0.15591627500,  0.6076837186,  0.3919573931,
						  0.15591627500,  0.6076837186,  0.3919573931 };

//Cartisian Cartesian angular values of the orbitals, in the following order: H1s H2s O1s O2s O2px O2py O2pz
int cart_ang[7][3] = { 0, 0, 0,
					   0, 0, 0, 
					   0, 0, 0,
					   0, 0, 0, 
					   1, 0, 0, 
					   0, 1, 0, 
					   0, 0, 1 };

int Z[3] = { 1,1,8 };

/*
#####################################################################################
##                         --Useful General Functions--                            ##
#####################################################################################
*/

int Factorial(int n)
{
	int factorial = 1;
	for (int i = 1; i <= n; i++) {
		factorial = factorial*i;
	}
	return factorial;
}

/*
#######################################################################################
##              --Functions specific to calculating overlap matrixes--               ##
#######################################################################################
*/
 
//Normalisation Value for overlap orbitals -- apparently i wrote this sometime ago and it works
double Norm(double alpha, int orbital)
{
	double Norm;
	Norm = pow((2 * alpha / M_PI), .75)*pow((pow(8 * alpha, orbital)*Factorial(orbital))
		/ (Factorial(2 * orbital)), .5);
	return Norm;
}

/*
##############################################################################################
##                                --Overlap Function--                                      ##
## Based on Algorithm found in the Mathematica Journal vol 16 Published Febuary 16, 2012    ##
## by MINHHUY HÔ, JULIO MANUEL, HERNÁNDEZ-PÉREZ                                             ##
##  http://www.mathematica-journal.com/2012/02/evaluation-of-gaussian-molecular-integrals/  ##
##############################################################################################
*/

void Calculate_Overlap()
{
	//Matrix size is pre-determined ATM!
	//initalise the Overlap matrix to 0.
	double OverlapMatrix[7][7];
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			OverlapMatrix[i][j] = 0;
		}
	}

	//interations for Atoms 
	for (int Atom1 = 0; Atom1 < 7; Atom1++) {
		for (int Atom2 = 0; Atom2 < 7; Atom2++) {

			//
			double Total_Overlap = 0;
			int orbital_Type_Atom1 = cart_ang[Atom1][0] + cart_ang[Atom1][1] + cart_ang[Atom1][2];
			int orbital_Type_Atom2 = cart_ang[Atom2][0] + cart_ang[Atom2][1] + cart_ang[Atom2][2];

			//Orbital indicates which set of orbital coeffcients and cut-offs are used
			for (int Orbital1 = 0; Orbital1 < 3; Orbital1++) {
				for (int Orbital2 = 0; Orbital2 < 3; Orbital2++) {

					// Sxyz is storage array for X, Y & Z componets of a give orbital overlap between two atoms
					double Sxyz[3];

					//Coordinate give X,Y,Z values.
					for (int Coordinate = 0; Coordinate < 3; Coordinate++) {

						//This array is used for storing the data generated from recussion
						//the size(a) here is set to 3 by inspection -- fix it!
						const int a = 5;
						double S[a][a];

						for (int i = 0; i < 3; i++) {
							for (int j = 0; j < 3; j++) {
								S[i][j] = 0;
							}
						}

						//initial conditions

						S[0][0] = 1;
						S[0][1] = 0;
						S[1][0] = -(Cart[Atom1][Coordinate] - ((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])));

						//calculate other Gaussian by regression.

						if (a>1)
							for (int b = 2; b<a; b++) {
								S[b][0] = -(Cart[Atom1][Coordinate] - ((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])))
									*S[b - 1][0] + ((b - 1) / (2 * (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])))*S[b - 2][0];
							}

						//Transfer equation

						for (int c = 0; c<a - 1; c++) {
							for (int d = 1; d<a; d++) {
								S[c][d] = S[c + 1][d - 1] + (Cart[Atom1][Coordinate] - Cart[Atom2][Coordinate])*S[c][d - 1];
							}
						}

						//Store calculated Cartesian component in Sxyz matrix

						Sxyz[Coordinate] = S[(cart_ang[Atom1][Coordinate])][(cart_ang[Atom2][Coordinate])];

					}

					//Calculate and normalise overlap for a given set of 2 atoms -- sum of a total of 9 seperate orbital combinations
					double EAB = exp(-1 * ((OrbCoeff[Atom1][Orbital1] * OrbCoeff[Atom2][Orbital2] / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2]))*
						         ((Cart[Atom1][0] - Cart[Atom2][0])*(Cart[Atom1][0] - Cart[Atom2][0]) + (Cart[Atom1][1] - Cart[Atom2][1])*(Cart[Atom1][1] - Cart[Atom2][1])
							     + (Cart[Atom1][2] - Cart[Atom2][2])*(Cart[Atom1][2] - Cart[Atom2][2]))));

					double Overlap = EAB*(pow((M_PI / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])), 1.5))*Sxyz[0] * Sxyz[1] * Sxyz[2];

					double Normed = Norm(OrbCoeff[Atom1][Orbital1], orbital_Type_Atom1)*Norm(OrbCoeff[Atom2][Orbital2], orbital_Type_Atom2)*PrimCut[Atom1][Orbital1] * PrimCut[Atom2][Orbital2];

					double NormOverlap = Normed*Overlap;

					//Total Overlap comprises of the sumation of 9 seperate normalised intergrals using the STO-3G Basis set
					Total_Overlap += NormOverlap;

					}
				}

			//Store final overlap value in the Overlap Matrix
			OverlapMatrix[Atom1][Atom2] = Total_Overlap;

			}
		}

	//clean up matrix ie. number smaller than 1e-8 = 0
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			if (OverlapMatrix[i][j] < 1.0e-8 && OverlapMatrix[i][j]> -1e-8) 
				OverlapMatrix[i][j] = 0;
		}
	}

	//Temporary print function to display Overlap Matrix
	std::cout << "---Overlap Matrix---" << std::endl << std::endl;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			std::cout << std::setw(10) << OverlapMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	}

/*
##############################################################################################
##                                --Kinetic Overlap Function--                              ##
## Based on Algorithm found in the Mathematica Journal vol 16 Published January 31, 2013    ##
## by MINHHUY HÔ, JULIO MANUEL, HERNÁNDEZ-PÉREZ                                             ##
## http://www.mathematica-journal.com/2013/01/evaluation-of-gaussian-molecular-integrals-2/ ##
##                                                                                          ##
## This function also calculates the Overlap Matrix (see above) as a by product             ##
##############################################################################################
*/

void Calculate_Kinetic_Overlap()
{
	//Matrix for the standard Overlap Matrix
	double OverlapMatrix[7][7];
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			OverlapMatrix[i][j] = 0;
		}
	}

	//Matrix used to store kinetic overlap values
	//Matrix size is pre-determined ATM!
	//initalise the Overlap matrix to 0.
	double Kinetic_Overlap_Matrix[7][7];
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			Kinetic_Overlap_Matrix[i][j] = 0;
		}
	}

	//interations for Atoms 
	for (int Atom1 = 0; Atom1 < 7; Atom1++) {
		for (int Atom2 = 0; Atom2 <= Atom1; Atom2++) {

			//
			double Total_Overlap = 0;
			double Total_Kinetic_Overlap = 0;

			int orbital_Type_Atom1 = cart_ang[Atom1][0] + cart_ang[Atom1][1] + cart_ang[Atom1][2];
			int orbital_Type_Atom2 = cart_ang[Atom2][0] + cart_ang[Atom2][1] + cart_ang[Atom2][2];

			//Orbital indicates which set of orbital coeffcients and cut-offs are used
			for (int Orbital1 = 0; Orbital1 < 3; Orbital1++) {
				for (int Orbital2 = 0; Orbital2 < 3; Orbital2++) {

					// Sxyz is storage array for X, Y & Z componets of a give orbital overlap between two atoms
					double Sxyz[3];
					// Kxyz is storage array for X, Y & Z componets of the Kinetic Overlap between two atoms
					double Kxyz[3];

					//Coordinate give X,Y,Z values.
					for (int Coordinate = 0; Coordinate < 3; Coordinate++) {

						//This array is used for storing the data generated from recussion
						//the size(a) here is set to 3 by inspection -- fix it!
						const int a = 4;
						double S[a][a];

						for (int i = 0; i < a; i++) {
							for (int j = 0; j < a; j++) {
								S[i][j] = 0;
							}
						}

						//Stores kinetic values a >= 4 by trial and error -- fix it!
						double K[a][a];
						for (int i = 0; i < a; i++) {
							for (int j = 0; j < a; j++) {
								K[i][j] = 0;
							}
						}

						//initial conditions

						S[0][0] = 1;
						S[0][1] = 0;
						S[1][0] = -(Cart[Atom1][Coordinate] - ((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])));

						//calculate other Gaussian by regression.

						if (a > 1)
							for (int b = 2; b < a; b++) {
								S[b][0] = -(Cart[Atom1][Coordinate] - ((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])))
									*S[b - 1][0] + ((b - 1) / (2 * (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])))*S[b - 2][0];
							}

						//Transfer equation

						for (int c = 0; c < a - 1; c++) {
							for (int d = 1; d < a; d++) {
								S[c][d] = S[c + 1][d - 1] + (Cart[Atom1][Coordinate] - Cart[Atom2][Coordinate])*S[c][d - 1];
							}
						}

						//Store calculated Cartesian component in Sxyz matrix

						Sxyz[Coordinate] = S[(cart_ang[Atom1][Coordinate])][(cart_ang[Atom2][Coordinate])];


						//Inital conditions for kinetic intergral

						K[0][0] = 2 * OrbCoeff[Atom1][Orbital1] * OrbCoeff[Atom2][Orbital2] * S[1][1];

						for (int e = 1; e<a - 1; e++) {
							K[e][0] = -e*OrbCoeff[Atom2][Orbital2] * S[e - 1][1] + 2 * OrbCoeff[Atom1][Orbital1] * OrbCoeff[Atom2][Orbital2] * S[e + 1][1];
						}


						for (int f = 1; f<a - 1; f++) {
							K[0][f] = -f*OrbCoeff[Atom1][Orbital1] * S[1][f - 1] + 2 * OrbCoeff[Atom1][Orbital1] * OrbCoeff[Atom2][Orbital2] * S[1][f + 1];
						}


						//fills in rest of matrix
						for (int g = 1; g<a - 1; g++) {
							for (int h = 1; h<a - 1; h++) {
								K[g][h] = ((g*h*S[g - 1][h - 1])
									- (2 * g*OrbCoeff[Atom2][Orbital2] * S[g - 1][h + 1])
									- (2 * h*OrbCoeff[Atom1][Orbital1] * S[g + 1][h - 1])
									+ (4 * OrbCoeff[Atom1][Orbital1] * OrbCoeff[Atom2][Orbital2] * S[g + 1][h + 1])) / 2;
							}
						}


						Kxyz[Coordinate] = K[(cart_ang[Atom1][Coordinate])][(cart_ang[Atom2][Coordinate])];

					}

					//Calculate both kinetic and Overlap matrix values

					double EAB = exp(-1 * ((OrbCoeff[Atom1][Orbital1] * OrbCoeff[Atom2][Orbital2] / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2]))*
						((Cart[Atom1][0] - Cart[Atom2][0])*(Cart[Atom1][0] - Cart[Atom2][0]) + (Cart[Atom1][1] - Cart[Atom2][1])*(Cart[Atom1][1] - Cart[Atom2][1])
							+ (Cart[Atom1][2] - Cart[Atom2][2])*(Cart[Atom1][2] - Cart[Atom2][2]))));

					double Normanalisation_Const = Norm(OrbCoeff[Atom1][Orbital1], orbital_Type_Atom1)*Norm(OrbCoeff[Atom2][Orbital2], orbital_Type_Atom2)*PrimCut[Atom1][Orbital1] * PrimCut[Atom2][Orbital2];

					//Noralised Overlap matrix

					double Overlap = EAB*(pow((M_PI / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])), 1.5))*Sxyz[0] * Sxyz[1] * Sxyz[2];

					double Norm_Overlap = Normanalisation_Const*Overlap;

					//Total Overlap comprises of the sumation of 9 seperate normalised intergrals using the STO-3G Basis set
					Total_Overlap += Norm_Overlap;


					//Normalised Kinetic Matrix
					double Kinetic_Overlap = EAB*(pow((M_PI / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])), 1.5))*((Kxyz[0] * Sxyz[1] * Sxyz[2]) + (Sxyz[0] * Kxyz[1] * Sxyz[2]) + (Sxyz[0] * Sxyz[1] * Kxyz[2]));

					double Norm_Kinetic_Overlap = Normanalisation_Const*Kinetic_Overlap;
					Total_Kinetic_Overlap += Norm_Kinetic_Overlap;

				}
			}

			OverlapMatrix[Atom1][Atom2] = Total_Overlap;
			Kinetic_Overlap_Matrix[Atom1][Atom2] = Total_Kinetic_Overlap;
		}
	}

	//clean up matrix ie. number smaller than 1e-8 = 0
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			if (OverlapMatrix[i][j] < 1.0e-8 && OverlapMatrix[i][j]> -1e-8)
				OverlapMatrix[i][j] = 0;
		}
	}

	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			if (Kinetic_Overlap_Matrix[i][j] < 1.0e-8 && Kinetic_Overlap_Matrix[i][j]> -1e-8)
				Kinetic_Overlap_Matrix[i][j] = 0;
		}
	}

	//Temporary print function to display Overlap Matrix
	std::cout << "---Overlap Matrix---" << std::endl << std::endl;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			std::cout << std::setw(10) << OverlapMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;
	std::cout << "---Kinetic Overlap Matrix---" << std::endl << std::endl;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			std::cout << std::setw(10) << Kinetic_Overlap_Matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

/*
##############################################################################################
##                            --One Electron intergral Function--                           ##
## Based on Algorithm found in the Mathematica Journal vol 16 Published September 15, 2014  ##
## by MINHHUY HÔ, JULIO MANUEL, HERNÁNDEZ-PÉREZ                                             ##
## http://www.mathematica-journal.com/2014/09/evaluation-of-gaussian-molecular-integrals-3/ ##
##                                                                                          ##
## Gauss-Chebyshev quadrature numerical integration Algorithum based on fortran 77 code by  ##
## J. Pérez-Jorda and E. San-Fabián, Computer Physics Communications, 77(1), 1993 pp. 46-56 ##
##############################################################################################
*/

/*
##############################################################################################
## These values are for the limits of intergraion when using the Gauss-Chebyshev quadrature ##
## numerical integration Algorithum                                                         ##
##############################################################################################
*/
double abscissa(int n, int i)
{
	double A = (n + 1.0 - 2.0*i) / (n + 1.0) + (2.0 / M_PI)*(1.0 + (2.0 / 3.0)*pow(sin((i*M_PI) / (n + 1.0)), 2))*cos(i*M_PI / (n + 1.0))*sin(i*M_PI / (n + 1.0));
	return A;
}

double Omega(int n, int i)
{
	double A = (16.0 / (3.0*(n + 1.0)))*pow(sin(i*M_PI / (n + 1.0)), 4);
	return A;
}

/*
##############################################################################################
## Functions that are called to solve OEI.                                                  ##
## Fx= determines the function to be intergrated                                            ##
## CHP = Gauss-Chebyshev quadrature numerical integration Algorithum                        ##
##############################################################################################
*/

double Fx(double X, int Atom1, int Atom2, int rr, int Orbital1, int Orbital2)
{
	// Sxyz is storage array for X, Y & Z componets of a given orbital overlap between two atoms
	double Sxyz[3];

	//Coordinate give X,Y,Z values.
	for (int Coordinate = 0; Coordinate < 3; Coordinate++) {

		//This array is used for storing the data generated from recussion
		//the size(a) here is set to 3 by inspection -- fix it!
		const int a = 5;
		double S[a][a];

		for (int i = 0; i < a; i++) {
			for (int j = 0; j < a; j++) {
				S[i][j] = 0;
			}
		}

		//initial conditions

		S[0][0] = 1;
		S[0][1] = 0;
		S[1][0] = -1.0*(Cart[Atom1][Coordinate] - ((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2]))
			+ pow(((X + 1) / 2.0), 2) * (((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])) - RR[rr][Coordinate]));
	
		//Recurance index

		if (a>1)
			for (int b = 2; b<a; b++) {
				S[b][0] = -1.0*(Cart[Atom1][Coordinate] - ((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2]))
					+ pow(((X + 1) / 2.0), 2) * (((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2])) - RR[rr][Coordinate]))
					*S[b - 1][0] + ((b - 1) / (2 * (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2]))) * (1-pow(((X + 1) / 2.0), 2)) *S[b - 2][0];
			}

		//Transfere Equation

		for (int c = 0; c<a - 1; c++) {
			for (int d = 1; d<a; d++) {
				S[c][d] = S[c + 1][d - 1] + (Cart[Atom1][Coordinate] - Cart[Atom2][Coordinate]) * S[c][d - 1];
			}
		}

		//Store calculated Cartesian component in Sxyz matrix

		Sxyz[Coordinate] = S[(cart_ang[Atom1][Coordinate])][(cart_ang[Atom2][Coordinate])];
	}

	//std::cout << "Sx = " << Sxyz[0] << " Sy = " << Sxyz[1] << " Sz = " << Sxyz[2] << std::endl;

	double Dotproduct = 0;

	for (int Coordinate = 0; Coordinate < 3; Coordinate++) {

		//double Coord_Value = ((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2]) - Cart[Atom1][Coordinate]);
		double Coord_Value = ((OrbCoeff[Atom1][Orbital1] * Cart[Atom1][Coordinate] + OrbCoeff[Atom2][Orbital2] * Cart[Atom2][Coordinate]) / (OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2]) - RR[rr][Coordinate]);
		//std::cout << "coord = " << Coord_Value << std::endl;
		Coord_Value *= Coord_Value;
		Dotproduct += Coord_Value;
	}

	//std::cout << "Dotproduct = " << Dotproduct << std::endl;

	//What we are intergrating
	double Fx = 0.5 * (exp(-1.0*((OrbCoeff[Atom1][Orbital1] + OrbCoeff[Atom2][Orbital2]) * (pow(((X + 1) / 2.0), 2)) * Dotproduct)) *Sxyz[0] * Sxyz[1] * Sxyz[2]);
	
	return Fx;
}

double CHP(double eps, double M, int Atom1, int Atom2, int RR, int Alpha, int Beta)

//Numerical integration routine translated from Fortran 44 code printed in
//Integrates between 0 and 1.

{
	double S, C, xp;

	int n = 3;
	double err = 10.0;
	double C0 = cos(M_PI / 6.0);
	double S0 = sin(M_PI / 6.0);
	double C1 = S0;
	double S1 = C0;
	double q = (Fx(abscissa(2.0, 1.0), Atom1, Atom2, RR, Alpha, Beta) + Fx(-abscissa(2.0, 1.0), Atom1, Atom2, RR, Alpha, Beta))*Omega(2.0, 1.0);
	double p = Fx(0.0, Atom1, Atom2, RR, Alpha, Beta);
	double chp = q + p;

	int j = 0;

	while (err > eps && (2.0*n*(1 - j) + j * 4 * n / 3.0 - 1.0) < M) {

		j = 1 - j;
		C1 = j*C1 + (1 - j)*C0;
		S1 = j*S1 + (1 - j)*S0;
		C0 = j*C0 + (1 - j)*sqrt((1 + C0)*0.5);
		S0 = j*S0 + (1 - j)*S0 / (C0 + C0);
		C = C0;
		S = S0;

		for (double i = 1; i<n - 1; i = i + 2) {
			xp = 1 + 2.0 / (3.0*M_PI)*S*C*(3.0 + 2.0*S*S) - i / n;
			if (ceil(3 * (i + j + j) / 3.0)>i + j) {
				chp = chp + (Fx(-xp, Atom1, Atom2, RR, Alpha, Beta) + Fx(xp, Atom1, Atom2, RR, Alpha, Beta))*pow(S, 4);
			}
			xp = S;
			S = S*C1 + C*S1;
			C = C*C1 - xp*S1;

		}

		n = (1 + j)*n;

		p = p + (1 - j)*(chp - q);
		err = 16 * abs((1 - j)*(q - 3 * p / 2.0) + j*(chp - 2 * q)) / (3.0*n);
		q = (1 - j)*q + j*chp;

	}
	
	chp = 16 * q / (3 * n);

	//std::cout << "CHP = " << chp << std::endl;

	//Calculate and normalise overlap for a given set of 2 atoms -- sum of a total of 9 seperate orbital combinations
	double EAB = exp(-1 * ((OrbCoeff[Atom1][Alpha] * OrbCoeff[Atom2][Beta] / (OrbCoeff[Atom1][Alpha] + OrbCoeff[Atom2][Beta]))*
		((Cart[Atom1][0] - Cart[Atom2][0])*(Cart[Atom1][0] - Cart[Atom2][0]) + (Cart[Atom1][1] - Cart[Atom2][1])*(Cart[Atom1][1] - Cart[Atom2][1])
			+ (Cart[Atom1][2] - Cart[Atom2][2])*(Cart[Atom1][2] - Cart[Atom2][2]))));
	
	//Calculates resultant of intergration -- Un-normalised!!
	double result = EAB * ((2 * M_PI) / (OrbCoeff[Atom1][Alpha] + OrbCoeff[Atom2][Beta])) * chp;

	return result;
}

/*
#############################################################################################
## Driver function to pass in relavate combinations of atoms values to calculate OEI       ##
#############################################################################################
*/

void Calculate_TwoElectron_Overlap()
{
	
	std::vector<std::vector<double>> OverlapMatrix;

	//interations for Atoms 
	for (int Atom1 = 0; Atom1 < 2; Atom1++) {
		for (int Atom2 = 0; Atom2 <= Atom1; Atom2++) {
			for (int Atom3 = 0; Atom3 < 2; Atom3++) {
				for (int Atom4 = 0; Atom4 <=Atom3; Atom4++) {
					std::vector<double>  row;
					row.push_back(Atom1);
					row.push_back(Atom2);
					row.push_back(Atom3);
					row.push_back(Atom4);
					double RR_Total_OEI[3] = { 0,0,0 };
					double new_Total_OEI = 0;
					int orbital_Type_Atom1 = cart_ang[Atom1][0] + cart_ang[Atom1][1] + cart_ang[Atom1][2];
					int orbital_Type_Atom2 = cart_ang[Atom2][0] + cart_ang[Atom2][1] + cart_ang[Atom2][2];
					int orbital_Type_Atom3 = cart_ang[Atom3][0] + cart_ang[Atom3][1] + cart_ang[Atom3][2];
					int orbital_Type_Atom4 = cart_ang[Atom4][0] + cart_ang[Atom4][1] + cart_ang[Atom4][2];

					//RR running from 1-3 
					for (int RR = 0; RR < 3; RR++) {

						double Total_OEI = 0;

						//Orbital indicates which set of orbital coeffcients and cut-offs are used
						for (int Orbital1 = 0; Orbital1 < 3; Orbital1++) {

							for (int Orbital2 = 0; Orbital2 < 3; Orbital2++) {


								double Normed = Norm(OrbCoeff[Atom1][Orbital1], orbital_Type_Atom1)*Norm(OrbCoeff[Atom2][Orbital2], orbital_Type_Atom2)*PrimCut[Atom1][Orbital1] * PrimCut[Atom2][Orbital2];
								double Normed2 = Norm(OrbCoeff[Atom3][Orbital1], orbital_Type_Atom3)*Norm(OrbCoeff[Atom4][Orbital2], orbital_Type_Atom4)*PrimCut[Atom3][Orbital1] * PrimCut[Atom4][Orbital2];
								double OEI = Normed * CHP(10e-10, 50000, Atom1, Atom2, RR, Orbital1, Orbital2);
								double TEI = Normed2 * CHP(10e-10, 50000, Atom3, Atom4, RR, Orbital1, Orbital2);

								Total_OEI += (OEI + TEI);
							}
						}
						RR_Total_OEI[RR] = Total_OEI;//*-1.0;// *Z[RR];
					}
					row.push_back(RR_Total_OEI[0] + RR_Total_OEI[1] + RR_Total_OEI[2]);
					OverlapMatrix.push_back(row);
				}
			}
		}
		
	}

	//Temporary print function to display Overlap Matrix
	std::cout << "---Two Electron Intergrals---" << std::endl << std::endl;
	for (int i = 0; i < OverlapMatrix.size(); i++) {
		for (int j = 0; j < 5; j++) {
			std::cout << std::setw(10) << OverlapMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

}

/*
############################################################################################
##  Functions related to molecular geometry calculations                                  ##
############################################################################################
*/

std::vector<std::vector<double>> Bond_Lengths(std::vector<std::vector<double>> Input_Coordinates, int Number_of_atoms)
{
	int Tri = ((Number_of_atoms*(Number_of_atoms + 1)) / 2) - Number_of_atoms;
	std::vector<std::vector<double>> Bond_lengths(Tri, std::vector<double>(3));
	
	int line = 0;

	for (int i = 0; i < Number_of_atoms; i++) {
		for (int j = 0; j < i; j++) {
			Bond_lengths[line][0] = i;
			Bond_lengths[line][1] = j;
			Bond_lengths[line][2] = sqrt(((Input_Coordinates[i][1] - Input_Coordinates[j][1]) * (Input_Coordinates[i][1] - Input_Coordinates[j][1]))
				+ ((Input_Coordinates[i][2] - Input_Coordinates[j][2]) * (Input_Coordinates[i][2] - Input_Coordinates[j][2]))
				+ ((Input_Coordinates[i][3] - Input_Coordinates[j][3]) * (Input_Coordinates[i][3] - Input_Coordinates[j][3])));
			
			line++;

		}
	}

	Bond_lengths.shrink_to_fit();

	return Bond_lengths;
}

double bond(std::vector<std::vector<double>>Input_Coordinates, int i, int j)
{
	return sqrt(((Input_Coordinates[i][1] - Input_Coordinates[j][1])*(Input_Coordinates[i][1] - Input_Coordinates[j][1]))
		+ ((Input_Coordinates[i][2] - Input_Coordinates[j][2])*(Input_Coordinates[i][2] - Input_Coordinates[j][2]))
		+ ((Input_Coordinates[i][3] - Input_Coordinates[j][3])*(Input_Coordinates[i][3] - Input_Coordinates[j][3])));
}

double unit(std::vector<std::vector<double>>Input_Coordinates, int cart, int a, int b)
{
	return -(Input_Coordinates[a][cart] - Input_Coordinates[b][cart]) / bond(Input_Coordinates, a, b);
}

double angle(std::vector<std::vector<double>>Input_Coordinates, int a, int b, int c)
{
	return acos(unit(Input_Coordinates, 1, b, a) * unit(Input_Coordinates, 1, b, c)
		+ unit(Input_Coordinates, 2, b, a) * unit(Input_Coordinates, 2, b, c)
		+ unit(Input_Coordinates, 3, b, a) * unit(Input_Coordinates, 3, b, c));
}

std::vector<std::vector<double>> Bond_Angles(std::vector<std::vector<double>>Input_Coordinates, int Number_of_atoms)
{
	std::vector<std::vector<double>> Angle(Number_of_atoms, std::vector<double>(4));

	int line = 0;
	for (int i = 0; i < Number_of_atoms; i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < j; k++) {
				if (bond(Input_Coordinates, i, j) < 4 && bond(Input_Coordinates, j, k) < 4) {
					Angle[line][0] = i;
					Angle[line][1] = j;
					Angle[line][2] = k;
					Angle[line][3] = (acos(unit(Input_Coordinates, 1, j, i) * unit(Input_Coordinates, 1, j, k)
						+ unit(Input_Coordinates, 2, j, i) * unit(Input_Coordinates, 2, j, k)
						+ unit(Input_Coordinates, 3, j, i) * unit(Input_Coordinates, 3, j, k)))*(180.0 / acos(-1.0));

					line++;
				}
			}
		}
	}

	return Angle;
}

std::vector<std::vector<double>> Out_Of_Plane_Angles(std::vector<std::vector<double>> Input_Data, int Number_of_atoms)
{
	std::vector<std::vector<double>> OOP(Number_of_atoms*Number_of_atoms, std::vector<double>(5));

	int line = 0;
	for (int i = 0; i < Number_of_atoms; i++) {
		for (int k = 0; k < Number_of_atoms; k++) {
			for (int j = 0; j < Number_of_atoms; j++) {
				for (int l = 0; l < j; l++) {
					if (i != j && i != k && i != l && j != k && k != l && bond(Input_Data, i, k) < 4 && bond(Input_Data, k, j) < 4 && bond(Input_Data, k, l) < 4) {
						OOP[line][0] = i;
						OOP[line][1] = j;
						OOP[line][2] = k;
						OOP[line][3] = l;

						double ox = (unit(Input_Data, 2, k, j) * unit(Input_Data, 3, k, l) - unit(Input_Data, 3, k, j) * unit(Input_Data, 2, k, l));
						double oy = (unit(Input_Data, 3, k, j) * unit(Input_Data, 1, k, l) - unit(Input_Data, 1, k, j) * unit(Input_Data, 3, k, l));
						double oz = (unit(Input_Data, 1, k, j) * unit(Input_Data, 2, k, l) - unit(Input_Data, 2, k, j) * unit(Input_Data, 1, k, l));

						double oxx = ox*unit(Input_Data, 1, k, i);
						double oyy = oy*unit(Input_Data, 2, k, i);
						double ozz = oz*unit(Input_Data, 3, k, i);

						double theta = (oxx + oyy + ozz) / sin(angle(Input_Data, j, k, l));

						if (theta < -1.0) {
							theta = asin(-1.0);
						}
						else if (theta > 1.0){
							theta = asin(1.0);
						}
						else {
							theta = asin(theta);
						}

						OOP[line][4] = theta*(180.0 / acos(-1.0));

						line++;

					}

				}

			}
		}
	}

	return OOP;
}

int main()
{
	std::cout <<
			"##############################################################################" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"## Quantum Chemistry Program Version 1.0                    Alan Faulkner   ##" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"## This is based in part on tutorials presented at :                        ##" << std::endl <<
			"## http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming           ##" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"## Code for formation of Overlap Matrixes is base on:                       ##" << std::endl <<
			"## Mathematica Journal vol 16 Feb 16, 2012                                  ##" << std::endl <<
			"## Mathematica Journal vol 16 Jan 31, 2013                                  ##" << std::endl <<
			"## Mathematica Journal vol 16 Sep 15, 2014                                  ##" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"## The aim of this project is to take a formatted input file that contains  ##" << std::endl <<
			"## settings and coordinate information for a given small organic molecule   ##" << std::endl <<
			"## Store this information in a series of arrays. use this information to    ##" << std::endl <<
			"## calculate basic geometical pramaters such as bond length, angles etc.    ##" << std::endl <<
			"## use the provided coordinates to build the overlap matrix, the kinetic    ##" << std::endl <<
			"## overlap matirix, nuclear attraction intergrals and two electon           ##" << std::endl <<
			"## intergrals using the STO-3G basis set information provided in an         ##" << std::endl <<
			"## external file. Finally use these results to calculate the single point   ##" << std::endl <<
			"## energy of the molecule using the Hartree_Fock Method                     ##" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"##############################################################################" << std::endl << std::endl;

		std::ofstream Output("Output.txt");

		Output <<
			"##############################################################################" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"## Quantum Chemistry Program Version 1.0                    Alan Faulkner   ##" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"## This is based in part on tutorials presented at :                        ##" << std::endl <<
			"## http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming           ##" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"## Code for formation of Overlap Matrixes is base on:                       ##" << std::endl <<
			"## Mathematica Journal vol 16 Feb 16, 2012                                  ##" << std::endl <<
			"## Mathematica Journal vol 16 Jan 31, 2013                                  ##" << std::endl <<
			"## Mathematica Journal vol 16 Sep 15, 2014                                  ##" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"## The aim of this project is to take a formatted input file that contains  ##" << std::endl <<
			"## settings and coordinate information for a given small organic molecule   ##" << std::endl <<
			"## Store this information in a series of arrays. use this information to    ##" << std::endl <<
			"## calculate basic geometical pramaters such as bond length, angles etc.    ##" << std::endl <<
			"## use the provided coordinates to build the overlap matrix, the kinetic    ##" << std::endl <<
			"## overlap matirix, nuclear attraction intergrals and two electon           ##" << std::endl <<
			"## intergrals using the STO-3G basis set information provided in an         ##" << std::endl <<
			"## external file. Finally use these results to calculate the single point   ##" << std::endl <<
			"## energy of the molecule using the Hartree_Fock Method                     ##" << std::endl <<
			"##                                                                          ##" << std::endl <<
			"##############################################################################" << std::endl << std::endl;
		Output.close();

		Input Input_data;

		int Number_of_atoms = Input_data.Get_number_of_Atoms("Input.dat");
		std::vector<std::vector<double>> Coordinates = Input_data.Get_Coordinates("Input.dat", Number_of_atoms);
		int Print_Level = Input_data.Get_print_level("Input.dat");

		/*Outputs Out;
		Out.Print_molecular_info();
		Out.Print_Geometry();
		std::cout << std::defaultfloat << std::endl;*/
	Orbital_Overlaps s;
	s.build_Basis_Set_Data(Coordinates, Number_of_atoms);
	//std::vector<std::vector<double>> Overlap = s.Overlap_Matrix();
	//std::vector<std::vector<double>> Kinetic = s.Calculate_Kinetic_Overlap();
 	//std::vector<std::vector<double>> OEI = s.Calculate_OneElectron_Overlap();
	Calculate_TwoElectron_Overlap();
	return 0;
}