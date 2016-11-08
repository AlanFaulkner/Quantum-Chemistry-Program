#include "Orbital_Overlaps.h"

int Orbital_Overlaps::Return_Line_Number(std::string Filename, std::string Search_Phrase)
{
	std::ifstream File(Filename);
	std::string Search = Search_Phrase;
	std::string Line;
	int Line_Number = 0;

	for (unsigned int curLine = 0; getline(File, Line); curLine++) {
		if (Line.find(Search) != std::string::npos) {
			Line_Number = curLine;
		}
	}

	File.close();

	return Line_Number;
}

//inital set up of required coordinate matrixes and orbital info
void Orbital_Overlaps::build_Basis_Set_Data(std::vector<std::vector<double>> Input_Data, int Number_of_Atoms)
{
	/*
	####################################################
	##   This Builds the RR Matrix                    ##
	##   i.e clone input coordinates - first column   ##
	####################################################
	*/

	for (int i = 0; i < Number_of_Atoms; i++) {
		std::vector<double> row;
		for (int j = 1; j < 4; j++) {
			row.push_back(Input_Data[i][j]);
		}
		RR.push_back(row);
	}
	
	/*
	####################################################
	##   This Builds the z Array                      ##
	##   i.e Nuclear charge on atom (atom no.)        ##
	####################################################
	*/
	for (int i = 0; i < Number_of_Atoms; i++)
		Orbital_Overlaps::Z.push_back(Input_Data[i][0]);

	/*
	####################################################
	##   This Builds the Orbital List Array           ##
	##   i.e The type of orbital S, P etc             ##
	####################################################
	*/

	for (int i = 0; i < Number_of_Atoms; i++) {
		if (Input_Data[i][0] <= 2)
			Orbital_List.push_back(1);

		else if (Input_Data[i][0] <= 4 && Input_Data[i][0] > 2) {
			Orbital_List.push_back(1);
			Orbital_List.push_back(1);
		}

		else if (Input_Data[i][0] <= 10 && Input_Data[i][0] > 4) {
			Orbital_List.push_back(1);
			Orbital_List.push_back(1);
			Orbital_List.push_back(3);
		}

	}

	//Calculate the total number of orbitals.
	Total_Number_Of_Orbitals = 0;
	for (int i = 0; i < Orbital_List.size(); i++)
		Total_Number_Of_Orbitals += Orbital_List[i];

	/*
	####################################################
	##   This Builds the Orbital Description Matrix   ##
	####################################################
	*/

	for (int i = 0; i < Orbital_Overlaps::Orbital_List.size(); i++) {

	if (Orbital_List[i]==1){
		std::vector<int> row;
			row.push_back(0);
			row.push_back(0);
			row.push_back(0);
			Orbital_Overlaps::Orbital_Descriptions.push_back(row);
		}

	else if (Orbital_List[i]==3){
		std::vector<int> row;
			row.push_back(1);
			row.push_back(0);
			row.push_back(0);
			Orbital_Overlaps::Orbital_Descriptions.push_back(row);
			row[0] = 0;
			row[1] = 1;
			row[2] = 0;
			Orbital_Overlaps::Orbital_Descriptions.push_back(row);
			row[0] = 0;
			row[1] = 0;
			row[2] = 1;
			Orbital_Overlaps::Orbital_Descriptions.push_back(row);
		}
	}

	/*
	##################################################################
	## Builds extended Coordinats e.g possition the of each orbital ##
	##################################################################
	*/

	for (int i = 0; i < Number_of_Atoms; i++) {
		if (Input_Data[i][0] <= 4) {
			std::vector<double> row;
			row.push_back(Input_Data[i][1]);
			row.push_back(Input_Data[i][2]);
			row.push_back(Input_Data[i][3]);
			Orbital_Overlaps::Extended_Coordinates.push_back(row);
		}

		else if (Input_Data[i][0] <= 10 && Input_Data[i][0] > 4) {
			for (int k = 0; k < 5; k++) {
				std::vector<double> row;
				row.push_back(Input_Data[i][1]);
				row.push_back(Input_Data[i][2]);
				row.push_back(Input_Data[i][3]);
				Orbital_Overlaps::Extended_Coordinates.push_back(row);
			}
		}
	}

	/*
	###########################################################
	##  Build orbital extinction and cut-off matrixes        ##
	###########################################################
	*/

	for (int i = 0; i < Number_of_Atoms; i++) {
		int orbit_count = 0;
		int copy = 1;
		int Atom_name = Z[i];
		int l = Return_Line_Number("Basis.dat", Elements[Atom_name]);

		std::ifstream Input("Basis.dat");
		std::string line;

		for (int i = 0; i <= l; i++)
			getline(Input, line);

		std::vector<double> Orb_coeff_row(3);
		std::vector<double> ext_row(3);

		while (std::getline(Input, line)) {
			if (line == "") { break; }
			for (int i = 0; i < 3; i++) {
				std::istringstream ss(line);
				ss >> Orb_coeff_row[i] >> ext_row[i];
				std::getline(Input, line);
				if (line == "") { break; }
			}
			orbit_count++;
			if (orbit_count == 2) { copy = 4; }
			for (int j = 0; j < copy; j++) {
				Orbital_Coefficents.push_back(Orb_coeff_row);
				Prime_Cutoffs.push_back(ext_row);
			}
		}
		Input.close();
	}

	//output checks

	std::cout << "-RR Matix--" << std::endl;
	for (int i = 0; i < RR.size(); i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << RR[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << "--Z matrix--" << std::endl;
	for (int i = 0; i < Z.size(); i++) {
		std::cout << Z[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "--Orbital lists--" << std::endl;
	for (int i = 0; i < Orbital_List.size(); i++) {
		std::cout << Orbital_List[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "--Orbital descriptions from class--" << std::endl;
	for (int i = 0; i < Orbital_Descriptions.size(); i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << Orbital_Descriptions[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << "--extended coordinates from class--" << std::endl;
	for (int i = 0; i < Extended_Coordinates.size(); i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << Extended_Coordinates[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << "--Orbital Coefficents from class--" << std::endl;
	for (int i = 0; i < Orbital_Coefficents.size(); i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << Orbital_Coefficents[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << "--Prime Cut-offs from class--" << std::endl;
	for (int i = 0; i < Prime_Cutoffs.size(); i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << Prime_Cutoffs[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


/*
#####################################################################################
##                         --Useful General Functions--                            ##
#####################################################################################
*/

int Orbital_Overlaps::Factorial(int n)
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
double Orbital_Overlaps::Norm(double alpha, int orbital)
{
	double Norm;
	Norm = pow((2 * alpha / M_PI), .75)*pow((pow(8 * alpha, orbital)*Orbital_Overlaps::Factorial(orbital))
		/ (Orbital_Overlaps::Factorial(2 * orbital)), .5);
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

std::vector<std::vector<double>> Orbital_Overlaps::Overlap_Matrix()
{
	std::vector<std::vector<double>> OverlapMatrix;
	
	//interations for Atoms 
			for (int Atom1 = 0; Atom1 < Total_Number_Of_Orbitals; Atom1++) {
				std::vector<double> row;
				for (int Atom2 = 0; Atom2 < Total_Number_Of_Orbitals; Atom2++) {

					//
					double Total_Overlap = 0;
					int orbital_Type_Atom1 = Orbital_Descriptions[Atom1][0] + Orbital_Descriptions[Atom1][1] + Orbital_Descriptions[Atom1][2];
					int orbital_Type_Atom2 = Orbital_Descriptions[Atom2][0] + Orbital_Descriptions[Atom2][1] + Orbital_Descriptions[Atom2][2];

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
								S[1][0] = -(Extended_Coordinates[Atom1][Coordinate] - ((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])));

								//calculate other Gaussian by regression.

								if (a > 1)
									for (int b = 2; b < a; b++) {
										S[b][0] = -(Extended_Coordinates[Atom1][Coordinate] - ((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])))
											*S[b - 1][0] + ((b - 1) / (2 * (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])))*S[b - 2][0];
									}

								//Transfer equation

								for (int c = 0; c < a - 1; c++) {
									for (int d = 1; d < a; d++) {
										S[c][d] = S[c + 1][d - 1] + (Extended_Coordinates[Atom1][Coordinate] - Extended_Coordinates[Atom2][Coordinate])*S[c][d - 1];
									}
								}

								//Store calculated Cartesian component in Sxyz matrix

								Sxyz[Coordinate] = S[(Orbital_Descriptions[Atom1][Coordinate])][(Orbital_Descriptions[Atom2][Coordinate])];

							}

							//Calculate and normalise overlap for a given set of 2 atoms -- sum of a total of 9 seperate orbital combinations
							double EAB = exp(-1 * ((Orbital_Coefficents[Atom1][Orbital1] * Orbital_Coefficents[Atom2][Orbital2] / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2]))*
								((Extended_Coordinates[Atom1][0] - Extended_Coordinates[Atom2][0])*(Extended_Coordinates[Atom1][0] - Extended_Coordinates[Atom2][0]) + (Extended_Coordinates[Atom1][1] - Extended_Coordinates[Atom2][1])*(Extended_Coordinates[Atom1][1] - Extended_Coordinates[Atom2][1])
									+ (Extended_Coordinates[Atom1][2] - Extended_Coordinates[Atom2][2])*(Extended_Coordinates[Atom1][2] - Extended_Coordinates[Atom2][2]))));

							double Overlap = EAB*(pow((M_PI / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])), 1.5))*Sxyz[0] * Sxyz[1] * Sxyz[2];

							double Normed = Norm(Orbital_Coefficents[Atom1][Orbital1], orbital_Type_Atom1)*Norm(Orbital_Coefficents[Atom2][Orbital2], orbital_Type_Atom2)*Prime_Cutoffs[Atom1][Orbital1] * Prime_Cutoffs[Atom2][Orbital2];

							double NormOverlap = Normed*Overlap;

							//Total Overlap comprises of the sumation of 9 seperate normalised intergrals using the STO-3G Basis set
							Total_Overlap += NormOverlap;

						}
					}

					//Store final overlap value in the Overlap Matrix
					row.push_back(Total_Overlap);

				}

		OverlapMatrix.push_back(row);
	}

	//clean up matrix ie. number smaller than 1e-8 = 0
	for (int i = 0; i < OverlapMatrix.size(); i++) {
		for (int j = 0; j < OverlapMatrix.size(); j++) {
			if (OverlapMatrix[i][j] < 1.0e-8 && OverlapMatrix[i][j]> -1e-8)
				OverlapMatrix[i][j] = 0;
		}
	}

	//Temporary print function to display Overlap Matrix
	std::cout << "---Overlap Matrix---" << std::endl << std::endl;
	for (int i = 0; i < OverlapMatrix.size(); i++) {
		for (int j = 0; j < OverlapMatrix.size(); j++) {
			std::cout << std::setw(10) << OverlapMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	return OverlapMatrix;
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

std::vector<std::vector<double>> Orbital_Overlaps::Calculate_Kinetic_Overlap()
{
	std::vector<std::vector<double>> OverlapMatrix;
	std::vector<std::vector<double>> Kinetic_Overlap_Matrix;

	//interations for Atoms 
	for (int Atom1 = 0; Atom1 < Total_Number_Of_Orbitals; Atom1++) {
		std::vector<double> row_Overlap;
		std::vector<double> row_Kinetic;
		for (int Atom2 = 0; Atom2 < Total_Number_Of_Orbitals; Atom2++) {

			//
			double Total_Overlap = 0;
			double Total_Kinetic_Overlap = 0;

			int orbital_Type_Atom1 = Orbital_Descriptions[Atom1][0] + Orbital_Descriptions[Atom1][1] + Orbital_Descriptions[Atom1][2];
			int orbital_Type_Atom2 = Orbital_Descriptions[Atom2][0] + Orbital_Descriptions[Atom2][1] + Orbital_Descriptions[Atom2][2];

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
						S[1][0] = -(Extended_Coordinates[Atom1][Coordinate] - ((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])));

						//calculate other Gaussian by regression.

						if (a > 1)
							for (int b = 2; b < a; b++) {
								S[b][0] = -(Extended_Coordinates[Atom1][Coordinate] - ((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])))
									*S[b - 1][0] + ((b - 1) / (2 * (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])))*S[b - 2][0];
							}

						//Transfer equation

						for (int c = 0; c < a - 1; c++) {
							for (int d = 1; d < a; d++) {
								S[c][d] = S[c + 1][d - 1] + (Extended_Coordinates[Atom1][Coordinate] - Extended_Coordinates[Atom2][Coordinate])*S[c][d - 1];
							}
						}

						//Store calculated Cartesian component in Sxyz matrix

						Sxyz[Coordinate] = S[(Orbital_Descriptions[Atom1][Coordinate])][(Orbital_Descriptions[Atom2][Coordinate])];


						//Inital conditions for kinetic intergral

						K[0][0] = 2 * Orbital_Coefficents[Atom1][Orbital1] * Orbital_Coefficents[Atom2][Orbital2] * S[1][1];

						for (int e = 1; e<a - 1; e++) {
							K[e][0] = -e*Orbital_Coefficents[Atom2][Orbital2] * S[e - 1][1] + 2 * Orbital_Coefficents[Atom1][Orbital1] * Orbital_Coefficents[Atom2][Orbital2] * S[e + 1][1];
						}


						for (int f = 1; f<a - 1; f++) {
							K[0][f] = -f*Orbital_Coefficents[Atom1][Orbital1] * S[1][f - 1] + 2 * Orbital_Coefficents[Atom1][Orbital1] * Orbital_Coefficents[Atom2][Orbital2] * S[1][f + 1];
						}


						//fills in rest of matrix
						for (int g = 1; g<a - 1; g++) {
							for (int h = 1; h<a - 1; h++) {
								K[g][h] = ((g*h*S[g - 1][h - 1])
									- (2 * g*Orbital_Coefficents[Atom2][Orbital2] * S[g - 1][h + 1])
									- (2 * h*Orbital_Coefficents[Atom1][Orbital1] * S[g + 1][h - 1])
									+ (4 * Orbital_Coefficents[Atom1][Orbital1] * Orbital_Coefficents[Atom2][Orbital2] * S[g + 1][h + 1])) / 2;
							}
						}


						Kxyz[Coordinate] = K[(Orbital_Descriptions[Atom1][Coordinate])][(Orbital_Descriptions[Atom2][Coordinate])];

					}

					//Calculate both kinetic and Overlap matrix values

					double EAB = exp(-1 * ((Orbital_Coefficents[Atom1][Orbital1] * Orbital_Coefficents[Atom2][Orbital2] / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2]))*
						((Extended_Coordinates[Atom1][0] - Extended_Coordinates[Atom2][0])*(Extended_Coordinates[Atom1][0] - Extended_Coordinates[Atom2][0]) + (Extended_Coordinates[Atom1][1] - Extended_Coordinates[Atom2][1])*(Extended_Coordinates[Atom1][1] - Extended_Coordinates[Atom2][1])
							+ (Extended_Coordinates[Atom1][2] - Extended_Coordinates[Atom2][2])*(Extended_Coordinates[Atom1][2] - Extended_Coordinates[Atom2][2]))));

					double Normanalisation_Const = Norm(Orbital_Coefficents[Atom1][Orbital1], orbital_Type_Atom1)*Norm(Orbital_Coefficents[Atom2][Orbital2], orbital_Type_Atom2)*Prime_Cutoffs[Atom1][Orbital1] * Prime_Cutoffs[Atom2][Orbital2];

					//Noralised Overlap matrix

					double Overlap = EAB*(pow((M_PI / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])), 1.5))*Sxyz[0] * Sxyz[1] * Sxyz[2];

					double Norm_Overlap = Normanalisation_Const*Overlap;

					//Total Overlap comprises of the sumation of 9 seperate normalised intergrals using the STO-3G Basis set
					Total_Overlap += Norm_Overlap;


					//Normalised Kinetic Matrix
					double Kinetic_Overlap = EAB*(pow((M_PI / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])), 1.5))*((Kxyz[0] * Sxyz[1] * Sxyz[2]) + (Sxyz[0] * Kxyz[1] * Sxyz[2]) + (Sxyz[0] * Sxyz[1] * Kxyz[2]));

					double Norm_Kinetic_Overlap = Normanalisation_Const*Kinetic_Overlap;
					Total_Kinetic_Overlap += Norm_Kinetic_Overlap;

				}
			}

			row_Overlap.push_back(Total_Overlap);
			row_Kinetic.push_back(Total_Kinetic_Overlap);
			/*OverlapMatrix[Atom1][Atom2] = Total_Overlap;
			Kinetic_Overlap_Matrix[Atom1][Atom2] = Total_Kinetic_Overlap;*/
		}
		OverlapMatrix.push_back(row_Overlap);
		Kinetic_Overlap_Matrix.push_back(row_Kinetic);
	}

	//clean up matrix ie. number smaller than 1e-8 = 0
	for (int i = 0; i < Total_Number_Of_Orbitals; i++) {
		for (int j = 0; j < Total_Number_Of_Orbitals; j++) {
			if (OverlapMatrix[i][j] < 1.0e-8 && OverlapMatrix[i][j]> -1e-8)
				OverlapMatrix[i][j] = 0;
		}
	}

	for (int i = 0; i < Total_Number_Of_Orbitals; i++) {
		for (int j = 0; j < Total_Number_Of_Orbitals; j++) {
			if (Kinetic_Overlap_Matrix[i][j] < 1.0e-8 && Kinetic_Overlap_Matrix[i][j]> -1e-8)
				Kinetic_Overlap_Matrix[i][j] = 0;
		}
	}

	//Temporary print function to display Overlap Matrix
	std::cout << std::defaultfloat << std::endl << "---Overlap Matrix---" << std::endl << std::endl;
	for (int i = 0; i < Total_Number_Of_Orbitals; i++) {
		for (int j = 0; j < Total_Number_Of_Orbitals; j++) {
			std::cout << std::setw(10) << OverlapMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;
	std::cout << std::defaultfloat << std::endl << "---Kinetic Overlap Matrix---" << std::endl << std::endl;
	for (int i = 0; i < Total_Number_Of_Orbitals; i++) {
		for (int j = 0; j < Total_Number_Of_Orbitals; j++) {
			std::cout << std::setw(10) << Kinetic_Overlap_Matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	return Kinetic_Overlap_Matrix;
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
double Orbital_Overlaps::abscissa(int n, int i)
{
	double A = (n + 1.0 - 2.0*i) / (n + 1.0) + (2.0 / M_PI)*(1.0 + (2.0 / 3.0)*pow(sin((i*M_PI) / (n + 1.0)), 2))*cos(i*M_PI / (n + 1.0))*sin(i*M_PI / (n + 1.0));
	return A;
}

double Orbital_Overlaps::Omega(int n, int i)
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

double Orbital_Overlaps::Fx(double X, int Atom1, int Atom2, int rr, int Orbital1, int Orbital2)
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
		S[1][0] = -1.0*(Extended_Coordinates[Atom1][Coordinate] - ((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2]))
			+ pow(((X + 1) / 2.0), 2) * (((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])) - RR[rr][Coordinate]));

		//Recurance index

		if (a>1)
			for (int b = 2; b<a; b++) {
				S[b][0] = -1.0*(Extended_Coordinates[Atom1][Coordinate] - ((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2]))
					+ pow(((X + 1) / 2.0), 2) * (((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2])) - RR[rr][Coordinate]))
					*S[b - 1][0] + ((b - 1) / (2 * (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2]))) * (1 - pow(((X + 1) / 2.0), 2)) *S[b - 2][0];
			}

		//Transfere Equation

		for (int c = 0; c<a - 1; c++) {
			for (int d = 1; d<a; d++) {
				S[c][d] = S[c + 1][d - 1] + (Extended_Coordinates[Atom1][Coordinate] - Extended_Coordinates[Atom2][Coordinate]) * S[c][d - 1];
			}
		}

		//Store calculated Cartesian component in Sxyz matrix

		Sxyz[Coordinate] = S[(Orbital_Descriptions[Atom1][Coordinate])][(Orbital_Descriptions[Atom2][Coordinate])];
	}

	double Dotproduct = 0;

	for (int Coordinate = 0; Coordinate < 3; Coordinate++) {

		double Coord_Value = ((Orbital_Coefficents[Atom1][Orbital1] * Extended_Coordinates[Atom1][Coordinate] + Orbital_Coefficents[Atom2][Orbital2] * Extended_Coordinates[Atom2][Coordinate]) / (Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2]) - RR[rr][Coordinate]);
		Coord_Value *= Coord_Value;
		Dotproduct += Coord_Value;
	}

	//What we are intergrating
	double Fx = 0.5 * (exp(-1.0*((Orbital_Coefficents[Atom1][Orbital1] + Orbital_Coefficents[Atom2][Orbital2]) * (pow(((X + 1) / 2.0), 2)) * Dotproduct)) *Sxyz[0] * Sxyz[1] * Sxyz[2]);

	return Fx;
}

double Orbital_Overlaps::CHP(double eps, double M, int Atom1, int Atom2, int RR, int Alpha, int Beta)

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
	double q = (Orbital_Overlaps::Fx(Orbital_Overlaps::abscissa(2.0, 1.0), Atom1, Atom2, RR, Alpha, Beta) + Orbital_Overlaps::Fx(-Orbital_Overlaps::abscissa(2.0, 1.0), Atom1, Atom2, RR, Alpha, Beta))*Orbital_Overlaps::Omega(2.0, 1.0);
	double p = Orbital_Overlaps::Fx(0.0, Atom1, Atom2, RR, Alpha, Beta);
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
				chp = chp + (Orbital_Overlaps::Fx(-xp, Atom1, Atom2, RR, Alpha, Beta) + Orbital_Overlaps::Fx(xp, Atom1, Atom2, RR, Alpha, Beta))*pow(S, 4);
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

	//Calculate and normalise overlap for a given set of 2 atoms -- sum of a total of 9 seperate orbital combinations
	double EAB = exp(-1 * ((Orbital_Coefficents[Atom1][Alpha] * Orbital_Coefficents[Atom2][Beta] / (Orbital_Coefficents[Atom1][Alpha] + Orbital_Coefficents[Atom2][Beta]))*
		((Extended_Coordinates[Atom1][0] - Extended_Coordinates[Atom2][0])*(Extended_Coordinates[Atom1][0] - Extended_Coordinates[Atom2][0]) + (Extended_Coordinates[Atom1][1] - Extended_Coordinates[Atom2][1])*(Extended_Coordinates[Atom1][1] - Extended_Coordinates[Atom2][1])
			+ (Extended_Coordinates[Atom1][2] - Extended_Coordinates[Atom2][2])*(Extended_Coordinates[Atom1][2] - Extended_Coordinates[Atom2][2]))));

	//Calculates resultant of intergration -- Un-normalised!!
	double result = EAB * ((2 * M_PI) / (Orbital_Coefficents[Atom1][Alpha] + Orbital_Coefficents[Atom2][Beta])) * chp;

	return result;
}

/*
#############################################################################################
## Driver function to pass in relavate combinations of atoms values to calculate OEI       ##
#############################################################################################
*/

std::vector<std::vector<double>> Orbital_Overlaps::Calculate_OneElectron_Overlap()
{
	std::vector<std::vector<double>> OverlapMatrix;

	//interations for Atoms 
	for (int Atom1 = 0; Atom1 < Total_Number_Of_Orbitals; Atom1++) {
		std::vector<double> row;
		for (int Atom2 = 0; Atom2 < Total_Number_Of_Orbitals; Atom2++) {

			double RR_Total_OEI[3] = { 0,0,0 };
			double new_Total_OEI = 0;
			int orbital_Type_Atom1 = Orbital_Descriptions[Atom1][0] + Orbital_Descriptions[Atom1][1] + Orbital_Descriptions[Atom1][2];
			int orbital_Type_Atom2 = Orbital_Descriptions[Atom2][0] + Orbital_Descriptions[Atom2][1] + Orbital_Descriptions[Atom2][2];

			//RR running from 1-3 
			for (int RR = 0; RR < 3; RR++) {

				double Total_OEI = 0;

				//Orbital indicates which set of orbital coeffcients and cut-offs are used
				for (int Orbital1 = 0; Orbital1 < 3; Orbital1++) {

					for (int Orbital2 = 0; Orbital2 < 3; Orbital2++) {


						double Normed = Orbital_Overlaps::Norm(Orbital_Coefficents[Atom1][Orbital1], orbital_Type_Atom1)*Norm(Orbital_Coefficents[Atom2][Orbital2], orbital_Type_Atom2)*Prime_Cutoffs[Atom1][Orbital1] * Prime_Cutoffs[Atom2][Orbital2];

						double OEI = Normed * Orbital_Overlaps::CHP(10e-10, 50000, Atom1, Atom2, RR, Orbital1, Orbital2);

						Total_OEI += OEI;
					}
				}
				RR_Total_OEI[RR] = Total_OEI*-1.0*Z[RR];
			}

			double r = RR_Total_OEI[0] + RR_Total_OEI[1] + RR_Total_OEI[2];
			row.push_back(r);

		}
		OverlapMatrix.push_back(row);
	}

	//Temporary print function to display Overlap Matrix
	std::cout << std::endl;
	std::cout << std::defaultfloat << std::endl << "---One Electron Intergrals---" << std::endl << std::endl;
	for (int i = 0; i < Total_Number_Of_Orbitals; i++) {
		for (int j = 0; j < Total_Number_Of_Orbitals; j++) {
			std::cout << std::setw(10) << OverlapMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	return OverlapMatrix;
}