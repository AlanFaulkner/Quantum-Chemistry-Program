#include "Outputs.h"

//Public Functions.

void Outputs::Print_Geometry() {

	Input Data;
	int Number_of_atoms = Data.Get_number_of_Atoms("Input.dat");
	std::vector<std::vector<double>> Coordinates = Data.Get_Coordinates("Input.dat", Number_of_atoms);
	int Destination = Data.Get_print_level("input.dat");

	if (Destination == 0) {
		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##    --- Molecular Gemotery ---    ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;
		Outputs::Print_Bonds(Coordinates, Number_of_atoms, Destination);
		Outputs::Print_Angles(Coordinates, Number_of_atoms, Destination);
		Outputs::Print_OOP(Coordinates, Number_of_atoms, Destination);
		Outputs::Print_Torsions(Coordinates, Number_of_atoms, Destination);
		std::cout << std::defaultfloat << std::endl;
	}

	if (Destination == 1) {
		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##    --- Molecular Gemotery ---    ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;
		Outputs::Print_Bonds(Coordinates, Number_of_atoms, Destination);
		Outputs::Print_Angles(Coordinates, Number_of_atoms, Destination);
		Outputs::Print_OOP(Coordinates, Number_of_atoms, Destination);
		std::cout << std::defaultfloat << std::endl;
		Output.close();
	}

	if (Destination == 2) {
		std::ofstream Output("Output.txt", std::ios::app);
		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##    --- Molecular Gemotery ---    ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;
		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##    --- Molecular Gemotery ---    ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;
		Outputs::Print_Bonds(Coordinates, Number_of_atoms, Destination);
		Outputs::Print_Angles(Coordinates, Number_of_atoms, Destination);
		Outputs::Print_OOP(Coordinates, Number_of_atoms, Destination);
		std::cout << std::defaultfloat << std::endl;
		Output.close();
	}
	
}

void Outputs::Print_molecular_info()
{
	Input Data;
	int NoAtoms = Data.Get_number_of_Atoms("Input.dat");
	std::vector<std::vector<double>> Coords = Data.Get_Coordinates("Input.dat", NoAtoms);
	std::vector<double> CoM = Data.Calculate_centre_of_mass(Coords, NoAtoms);
	int Destination = Data.Get_print_level("input.dat");

	if (Destination == 0) {

		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##        --- Input Data ---        ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		std::cout << "The number of atoms in the molecule is: " << NoAtoms << std::endl << std::endl;
		std::cout << "The centre of mass for the molecule is:   X = " << CoM[0] << "   Y = " << CoM[1] << "   Z = " << CoM[2] << std::endl << std::endl;
		std::cout << "Input coordinates have been translated to reflect the centre of mass!" << std::endl << std::endl;
		Outputs::Print_Coordinates(Coords, Destination);
		std::cout << std::endl;
		Outputs::Print_Moment_of_inertia(Coords, NoAtoms, Destination);
	}

	else if (Destination == 1) {
		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##        --- Input Data ---        ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		Output << "The number of atoms in the molecule is: " << NoAtoms << std::endl << std::endl;
		Output << "The centre of mass for the molecule is:   X = " << CoM[0] << "   Y = " << CoM[1] << "   Z = " << CoM[2] << std::endl << std::endl;
		Output << "Input coordinates have been translated to reflect the centre of mass!" << std::endl << std::endl;
		Outputs::Print_Coordinates(Coords, Destination);
		Output << std::endl;
		Outputs::Print_Moment_of_inertia(Coords, NoAtoms, Destination);
		Output.close();
	}

	else if (Destination == 2) {

		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##        --- Input Data ---        ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		std::cout << "The number of atoms in the molecule is: " << NoAtoms << std::endl << std::endl;
		std::cout << "The centre of mass for the molecule is:   X = " << CoM[0] << "   Y = " << CoM[1] << "   Z = " << CoM[2] << std::endl << std::endl;
		std::cout << "Input coordinates have been translated to reflect the centre of mass!" << std::endl << std::endl;
		Outputs::Print_Coordinates(Coords, Destination);
		std::cout << std::endl;
		Outputs::Print_Moment_of_inertia(Coords, NoAtoms, Destination);

		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##        --- Input Data ---        ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		Output << "The number of atoms in the molecule is: " << NoAtoms << std::endl << std::endl;
		Output << "The centre of mass for the molecule is:   X = " << CoM[0] << "   Y = " << CoM[1] << "   Z = " << CoM[2] << std::endl << std::endl;
		Output << "Input coordinates have been translated to reflect the centre of mass!" << std::endl << std::endl;
		Output.close();
	}

}

//Private Functions.
void Outputs::Print_Coordinates(std::vector<std::vector<double>> Coordinates, int Destination)
{
	Input data;
	int NoAtoms = data.Get_number_of_Atoms("Input.dat");
	//Print to screen.
	if (Destination == 0) {
		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##  --- Molecular Coordinates ---   ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;
		for (int i = 0; i < NoAtoms; i++) {
			std::cout << std::setw(3) << Coordinates[i][0] << " ";
			for (int j = 1; j < 4; j++) {
				std::cout << std::fixed << std::setw(10) << std::setprecision(7) << std::right << Coordinates[i][j] << " ";
			}
			std::cout << std::defaultfloat << std::endl;
		}
	}

	//Print to file called "Output.txt".
	else if (Destination == 1) {
		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##  --- Molecular Coordinates ---   ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;
		for (int i = 0; i < NoAtoms; i++) {
			Output << std::setw(3) << Coordinates[i][0] << " ";
			for (int j = 1; j < 4; j++) {
				Output << std::fixed << std::setw(10) << std::setprecision(7) << std::right << Coordinates[i][j] << " ";
			}
			Output << std::defaultfloat << std::endl;
		}
		Output.close();
	}

	//Print to both screen and file called "Output.txt".
	else if (Destination == 2) {

		std::ofstream Output("Output.txt", std::ios::app);

		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##  --- Molecular Coordinates ---   ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;
		for (int i = 0; i < NoAtoms; i++) {
			std::cout << std::setw(3) << Coordinates[i][0] << " ";
			for (int j = 1; j < 4; j++) {
				std::cout << std::fixed << std::setw(10) << std::setprecision(7) << std::right << Coordinates[i][j] << " ";
			}
			std::cout << std::defaultfloat << std::endl;
		}

		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##  --- Molecular Coordinates ---   ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;
		for (int i = 0; i < NoAtoms; i++) {
			Output << std::setw(3) << Coordinates[i][0] << " ";
			for (int j = 1; j < 4; j++) {
				Output << std::fixed << std::setw(10) << std::setprecision(7) << std::right << Coordinates[i][j] << " ";
			}
			Output <<std::defaultfloat << std::endl; 
		}
		Output.close();
	}
}

void Outputs::Print_Moment_of_inertia(std::vector<std::vector<double>> Coordinates, int NoAtoms, int Destination)
{
	Input data;
	std::vector<std::vector<double>> I = data.Get_Moments_of_inertia(Coordinates, NoAtoms);

	if (Destination == 0) {

		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##    --- Moments of Inertia ---    ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		std::cout << "Moment of inertia tensor (amu bohr^2)" << std::endl <<
			"-------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				std::cout << std::fixed << std::setw(10) << std::setprecision(7) << std::right << I[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl << std::endl;

		std::cout << "Principal moments of inertia (amu * bohr^2)" << std::endl <<
			"-------------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			std::cout << I[3][i] << "  ";
		std::cout << std::endl << std::endl << std::endl;

		double conv = 0.529177249 * 0.529177249;
		std::cout << "Principal moments of inertia (amu * AA^2)" << std::endl <<
			"-----------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			std::cout << I[3][i] * conv << "  ";
		std::cout << std::endl << std::endl << std::endl;

		double conv2 = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
		std::cout << "Principal moments of inertia (g * cm^2)" << std::endl <<
			"---------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			std::cout << I[3][i] * conv2 << "  ";
		std::cout << std::endl << std::endl;

		// classify the rotor 

		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##      --- Type Of Rotor ---       ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		if (NoAtoms == 2) std::cout << "Molecule is diatomic." << std::endl << std::endl;
		else if (I[3][0] < 1e-4) std::cout << "Molecule is linear." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) < 1e-4) && (fabs(I[3][1] - I[3][2]) < 1e-4))
			std::cout << "Molecule is a spherical top." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) < 1e-4) && (fabs(I[3][1] - I[3][2]) > 1e-4))
			std::cout << "Molecule is an oblate symmetric top." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) > 1e-4) && (fabs(I[3][1] - I[3][2]) < 1e-4))
			std::cout << "Molecule is a prolate symmetric top." << std::endl << std::endl;
		else std::cout << "Molecule is an asymmetric top." << std::endl << std::endl;
	}

	else if (Destination == 1) {

		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##    --- Moments of Inertia ---    ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		Output << "Moment of inertia tensor (amu bohr^2)" << std::endl <<
			"-------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Output << std::fixed << std::setw(10) << std::setprecision(7) << std::right << I[i][j] << " ";
			}
			Output << std::endl;
		}
		Output << std::endl << std::endl;

		Output << "Principal moments of inertia (amu * bohr^2)" << std::endl <<
			"-------------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			Output << I[3][i] << "  ";
		Output << std::endl << std::endl << std::endl;

		double conv = 0.529177249 * 0.529177249;
		Output << "Principal moments of inertia (amu * AA^2)" << std::endl <<
			"-----------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			Output << I[3][i] * conv << "  ";
		Output << std::endl << std::endl << std::endl;

		double conv2 = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
		Output << "Principal moments of inertia (g * cm^2)" << std::endl <<
			"---------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			Output << I[3][i] * conv2 << "  ";
		Output << std::endl << std::endl;

		// classify the rotor 

		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##      --- Type Of Rotor ---       ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		if (NoAtoms == 2) Output << "Molecule is diatomic." << std::endl << std::endl;
		else if (I[3][0] < 1e-4) Output << "Molecule is linear." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) < 1e-4) && (fabs(I[3][1] - I[3][2]) < 1e-4))
			Output << "Molecule is a spherical top." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) < 1e-4) && (fabs(I[3][1] - I[3][2]) > 1e-4))
			Output << "Molecule is an oblate symmetric top." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) > 1e-4) && (fabs(I[3][1] - I[3][2]) < 1e-4))
			Output << "Molecule is a prolate symmetric top." << std::endl << std::endl;
		else Output << "Molecule is an asymmetric top." << std::endl << std::endl;
		Output.close();
	}

	else if (Destination == 2) {

		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##    --- Moments of Inertia ---    ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		std::cout << "Moment of inertia tensor (amu bohr^2)" << std::endl <<
			"-------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				std::cout << std::fixed << std::setw(10) << std::setprecision(7) << std::right << I[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl << std::endl;

		std::cout << "Principal moments of inertia (amu * bohr^2)" << std::endl <<
			"-------------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			std::cout << I[3][i] << "  ";
		std::cout << std::endl << std::endl << std::endl;

		double conv = 0.529177249 * 0.529177249;
		std::cout << "Principal moments of inertia (amu * AA^2)" << std::endl <<
			"-----------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			std::cout << I[3][i] * conv << "  ";
		std::cout << std::endl << std::endl << std::endl;

		double conv2 = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
		std::cout << "Principal moments of inertia (g * cm^2)" << std::endl <<
			"---------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			std::cout << I[3][i] * conv2 << "  ";
		std::cout << std::endl << std::endl;

		// classify the rotor 

		std::cout << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##      --- Type Of Rotor ---       ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		if (NoAtoms == 2) std::cout << "Molecule is diatomic." << std::endl << std::endl;
		else if (I[3][0] < 1e-4) std::cout << "Molecule is linear." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) < 1e-4) && (fabs(I[3][1] - I[3][2]) < 1e-4))
			std::cout << "Molecule is a spherical top." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) < 1e-4) && (fabs(I[3][1] - I[3][2]) > 1e-4))
			std::cout << "Molecule is an oblate symmetric top." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) > 1e-4) && (fabs(I[3][1] - I[3][2]) < 1e-4))
			std::cout << "Molecule is a prolate symmetric top." << std::endl << std::endl;
		else std::cout << "Molecule is an asymmetric top." << std::endl << std::endl;

		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##    --- Moments of Inertia ---    ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		Output << "Moment of inertia tensor (amu bohr^2)" << std::endl <<
			"-------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Output << std::fixed << std::setw(10) << std::setprecision(7) << std::right << I[i][j] << " ";
			}
			Output << std::endl;
		}
		Output << std::endl << std::endl;

		Output << "Principal moments of inertia (amu * bohr^2)" << std::endl <<
			"-------------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			Output << I[3][i] << "  ";
		Output << std::endl << std::endl << std::endl;

		Output << "Principal moments of inertia (amu * AA^2)" << std::endl <<
			"-----------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			Output << I[3][i] * conv << "  ";
		Output << std::endl << std::endl << std::endl;

		Output << "Principal moments of inertia (g * cm^2)" << std::endl <<
			"---------------------------------------" << std::endl << std::endl;
		for (int i = 0; i < 3; i++)
			Output << I[3][i] * conv2 << "  ";
		Output << std::endl << std::endl;

		// classify the rotor 

		Output << std::endl << "######################################" << std::endl
			<< "##                                  ##" << std::endl
			<< "##      --- Type Of Rotor ---       ##" << std::endl
			<< "##                                  ##" << std::endl
			<< "######################################" << std::endl << std::endl;

		if (NoAtoms == 2) Output << "Molecule is diatomic." << std::endl << std::endl;
		else if (I[3][0] < 1e-4) Output << "Molecule is linear." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) < 1e-4) && (fabs(I[3][1] - I[3][2]) < 1e-4))
			Output << "Molecule is a spherical top." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) < 1e-4) && (fabs(I[3][1] - I[3][2]) > 1e-4))
			Output << "Molecule is an oblate symmetric top." << std::endl << std::endl;
		else if ((fabs(I[3][0] - I[3][1]) > 1e-4) && (fabs(I[3][1] - I[3][2]) < 1e-4))
			Output << "Molecule is a prolate symmetric top." << std::endl << std::endl;
		else Output << "Molecule is an asymmetric top." << std::endl << std::endl;
		Output.close();
	}
}

double Outputs::bond(std::vector<std::vector<double>>Input_Coordinates, int i, int j)
{
	return sqrt(((Input_Coordinates[i][1] - Input_Coordinates[j][1])*(Input_Coordinates[i][1] - Input_Coordinates[j][1]))
		+ ((Input_Coordinates[i][2] - Input_Coordinates[j][2])*(Input_Coordinates[i][2] - Input_Coordinates[j][2]))
		+ ((Input_Coordinates[i][3] - Input_Coordinates[j][3])*(Input_Coordinates[i][3] - Input_Coordinates[j][3])));
}

double Outputs::unit(std::vector<std::vector<double>>Input_Coordinates, int cart, int a, int b)
{
	return -(Input_Coordinates[a][cart] - Input_Coordinates[b][cart]) / Outputs::bond(Input_Coordinates, a, b);
}

double Outputs::angle(std::vector<std::vector<double>>Input_Coordinates, int a, int b, int c)
{
	return acos(Outputs::unit(Input_Coordinates, 1, b, a) * Outputs::unit(Input_Coordinates, 1, b, c)
		+ Outputs::unit(Input_Coordinates, 2, b, a) * Outputs::unit(Input_Coordinates, 2, b, c)
		+ Outputs::unit(Input_Coordinates, 3, b, a) * Outputs::unit(Input_Coordinates, 3, b, c));
}

void Outputs::Print_Bonds(std::vector<std::vector<double>> Coordinates, int Number_of_atoms, int Destination) {

	if (Destination == 0) {
		std::cout << std::endl << " -- Bond Lengths (A) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				std::cout << std::setw(3) << i << " " << j << " ";
				double bond_length = Outputs::bond(Coordinates, i, j);
				std::cout << std::fixed << std::setw(10) << std::setprecision(6) << bond_length << std::endl;
			}
		}
	}

	else if (Destination == 1) {
		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << " -- Bond Lengths (A) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				Output << std::setw(3) << i << " " << j << " ";
				double bond_length = Outputs::bond(Coordinates, i, j);
				Output << std::fixed << std::setw(10) << std::setprecision(6) << bond_length << std::endl;
			}
		}
		Output.close();
	}

	else if (Destination == 2) {
		std::ofstream Output("Output.txt", std::ios::app);
		std::cout << std::endl << " -- Bond Lengths (A) --" << std::endl << std::endl;
		Output << std::endl << " -- Bond Lengths (A) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				std::cout << std::setw(3) << i << " " << j << " ";
				Output << std::setw(3) << i << " " << j << " ";
				double bond_length = Outputs::bond(Coordinates, i, j);
				std::cout << std::fixed << std::setw(10) << std::setprecision(6) << bond_length << std::endl;
				Output << std::fixed << std::setw(10) << std::setprecision(6) << bond_length << std::endl;
			}
		}
		Output.close();
	}
}

void Outputs::Print_Angles(std::vector<std::vector<double>> Coordinates, int Number_of_atoms, int Destination) {

	if (Destination == 0) {
		std::cout << std::endl << "-- Bond Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < j; k++) {
					if (bond(Coordinates, i, j) < 4 && bond(Coordinates, j, k) < 4) {
						std::cout << std::setw(3) << i << " " << j << " " << k << " ";
						double Angle = Outputs::angle(Coordinates, i, j, k)*(180.0 / acos(-1.0));
						std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right  << Angle << std::endl;
					}
				}
			}
		}
	}

	else if (Destination == 1) {
		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "-- Bond Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < j; k++) {
					if (bond(Coordinates, i, j) < 4 && bond(Coordinates, j, k) < 4) {
						Output << std::setw(3) << i << " " << j << " " << k << " ";
						double Angle = Outputs::angle(Coordinates, i, j, k)*(180.0 / acos(-1.0));
						Output << std::fixed << std::setw(7) << std::setprecision(2) << std::right << Angle << std::endl;
					}
				}
			}
		}
		Output.close();
	}

	else if (Destination == 2) {
		std::ofstream Output;
		Output.open("Output.txt", std::ios::app);
		std::cout << std::endl << "-- Bond Angles (deg) --" << std::endl << std::endl;
		Output << std::endl << "-- Bond Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < j; k++) {
					if (bond(Coordinates, i, j) < 4 && bond(Coordinates, j, k) < 4) {
						std::cout << std::setw(3) << i << " " << j << " " << k << " ";
						Output << std::setw(3) << i << " " << j << " " << k << " ";
						double Angle = Outputs::angle(Coordinates, i, j, k)*(180.0 / acos(-1.0));
						std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << Angle << std::endl;
						Output << std::fixed << std::setw(7) << std::setprecision(2) << std::right << Angle << std::endl;
					}
				}
			}
		}
		Output.close();
	}
}

void Outputs::Print_OOP(std::vector<std::vector<double>> Coordinates, int Number_of_atoms, int Destination) {
	if (Destination == 0) {
		std::cout << std::endl << "-- Out Of Plane Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int k = 0; k < Number_of_atoms; k++) {
				for (int j = 0; j < Number_of_atoms; j++) {
					for (int l = 0; l < j; l++) {
						if (i != j && i != k && i != l && j != k && k != l && Outputs::bond(Coordinates, i, k) < 4 && Outputs::bond(Coordinates, k, j) < 4 && Outputs::bond(Coordinates, k, l) < 4) {
							std::cout << std::setw(3) << i << " " << k << " " << j << " " << l << " ";

							double ox = (Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 3, k, l) - Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 2, k, l));
							double oy = (Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 1, k, l) - Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 3, k, l));
							double oz = (Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 2, k, l) - Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 1, k, l));

							double oxx = ox*Outputs::unit(Coordinates, 1, k, i);
							double oyy = oy*Outputs::unit(Coordinates, 2, k, i);
							double ozz = oz*Outputs::unit(Coordinates, 3, k, i);

							double theta = (oxx + oyy + ozz) / sin(Outputs::angle(Coordinates, j, k, l));

							if (theta < -1.0) { theta = asin(-1.0); }
							else if (theta > 1.0) { theta = asin(1.0); }
							else { theta = asin(theta); }

							theta = theta*(180.0 / acos(-1.0));
							std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << theta << std::endl;
						}
					}
				}
			}
		}
	}

	else if (Destination == 1) {
		std::ofstream Output;
		Output.open("Output.txt", std::ios::app);
		Output << std::endl << "-- Out Of Plane Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int k = 0; k < Number_of_atoms; k++) {
				for (int j = 0; j < Number_of_atoms; j++) {
					for (int l = 0; l < j; l++) {
						if (i != j && i != k && i != l && j != k && k != l && Outputs::bond(Coordinates, i, k) < 4 && Outputs::bond(Coordinates, k, j) < 4 && Outputs::bond(Coordinates, k, l) < 4) {
							Output << std::setw(3) << i << " " << k << " " << j << " " << l << " ";

							double ox = (Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 3, k, l) - Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 2, k, l));
							double oy = (Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 1, k, l) - Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 3, k, l));
							double oz = (Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 2, k, l) - Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 1, k, l));

							double oxx = ox*Outputs::unit(Coordinates, 1, k, i);
							double oyy = oy*Outputs::unit(Coordinates, 2, k, i);
							double ozz = oz*Outputs::unit(Coordinates, 3, k, i);

							double theta = (oxx + oyy + ozz) / sin(Outputs::angle(Coordinates, j, k, l));

							if (theta < -1.0) { theta = asin(-1.0); }
							else if (theta > 1.0) { theta = asin(1.0); }
							else { theta = asin(theta); }

							theta = theta*(180.0 / acos(-1.0));
							Output << std::fixed << std::setw(7) << std::setprecision(2) << std::right << theta << std::endl;
						}
					}
				}
			}
		}
		Output.close();
	}

	else if (Destination == 2) {
		std::ofstream Output;
		Output.open("Output.txt", std::ios::app);
		std::cout << std::endl << "-- Out Of Plane Angles (deg) --" << std::endl << std::endl;
		Output << std::endl << "-- Out Of Plane Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int k = 0; k < Number_of_atoms; k++) {
				for (int j = 0; j < Number_of_atoms; j++) {
					for (int l = 0; l < j; l++) {
						if (i != j && i != k && i != l && j != k && k != l && Outputs::bond(Coordinates, i, k) < 4 && Outputs::bond(Coordinates, k, j) < 4 && Outputs::bond(Coordinates, k, l) < 4) {
							std::cout << std::setw(3) << i << " " << k << " " << j << " " << l << " ";
							Output << std::setw(3) << i << " " << k << " " << j << " " << l << " ";

							double ox = (Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 3, k, l) - Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 2, k, l));
							double oy = (Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 1, k, l) - Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 3, k, l));
							double oz = (Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 2, k, l) - Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 1, k, l));

							double oxx = ox*Outputs::unit(Coordinates, 1, k, i);
							double oyy = oy*Outputs::unit(Coordinates, 2, k, i);
							double ozz = oz*Outputs::unit(Coordinates, 3, k, i);

							double theta = (oxx + oyy + ozz) / sin(Outputs::angle(Coordinates, j, k, l));

							if (theta < -1.0) { theta = asin(-1.0); }
							else if (theta > 1.0) { theta = asin(1.0); }
							else { theta = asin(theta); }

							theta = theta*(180.0 / acos(-1.0));
							std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << theta << std::endl;
							Output << std::fixed << std::setw(7) << std::setprecision(2) << std::right << theta << std::endl;
						}
					}
				}
			}
		}
		Output.close();
	}
}

void Outputs::Print_Torsions(std::vector<std::vector<double>> Coordinates, int Number_of_atoms, int Destination) {
	if (Destination == 0) {
		std::cout << std::endl << "-- Torsion Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < j; k++) {
					for (int l = 0; l < k; l++) {
						if (Outputs::bond(Coordinates, i, j) < 4.0 && Outputs::bond(Coordinates, j, k) < 4.0 && Outputs::bond(Coordinates, k, l) < 4.0) {
							std::cout << std::setw(3) << i << " " << j << " " << k << " " << l << " ";
							double eabc_x = (Outputs::unit(Coordinates, 2, j, i) * Outputs::unit(Coordinates, 3, j, k) - Outputs::unit(Coordinates, 3, j, i) * Outputs::unit(Coordinates, 2, j, k));
							double eabc_y = (Outputs::unit(Coordinates, 3, j, i) * Outputs::unit(Coordinates, 1, j, k) - Outputs::unit(Coordinates, 1, j, i) * Outputs::unit(Coordinates, 3, j, k));
							double eabc_z = (Outputs::unit(Coordinates, 1, j, i) * Outputs::unit(Coordinates, 2, j, k) - Outputs::unit(Coordinates, 2, j, i) * Outputs::unit(Coordinates, 1, j, k));

							double ebcd_x = (Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 3, k, l) - Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 2, k, l));
							double ebcd_y = (Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 1, k, l) - Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 3, k, l));
							double ebcd_z = (Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 2, k, l) - Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 1, k, l));
							
							double exx = eabc_x * ebcd_x;
							double eyy = eabc_y * ebcd_y;
							double ezz = eabc_z * ebcd_z;

							double tau = (exx + eyy + ezz) / (sin(Outputs::angle(Coordinates, i, j, k)) * sin(Outputs::angle(Coordinates, j, k, l)));

							if (tau < -1.0) tau = acos(-1.0);
							else if (tau > 1.0) tau = acos(1.0);
							else tau = acos(tau);

							// Compute the sign of the torsion 
							double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
							double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
							double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
							double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
							cross_x /= norm;
							cross_y /= norm;
							cross_z /= norm;
							double sign = 1.0;
							double dot = cross_x*Outputs::unit(Coordinates, 0, j, k) + cross_y*Outputs::unit(Coordinates, 1, j, k) + cross_z*Outputs::unit(Coordinates, 2, j, k);
							if (dot < 0.0) sign = -1.0;

							tau = tau*sign*(180.0 / acos(-1.0));
							std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << tau << std::endl;
						}	
					}
				}
			}
		}
	}
	
	if (Destination == 1) {
		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "-- Torsion Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < j; k++) {
					for (int l = 0; l < k; l++) {
						if (Outputs::bond(Coordinates, i, j) < 4.0 && Outputs::bond(Coordinates, j, k) < 4.0 && Outputs::bond(Coordinates, k, l) < 4.0) {
							Output << std::setw(3) << i << " " << j << " " << k << " " << l << " ";
							double eabc_x = (Outputs::unit(Coordinates, 2, j, i) * Outputs::unit(Coordinates, 3, j, k) - Outputs::unit(Coordinates, 3, j, i) * Outputs::unit(Coordinates, 2, j, k));
							double eabc_y = (Outputs::unit(Coordinates, 3, j, i) * Outputs::unit(Coordinates, 1, j, k) - Outputs::unit(Coordinates, 1, j, i) * Outputs::unit(Coordinates, 3, j, k));
							double eabc_z = (Outputs::unit(Coordinates, 1, j, i) * Outputs::unit(Coordinates, 2, j, k) - Outputs::unit(Coordinates, 2, j, i) * Outputs::unit(Coordinates, 1, j, k));

							double ebcd_x = (Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 3, k, l) - Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 2, k, l));
							double ebcd_y = (Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 1, k, l) - Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 3, k, l));
							double ebcd_z = (Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 2, k, l) - Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 1, k, l));

							double exx = eabc_x * ebcd_x;
							double eyy = eabc_y * ebcd_y;
							double ezz = eabc_z * ebcd_z;

							double tau = (exx + eyy + ezz) / (sin(Outputs::angle(Coordinates, i, j, k)) * sin(Outputs::angle(Coordinates, j, k, l)));

							if (tau < -1.0) tau = acos(-1.0);
							else if (tau > 1.0) tau = acos(1.0);
							else tau = acos(tau);

							// Compute the sign of the torsion 
							double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
							double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
							double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
							double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
							cross_x /= norm;
							cross_y /= norm;
							cross_z /= norm;
							double sign = 1.0;
							double dot = cross_x*Outputs::unit(Coordinates, 0, j, k) + cross_y*Outputs::unit(Coordinates, 1, j, k) + cross_z*Outputs::unit(Coordinates, 2, j, k);
							if (dot < 0.0) sign = -1.0;

							tau = tau*sign*(180.0 / acos(-1.0));
							Output << std::fixed << std::setw(7) << std::setprecision(2) << std::right << tau << std::endl;
						}
					}
				}
			}
		}
		Output.close();
	}

	if (Destination == 2) {
		std::cout << std::endl << "-- Torsion Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < j; k++) {
					for (int l = 0; l < k; l++) {
						if (Outputs::bond(Coordinates, i, j) < 4.0 && Outputs::bond(Coordinates, j, k) < 4.0 && Outputs::bond(Coordinates, k, l) < 4.0) {
							std::cout << std::setw(3) << i << " " << j << " " << k << " " << l << " ";
							double eabc_x = (Outputs::unit(Coordinates, 2, j, i) * Outputs::unit(Coordinates, 3, j, k) - Outputs::unit(Coordinates, 3, j, i) * Outputs::unit(Coordinates, 2, j, k));
							double eabc_y = (Outputs::unit(Coordinates, 3, j, i) * Outputs::unit(Coordinates, 1, j, k) - Outputs::unit(Coordinates, 1, j, i) * Outputs::unit(Coordinates, 3, j, k));
							double eabc_z = (Outputs::unit(Coordinates, 1, j, i) * Outputs::unit(Coordinates, 2, j, k) - Outputs::unit(Coordinates, 2, j, i) * Outputs::unit(Coordinates, 1, j, k));

							double ebcd_x = (Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 3, k, l) - Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 2, k, l));
							double ebcd_y = (Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 1, k, l) - Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 3, k, l));
							double ebcd_z = (Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 2, k, l) - Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 1, k, l));

							double exx = eabc_x * ebcd_x;
							double eyy = eabc_y * ebcd_y;
							double ezz = eabc_z * ebcd_z;

							double tau = (exx + eyy + ezz) / (sin(Outputs::angle(Coordinates, i, j, k)) * sin(Outputs::angle(Coordinates, j, k, l)));

							if (tau < -1.0) tau = acos(-1.0);
							else if (tau > 1.0) tau = acos(1.0);
							else tau = acos(tau);

							// Compute the sign of the torsion 
							double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
							double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
							double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
							double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
							cross_x /= norm;
							cross_y /= norm;
							cross_z /= norm;
							double sign = 1.0;
							double dot = cross_x*Outputs::unit(Coordinates, 0, j, k) + cross_y*Outputs::unit(Coordinates, 1, j, k) + cross_z*Outputs::unit(Coordinates, 2, j, k);
							if (dot < 0.0) sign = -1.0;

							tau = tau*sign*(180.0 / acos(-1.0));
							std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << tau << std::endl;
						}
					}
				}
			}
		}

		std::ofstream Output("Output.txt", std::ios::app);
		Output << std::endl << "-- Torsion Angles (deg) --" << std::endl << std::endl;
		for (int i = 0; i < Number_of_atoms; i++) {
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < j; k++) {
					for (int l = 0; l < k; l++) {
						if (Outputs::bond(Coordinates, i, j) < 4.0 && Outputs::bond(Coordinates, j, k) < 4.0 && Outputs::bond(Coordinates, k, l) < 4.0) {
							Output << std::setw(3) << i << " " << j << " " << k << " " << l << " ";
							double eabc_x = (Outputs::unit(Coordinates, 2, j, i) * Outputs::unit(Coordinates, 3, j, k) - Outputs::unit(Coordinates, 3, j, i) * Outputs::unit(Coordinates, 2, j, k));
							double eabc_y = (Outputs::unit(Coordinates, 3, j, i) * Outputs::unit(Coordinates, 1, j, k) - Outputs::unit(Coordinates, 1, j, i) * Outputs::unit(Coordinates, 3, j, k));
							double eabc_z = (Outputs::unit(Coordinates, 1, j, i) * Outputs::unit(Coordinates, 2, j, k) - Outputs::unit(Coordinates, 2, j, i) * Outputs::unit(Coordinates, 1, j, k));

							double ebcd_x = (Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 3, k, l) - Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 2, k, l));
							double ebcd_y = (Outputs::unit(Coordinates, 3, k, j) * Outputs::unit(Coordinates, 1, k, l) - Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 3, k, l));
							double ebcd_z = (Outputs::unit(Coordinates, 1, k, j) * Outputs::unit(Coordinates, 2, k, l) - Outputs::unit(Coordinates, 2, k, j) * Outputs::unit(Coordinates, 1, k, l));

							double exx = eabc_x * ebcd_x;
							double eyy = eabc_y * ebcd_y;
							double ezz = eabc_z * ebcd_z;

							double tau = (exx + eyy + ezz) / (sin(Outputs::angle(Coordinates, i, j, k)) * sin(Outputs::angle(Coordinates, j, k, l)));

							if (tau < -1.0) tau = acos(-1.0);
							else if (tau > 1.0) tau = acos(1.0);
							else tau = acos(tau);

							// Compute the sign of the torsion 
							double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
							double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
							double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
							double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
							cross_x /= norm;
							cross_y /= norm;
							cross_z /= norm;
							double sign = 1.0;
							double dot = cross_x*Outputs::unit(Coordinates, 0, j, k) + cross_y*Outputs::unit(Coordinates, 1, j, k) + cross_z*Outputs::unit(Coordinates, 2, j, k);
							if (dot < 0.0) sign = -1.0;

							tau = tau*sign*(180.0 / acos(-1.0));
							Output << std::fixed << std::setw(7) << std::setprecision(2) << std::right << tau << std::endl;
						}
					}
				}
			}
		}
		Output.close();
	}
}

Outputs::Outputs()
{
}

Outputs::~Outputs()
{
}
