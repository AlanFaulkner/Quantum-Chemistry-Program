#include "Input.h"

int Input::Get_number_of_Atoms(std::string Filename)
{
	std::ifstream Input(Filename);
	int Number_of_Atoms;
	Input >> Number_of_Atoms;
	Input.close();
	return Number_of_Atoms;
}

int Input::Get_print_level(std::string Filename)
{
	int Line_number = Return_Line_Number(Filename, "Print Level");

	std::ifstream Input(Filename);
	std::string line;
	for (int i = 0; i <= Line_number; i++)
		getline(Input, line);

	int Print_level;
	Input >> Print_level;

	return Print_level;
}

std::vector<std::vector<double>> Input::Get_Moments_of_inertia(std::vector<std::vector<double>> Coordinates, int Number_of_atoms)
{
	Matrix I(3, 3);
	double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
	double Mass;
	for (int i = 0; i < Number_of_atoms; i++) {
		int Atom_Number = Coordinates[i][0];
		Mass = Input::Atomic_Masses[Atom_Number];
		a += Mass * (Coordinates[i][2] * Coordinates[i][2] + Coordinates[i][3] * Coordinates[i][3]);
		b += Mass * (Coordinates[i][1] * Coordinates[i][1] + Coordinates[i][3] * Coordinates[i][3]);
		c += Mass * (Coordinates[i][1] * Coordinates[i][1] + Coordinates[i][2] * Coordinates[i][2]);
		d += Mass * Coordinates[i][1] * Coordinates[i][2];
		e += Mass * Coordinates[i][1] * Coordinates[i][3];
		f += Mass * Coordinates[i][2] * Coordinates[i][3];
	}

	I(0, 0) = a;
	I(1, 1) = b;
	I(2, 2) = c;
	I(0, 1) = d;
	I(0, 2) = e;
	I(1, 2) = f;
	I(1, 0) = I(0, 1);
	I(2, 0) = I(0, 2);
	I(2, 1) = I(1, 2);

	// find the principal moments
	Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
	Matrix evecs = solver.eigenvectors();
	Matrix evals = solver.eigenvalues();

	std::vector<std::vector<double>> MoI;
	
	std::vector<double> eval;
	for (int i = 0; i < 3; i++) {
		std::vector<double> row;
		for (int j = 0; j < 3; j++) {
			row.push_back(I(i, j));
		}
		MoI.push_back(row);
	}

	for (int i = 0; i < 3; i++)
		eval.push_back(evals(i));
	MoI.push_back(eval);
	return MoI;
	
}

//Search for deseried phrase in a file and stores the line number on which it is found.
int Input::Return_Line_Number(std::string Filename, std::string Search_Phrase)
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

std::vector<double> Input::Calculate_centre_of_mass(std::vector<std::vector<double>> Coordinates, int Number_of_atoms)
{
	double M = 0.0;
	double xcm = 0.0;
	double ycm = 0.0;
	double zcm = 0.0;
	double mi;

	for (int i = 0; i < Number_of_atoms; i++) {
		int Atom_Number = Coordinates[i][0];
		M += Input::Atomic_Masses[Atom_Number];
		mi = Input::Atomic_Masses[Atom_Number];
		xcm += mi * Coordinates[i][1];
		ycm += mi * Coordinates[i][2];
		zcm += mi * Coordinates[i][3];
	}

	xcm /= M;
	ycm /= M;
	zcm /= M;

	std::vector<double> Translation_Coords;

	Translation_Coords.push_back(xcm);
	Translation_Coords.push_back(ycm);
	Translation_Coords.push_back(zcm);

	return Translation_Coords;
}

std::vector<std::vector<double>> Input::Get_Coordinates(std::string Filename, int Number_of_Atoms)
{
	int Line_number = Input::Return_Line_Number(Filename, "Cart");

	std::ifstream Input(Filename);
	std::string line;
	for (int i = 0; i <= Line_number; i++)
		getline(Input, line);

	//sets up and initalises a vector array for storing atom number, x y z coordinates -- varible size
	std::vector<std::vector <double>> Input_Data(Number_of_Atoms, std::vector<double>(4));

	for (int i = 0; i < Number_of_Atoms; i++)
		Input >> Input_Data[i][0] >> Input_Data[i][1] >> Input_Data[i][2] >> Input_Data[i][3];

	Input.close();

	//sorts atoms decending 
	std::sort(Input_Data.rbegin(), Input_Data.rend());

	std::vector<double> Translation_coords = Input::Calculate_centre_of_mass(Input_Data, Number_of_Atoms);

	std::vector<std::vector<double>> Mass_centered_Coordinates(Number_of_Atoms, std::vector<double>(4));
	for (int i = 0; i < Number_of_Atoms; i++) {
		Mass_centered_Coordinates[i][0] = Input_Data[i][0];
		Mass_centered_Coordinates[i][1] = Input_Data[i][1] - Translation_coords[0];
		Mass_centered_Coordinates[i][2] = Input_Data[i][2] - Translation_coords[1];
		Mass_centered_Coordinates[i][3] = Input_Data[i][3] - Translation_coords[2];
	}

	return Mass_centered_Coordinates;
}

Input::Input()
{
}


Input::~Input()
{
}
