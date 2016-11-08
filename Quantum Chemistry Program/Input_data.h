#pragma once

#include <iostream>       // STD C++ Lib
#include <fstream>		  // read and wirte to/from file
#include <vector>		  // allows use of vectors - varible size arrys
#include <algorithm>	  // Sort function
#include <string>         // for use of character strings
#include <istream>		  // for getline

class Input_data
{
public:

	int Get_number_of_Atoms(std::string Filename);
	std::vector<std::vector<double>> Load_Input(std::string Filename, int Number_of_Atoms);
	
	Input_data();
	~Input_data();

private:
	int Line_number;
	int Return_Line_Number(std::string Filename, std::string Search_Phrase);

};

