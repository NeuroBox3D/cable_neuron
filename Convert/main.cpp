/*
 * main.cpp
 *
 *  Created on: 30.09.2013
 *      Author: unheiliger
 */

#include "Converter.h"
#include <iostream>



int main(int argn, char* argv[]) {
if (argn == 2)
{
	std::cout << "Programm gestartet" << std::endl;
	Converter test = Converter();
	std::vector<string> Zeilen;
	Zeilen = test.Openfile(argv[1]);
	std::cout << "Datei ist eingelesn" << std::endl;
	std::vector<pair<int, int> > Blocks;
	Blocks = test.FindBlocks(Zeilen);
	std::cout << "Blocks wurden erstellt" << std::endl;

	string tester(argv[1]);

	string file = tester.substr(0, tester.find("."));
	file = file + "_converted_UG";
	std::cout << file << std::endl;

	test.WriteStart(file, Blocks, Zeilen);

	std::cout << "Files wurden geschrieben" << std::endl;
	std::cout << test.Unit_Conv_Value("(mm)") << std::endl;
	std::cout << "neues: " << std::endl;
	std::cout << test.Unit_Conv_All("(mho/cm2)") << std::endl; // /10000
}
else
{
	std::cout << "Give only 1 Argument which should be the Nmodl-File" << std::endl;
}



return 0;
}

