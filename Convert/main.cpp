/*
 * main.cpp
 *
 *  Created on: 30.09.2013
 *      Author: unheiliger
 */

#include "Converter.h"
#include <iostream>
#include <stdexcept>



int main(int argn, char* argv[]) {
if (argn == 3)
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
	string forceS(argv[2]);
	string added;
	bool test1;
	if (forceS=="false")
	{
		test1 = false;
		added = "standard";
	}
	if (forceS=="true")
	{
		test1 = true;
		added = "allNernst";
	}

	if (forceS=="false" || forceS=="true")
	{

		string file = tester.substr(0, tester.find("."));
		file = file + "_converted_" + added + "_UG";
		std::cout << file << std::endl;

		test.WriteStart(file, Blocks, Zeilen, test1);

	/*std::cout << "Files wurden geschrieben" << std::endl;
	std::cout << test.Unit_Conv_Value("(mm)") << std::endl;
	std::cout << "neues: " << std::endl;
	std::cout << test.Unit_Conv_All("(pS/um2)") << std::endl; // /10000*/

	} else
	{
		std::cout << "second Argument has to be false/true" << std::endl;
	}

}
else
{
	std::cout << "First Argument should be the Nmodl-File, second argument true or false" << std::endl;
}



return 0;
}

