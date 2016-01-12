/*
 * main.cpp
 *
 *  Created on: 30.09.2013
 *      Author: unheiliger
 */

#include <iostream>
#include <stdexcept>
#include "Converter.h"



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

		std::cout << "Channel file written" << std::endl;

		std::vector<string> needed_files;
		needed_files = test.WriteChannelFile(file, "channels.cpp");

		std::cout << "Channels file written" << std::endl;

		test.WriteInclude_List(needed_files, "includefile.h", "channel_sources");

		std::cout << "Files wurden geschrieben" << std::endl;

		test.WriteInPlugin("channel_sources", "../../CMakeLists.txt");


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
