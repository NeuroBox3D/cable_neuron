/*
 * main.cpp
 *
 *  Created on: 30.09.2013
 *      Author: unheiliger
 */

#include "Converter.h"
#include <iostream>



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

	test.WriteStart(argv[2], Blocks, Zeilen);

	std::cout << "Files wurden geschrieben" << std::endl;
	std::cout << test.Unit_Conv("(mm)") << std::endl;
}
else
{
	std::cout << "Bitte 2 Argumente angeben, das erste Argument sollte die NModl File sein, das 2. der Dateiname fuer die Ausgabe" << std::endl;
}



return 0;
}

