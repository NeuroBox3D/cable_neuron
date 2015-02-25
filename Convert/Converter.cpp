/*
 * Converter.cpp
 *
 *  Created on: 04.02.2015
 *      Author: unheiliger
 */

#include "Converter.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>



Converter::Converter() {
	// TODO Auto-generated constructor stub

}

Converter::~Converter() {
	// TODO Auto-generated destructor stub
}

double Converter::Unit_Conv(string s)
{
	double erg = 0;
	//needed units /c /m^2 /mV /ms
	// existing units amper: pA, mA, A
	// existing units volt: pV, mV, V
	// mm2 cm2 dm2 m2 km2
	// existing units time: ms, s,
	// working with Siemens!!
	// existing unigs simens: pS, mS, S
	// existing units concentrats: mM,
	// existing units temp: degC
	// (um) = (micron) ????
	// mho/cm2 mA/cm2 S/cm2
	// 2 standing for ²


	// In our model we have mA and mV
	if ((s.find("V")!=s.npos) && (s.find("A")!=s.npos))
	{
		if ((s.find("pV")!=s.npos) || (s.find("pA")!=s.npos))
		{
			if (erg==0)
				erg = 1e-9;
			else erg *= 1e-9;
		}

		if ((s.find("nV")!=s.npos) || (s.find("nA")!=s.npos))
		{
			if (erg==0)
				erg = 1e-6;
			else erg *= 1e-6;
		}

		if ((s.find("V)")!=s.npos) || (s.find("A)")!=s.npos))
		{
			if (erg==0)
				erg = 1e3;
			else erg *= 1e3;
		}
	}

	//special mho is siemens
	// In our model we need Siemens
	if ((s.find("S")!=s.npos))
	{
		if ((s.find("pS")!=s.npos))
		{
			if (erg==0)
				erg = 1e12;
			else erg *= 1e12;
		}

		if ((s.find("nS")!=s.npos))
		{
			if (erg==0)
				erg = 1e9;
			else erg *= 1e9;
		}



		if ((s.find("mS")!=s.npos))
		{
			if (erg==0)
				erg = 1000;
			else erg *= 1000;
		}
	}





	// In our model we have m/m2/m3/m4
	string dims[] = {"", "2", "3", "4"};
	int dimi[] = {1, 100, 10000, 1000000};
	int factor = 0;
	size_t pos;

	for (size_t i = 0; i<4; i++)
	{
		pos = s.find("m"+dims[i]);
		if (pos!=s.npos)
		{
			for (size_t j = 1; j<4; j++)
			{
				if (s.find("m"+dims[j])!=s.npos)
				{
					i = j;
				}
			}
			factor = i+1;

			if (s.find("um"+dims[i])!=s.npos)
			{
				if (erg==0)
					erg = std::pow(10, factor)*100000*dimi[i];
				else erg *= std::pow(10, factor)*100000*dimi[i];
			}


			if (s.find("mm"+dims[i])!=s.npos)
			{
				if (erg==0)
					erg = std::pow(10, factor)*100*dimi[i];
				else erg *= std::pow(10, factor)*100*dimi[i];
			}

			if (s.find("cm")!=s.npos)
			{
				if (erg==0)
					erg = std::pow(10, factor)*10*dimi[i];
				else erg *= std::pow(10, factor)*10*dimi[i];
			}

			if (s.find("dm")!=s.npos)
			{
				if (erg==0)
					erg = std::pow(10, factor)*dimi[i];
				else erg *= std::pow(10, factor)*dimi[i];
			}

			if (s.find("km")!=s.npos)
			{
				if (erg==0)
					erg = std::pow(0.1, factor)*0.01*(1/dimi[i]);
				else erg *= std::pow(0.1, factor)*0.01*(1/dimi[i]);
			}

		}

	}



	// In our model we have ms

	return erg;
}


std::vector<string> Converter::GetProcEqualString(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s)
{
	std::vector<string>erg;
	std::vector<string>PROCEDURE = GetBlock(Pairs, Zeilen, "PROCEDURE");

	size_t beg = 0;
	size_t end = 0;

	for (size_t i=0; i<PROCEDURE.size(); i++)
	{
		string test;
		test = "PROCEDURE " + s;
		//std::cout << "vergleich von " << test << "mit: " << std::endl;
		std::cout << PROCEDURE[i] << std::endl;
		if (PROCEDURE[i].find(test)!=PROCEDURE[i].npos)
		{
			beg = i;
			std::cout << "gleich" << std::endl;
		}
		if ((PROCEDURE[i].find("PROCEDURE")!=PROCEDURE[i].npos) && (i!=beg))
		{
			end = i;
			i = PROCEDURE.size();
		}
	}
	if (end == 0)
		end = PROCEDURE.size()-1;

	std::cout << "begin: " << beg << " ende: " << end << std::endl;

	//if (beg != 0)
	//{
		std::cout << "GetProcEqual does anything" << std::endl;
		// Rekursion könnte hier mehrfachaufrufe von prozeduren möglich machen
		for (size_t i=beg; i<end; i++)
		{
			// deletes all not needed parts
			PROCEDURE[i] = writing_starts(PROCEDURE[i]);
			if ((PROCEDURE[i]!="") && (PROCEDURE[i]!="\t"))
				erg.push_back(PROCEDURE[i]);
		}
	//}

	return erg;
}

bool Converter::Only_Read(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s)
{
	bool erg = false;

	std::vector<string> NEURON = GetBlock(Pairs, Zeilen, "NEURON");
	for (size_t i=0; i<NEURON.size(); i++)
	{
		if (NEURON[i].find("USEION " + s)!=NEURON[i].npos)
		{
			if (NEURON[i].find("WRITE")==NEURON[i].npos)
			{
				erg = true;
			}
		}
	}

	return erg;
}

string Converter::Write_Only_Read(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s)
{
	string erg = "double " + s + " = aa";
	size_t pos;
	size_t end;

	std::vector<string> NEURON = GetBlock(Pairs, Zeilen, "NEURON");
	for (size_t i=0; i<NEURON.size(); i++)
	{
		pos = NEURON[i].find("USEION");
		if (NEURON[i].find(s)!=NEURON[i].npos)
		{
			end = NEURON[i].find("READ");
			if (end!=NEURON[i].npos)
			{
				//aacaGate[*iter]
				std::cout << "beg: " << pos+7 << " dur: " << end-pos+7 << "end: "<< end << std::endl;
				erg = erg + NEURON[i].substr(pos+7, (end-(pos+7))-1)+"Gate[*iter]";
			}
		}
	}
	return (erg + "; \n");
}


bool Converter::begG(string s)
{
	//std::cout << "fehler in beg" << std::endl;
	bool erg = false;
	const string KeyWord[] = {"TABLE", "DEPEND", "FROM", "TO", "HOC"};
	for ( size_t i=0 ; i<5; i++ )
	{
		if (s.find(KeyWord[i])!=s.npos)
			erg = true;
		//std::cout << "anzahl i's: " << i << std::endl;
	}
	//std::cout << "vor rueckgabe" << std::endl;
	return erg;
}


// writing Nernst Equatation for Ions with write
std::vector<string> Converter::equali(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen)
{
	std::vector<string> out;
	out.push_back("const number helpV = (m_R*m_T)/m_F;");

	size_t Ion, IonRead, IonRend;
	string IonS;
	string PartR, PartW;

	size_t R_beg, R_end;


	std::vector<string> ListIonRead;
	std::vector<string> ListIon;
	std::vector<bool> version;

	std::vector<string> NEURON = GetBlock(Pairs, Zeilen, "NEURON");

	for (size_t i=0; i<NEURON.size(); i++)
	{
		Ion = NEURON[i].find("USEION");
		if (Ion!=NEURON[i].npos)
		{
			// needs two different version one with "," one without ","
			// first verion with ","
			if (NEURON[i].find(",")!=NEURON[i].npos)
			{
				version.push_back(true);
				out.push_back("0");
				// 7 for USEION
				PartR = NEURON[i].substr(7, NEURON[i].find(",")-7);
				PartW = NEURON[i].substr(NEURON[i].find(",")+1, NEURON[i].npos-(NEURON[i].find(",")+1));
				// PartR means cai = ca
				// PartW means: ica = cao means ca_out ica changes conc out
				// if 0 needed extra ion
				R_beg = PartR.find("READ")+5;
				R_end = PartR.find("READ")-1;
				out.push_back(PartR.substr(0, R_end));
				out.push_back("number " + PartR.substr(R_beg, PartR.npos) + " = " + PartR.substr(0, R_end) + ";");

			} else // version without ","
			{
				version.push_back(false);
				out.push_back("1");
				IonRead = NEURON[i].find("READ", Ion+7);
				IonRend = NEURON[i].find(" ", IonRead+5);
				std::cout << IonRend << ", " << Ion << std::endl;
				IonS = NEURON[i].substr(Ion+7, IonRead-1-(Ion+7));
				std::cout << "Subion: "<< IonS << std::endl;
				ListIonRead.push_back(NEURON[i].substr(IonRead+5 , IonRend-(IonRead+5)));
				// Only if there is a WRITE in
				if (NEURON[i].find("WRITE")!=NEURON[i].npos)
						ListIon.push_back(IonS);
			}
		}
	}

		for (size_t j=0; j<ListIon.size(); j++)
		{
			if (version[j]==false)
			{
			std::cout << "read ion: " + ListIonRead[j] + "List ion: " + ListIon[j] << std::endl;
			out.push_back("number " + ListIonRead[j] + " = helpV*(log(" + ListIon[j] + "_out/" + ListIon[j] + "));");
			}
		}


		for (size_t j=0; j<ListIon.size(); j++)
		{
			if (version[j]==true)
			{
			std::cout << "read ion: " + ListIonRead[j] + "List ion: " + ListIon[j] << std::endl;
			out.push_back("number " + ListIonRead[j] + " = helpV*(log(" + ListIon[j] + "_out/" + ListIon[j] + "));");
			}
		}



	return out;
}


string Converter::func_head(string s)
{
	string fname;
	size_t startf, endf;
	std::vector<string> vars;
	startf = s.find(" ");
	endf = s.find("(");
	fname = s.substr(startf+1, endf-startf-1);
	std::cout << "fname: " << fname << std::endl;
	return fname;
}

std::vector<vector<string> > Converter::write_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen)
{
	std::vector<vector <string> > Proc_blocks;
	std::vector<string> vars;
	size_t unit_beg;

	std::vector<string>  PROCEDURE = GetBlock(Pairs, Zeilen, "PROCEDURE");


	for (size_t j = 0; j<PROCEDURE.size(); j++)
	{
		if (PROCEDURE[j].find("PROCEDURE") != PROCEDURE[j].npos)
		{
			vars = build_proc_vars(PROCEDURE[j]);
			for (size_t i = 0; i<vars.size(); i++)
			{
				//TODO write Unit -> changes

				// deletes Unit
				unit_beg = vars[i].find("(");
				vars[i] = vars[i].substr(0, unit_beg);
			}

			Proc_blocks.push_back(vars);
		}
	}


	return Proc_blocks;
}


std::vector<vector<string> > Converter::write_derivative_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen)
{
	std::vector<vector <string> > Proc_blocks;
	std::vector<string> vars;
	size_t unit_beg;

	std::vector<string>  DERIVATIVE = GetBlock(Pairs, Zeilen, "DERIVATIVE");


	for (size_t j = 0; j<DERIVATIVE.size(); j++)
	{
		if (DERIVATIVE[j].find("DERIVATIVE") != DERIVATIVE[j].npos)
		{
			vars = build_proc_vars(DERIVATIVE[j+1]);
			for (size_t i = 0; i<vars.size(); i++)
			{
				//TODO write Unit -> changes

				// deletes Unit
				unit_beg = vars[i].find("(");
				vars[i] = vars[i].substr(0, unit_beg);
			}

			Proc_blocks.push_back(vars);
		}
	}

	return Proc_blocks;
}




std::vector<string> Converter::get_local_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name)
{
	std::vector<string> locals;
	std::vector<string> vars;
	size_t beg;

	std::vector<string>  PROCEDURE = GetBlock(Pairs, Zeilen, "PROCEDURE");

	for (size_t i = 0; i<PROCEDURE.size(); i++)
	{
		beg = PROCEDURE[i].find(name);
		if (beg!=PROCEDURE[i].npos)
		{
			for (size_t j = i; j<PROCEDURE.size(); j++)
				{
				if (PROCEDURE[j].find("LOCAL")!=PROCEDURE[j].npos)
				{
					PROCEDURE[j].replace(PROCEDURE[j].find("LOCAL"), 5, "");
					locals.push_back(PROCEDURE[j]);
				}
			}
		}
	}
	return locals;
}


string Converter::writing_starts(string s)
{
	string test[] = {"UNITSOFF" ,"TABLE", "DEPEND", "FROM", "TO", "HOC, LOCAL"};

	for (size_t i = 0; i<6; i++)
	{
		if (s.find(test[i])!=s.npos)
			s = "";
	}
	return s;
}




std::vector<string> Converter::writer_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name)
{
	std::vector<string> fileInfo;
	std::vector<string> vars;
	std::vector<string> fileInfo_norm;
	std::vector<string> fileInfo_rec;
	size_t beg;
	size_t recs = 0;

	std::vector<vector <string> > Proc_list;

	std::vector<string>  PROCEDURE = GetBlock(Pairs, Zeilen, "PROCEDURE");

	// list with all proc funcs

	// only [j][0] have functionnames
	Proc_list = write_proc_block(Pairs, Zeilen);

	for (size_t i = 0; i<PROCEDURE.size(); i++)
	{
		beg = PROCEDURE[i].find("PROCEDURE " + name);
		if (beg!=PROCEDURE[i].npos)
		{
			for (size_t j = i; j<PROCEDURE.size(); j++)
				{
					PROCEDURE[j] = writing_starts(PROCEDURE[j]);
				}

			for (size_t k = i+1; k<PROCEDURE.size(); k++)
				{
					if (PROCEDURE[k]=="}")
						k = PROCEDURE.size()-1;
					bool rec = false;
					//if again PROCEDURE action than start new this one here
					for (size_t l=0; l<Proc_list.size(); l++)
					{
						if (PROCEDURE[k].find(Proc_list[l][0])!=PROCEDURE[k].npos)
						{
							std::cout << "zeug mal wieder: " << std::endl;
							std::cout << PROCEDURE[k].substr(PROCEDURE[k].find("(")+1, PROCEDURE[k].find(")")-PROCEDURE[k].find("(")-1) << std::endl;
							std::cout << Proc_list[l][1] << std::endl;
							fileInfo.push_back(Proc_list[l][1] + " = " + PROCEDURE[k].substr(PROCEDURE[k].find("(")+1, PROCEDURE[k].find(")")-PROCEDURE[k].find("(")-1));
							fileInfo_rec = writer_proc_block(Pairs, Zeilen, Proc_list[l][0]);
							//fileInfo = writer_proc_block(Pairs, Zeilen, Proc_list[l][0]);
							rec = true;
							recs += 1;

						}
						if (rec==true)
						{
							// adding all out of recursion in fileInfo output
							for (size_t l=0; l<fileInfo_rec.size(); l++)
							{
								fileInfo.push_back(fileInfo_rec[l]);
							}
						}
					}




					// writing everything that is not recursiv
					if (rec==false)
						//fileInfo.push_back(PROCEDURE[k]);
						fileInfo_norm.push_back(PROCEDURE[k]);

				}
		}
	}



	// adding all of not recursiv
	for (size_t l=0; l<fileInfo_norm.size(); l++)
	{
		fileInfo.push_back(fileInfo_norm[l]);
	}

	return fileInfo;
}



std::vector<string> Converter::build_proc_vars(string s)
{
	string fname, var;
	size_t startf, endf, endvars, buf, buf1;
	std::vector<string> vars;
	startf = s.find(" ");
	endf = s.find("(");
	fname = s.substr(startf+1, endf-startf-1);
	//std::cout << fname << std::endl;
	endvars = s.find_last_of(")");
	//std::cout << "beg " << endf << " end " << endvars << std::endl;
	var = s.substr(endf+1, endvars-endf-1);
	//std::cout << var << std::endl;
	buf = 0;
	// writs all vars out
	vars.push_back(fname.substr(count_beg(fname), fname.npos-count_beg(fname)));
	while (var.find(",", buf)!=var.npos)
	{
		buf1 = var.find(",", buf);
		vars.push_back(var.substr(buf, buf1-buf));
		buf = buf1+1;
	}
	vars.push_back(var.substr(buf, var.npos));
	return vars;
}

string Converter::build_func_head(string s)
{
	string out = "double " + m_class_name + "::";
	string fname, var;
	size_t startf, endf, endvars, buf, buf1;
	std::vector<string> vars;
	startf = s.find(" ");
	endf = s.find("(");
	fname = s.substr(startf+1, endf-startf-1);
	//std::cout << fname << std::endl;
	endvars = s.find_last_of(")");
	//std::cout << "beg " << endf << " end " << endvars << std::endl;
	var = s.substr(endf+1, endvars-endf);
	// if there is also a unit del it
	std::cout << "var: "<< var << std::endl;
	buf = 0;
	string buf2;
	// writs all vars out
	while (var.find(",", buf)!=var.npos)
	{
		buf1 = var.find(",", buf);
		buf2 = var.substr(buf, buf1-buf);
		if (buf2.find("(")!=var.npos)
		{
			buf2 = var.substr(buf, buf2.find("("));
			std::cout << "buf2: " << buf2 << std::endl;
		}
		vars.push_back(buf2);
		buf = buf1+1;
	}
	std::cout << "buf:  " << buf <<"bis: " << var.find("(", buf)-buf << std::endl;
	buf2 = var.substr(buf, var.find("(", buf)-buf);
	vars.push_back(buf2);

	out += fname + "(";

	for (size_t i=0; i<vars.size()-1; i++)
	{
		out += "double " + vars[i] + ", ";
	}
	out += "double " + vars[vars.size()-1] + ")";
	if (out.find("))")!=out.npos)
		out.replace(out.find("))"), 2, ")");
    return out;
}


size_t Converter::count_beg(string beg)
{
	size_t begi = 0;
	if (beg.find("\t")!=beg.npos)
	{
		begi = 1;
	}

	while (beg.find(" ", begi, 1)!=beg.npos)
	{
		begi += 1;
	}
	return begi;
}

size_t Converter::find_beg(string beg)
{
	size_t begi;
	begi = beg.find("{");
	if (begi!=beg.npos)
	{
		return (begi+1);
	}

	begi = beg.find(" ");
	if (begi!=beg.npos)
	{
		return (begi);
	}
	return begi;

}

size_t Converter::find_begg(string beg)
{
	size_t begi1, begi2;
	// deletes " " at the begining afterwards

	begi1 = beg.find(" ");
	if (begi1!=beg.npos)
	{
		begi2 = beg.find(" ", begi1+1);
		while (begi2==begi1+1)
		{
			begi1 +=1;
			begi2 = beg.find(" ", begi1+1);
		}
	}
	return begi1+1;

}


size_t Converter::number_(size_t pos, string s)
{
	size_t counter = 0;
	size_t begin, breakk;

	breakk = s.find("\n");
	while (breakk!=s.npos)
	{
		s.replace(breakk, 1, "");
		breakk = s.find("\n");
	}


	begin = s.find(" ", pos);
	while (begin !=s.npos)
	{
		counter += 1;
		begin = s.find(" ", begin+1);
	}
	if (counter > 1)
		counter +=1;
	return counter;
}


size_t Converter::pos_letter(string s)
{
	vector<size_t> pos;
	size_t erg;

	size_t beg;
	if (s.find("STATE")!=s.npos)
	{
		beg = s.find("STATE")+6;
	} else
	{
		beg = 0;
	}

	//std::cout << "pos_letter does anything" << std::endl;
	// little letters
	for ( char i( 97 ); i < 122; ++i )
	{
		//97-122
		if (s.find(i, beg)!=s.npos)
			pos.push_back(s.find(i));
	}
	// big letters
	for ( char i( 65 ); i < 90; ++i )
	{
		//97-122
		if (s.find(i, beg)!=s.npos)
			pos.push_back(s.find(i));
	}

	if (pos.size()!= 0)
		erg = *std::min_element(pos.begin(), pos.end());
	else return 1000;

	return erg;
}

bool Converter::pos_letterb(string s)
{
	vector<size_t> pos;
	std::cout << "pos_letter does anything " << s << std::endl;
	size_t beg;
	if (s.find("STATE")!=s.npos)
	{
		beg = s.find("STATE")+6;
	} else
	{
		beg = 0;
	}


	// little letters
	for ( char i( 97 ); i < 122; ++i )
	{
		//97-122
		if (s.find(i, beg)!=s.npos)
			return false;
	}

	// big letters
	for ( char i( 65 ); i < 90; ++i )
	{
		//97-122
		if (s.find(i, beg)!=s.npos)
			return false;
	}

	return true;
}


std::vector<string> Converter::Openfile(string filename)
{

	std::vector<string> zeilen;
	//try {

	const char* filenamechar = filename.c_str();
	ifstream ModFile(filenamechar, std::ios::out);

	string zeile;
		while (!ModFile.eof())          // Solange noch Daten vorliegen
	{
		//std::cout << "while" << std::endl;
	    getline(ModFile, zeile);        // Lese eine Zeile
	    /*if (zeile.find(":")!=zeile.npos)
	    {
	    	std::cout << "replace happens" << std::endl;
	    	zeile.replace(zeile.find(":"), 1,";//");
	    }// replace : with // as comment*/
	    zeilen.push_back(zeile);    // Zeige sie auf dem Bildschirm
	}
	ModFile.close();    //}
	  //catch (std::ifstream::failure e) {
	    //std::cerr << "Exception opening/reading/closing file\n";
	  //}

	std::cout << "einlesen erfolgreich" << std::endl;
	return zeilen;
}

std::vector<std::pair<int, int> > Converter::FindBlocks(std::vector<string> Zeilen)
{
	std::vector<std::pair<int, int> > Blocks;
	size_t start, ende;
	start = 0;
	ende = 0;
	pair<int, int> block;
	size_t blockH;
	bool works = false;


	for (size_t i=0; i<Zeilen.size(); i++)
	{

		start = Zeilen[i].find("{");
		ende = Zeilen[i].find("}");
		if (start!=Zeilen[i].npos)
		{
			if (works==false)
			{
				blockH = i;
				works = true;
			}
		}
		if (ende!=Zeilen[i].npos)
		{
			Blocks.push_back(std::make_pair(blockH, i+1));
			works = false;
		}

	}
	return Blocks;

}


std::vector<string> Converter::GetBlock(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name)
{
	std::vector<string> BlockName;
	bool Verba = false;

	for (size_t i=0; i<Pairs.size(); i++)
	{
		//std::cout << "first" << Pairs.size() << std::endl;
		// first find of right name and second find that not outsurced
		if ((Zeilen[Pairs[i].first].find(name)!=Zeilen[Pairs[i].first].npos) && (Zeilen[Pairs[i].first].find(":"+name)==Zeilen[Pairs[i].first].npos))
		{
			//std::cout << Zeilen[Pairs[i].first] << std::endl;
			//überprüft wenn weitere klammern da sind, dass nur die letzte verwendet wird;

			while (Pairs[i].first == Pairs[i+1].first)
			{
				//std::cout << "in whhile by i: " << i << std::endl;
				i = i+1;
			}
			//überspringt doppelte Pairs
			//std:: cout << counter << " , " << Pairs[i+counter].second <<std::endl;

			size_t breaks;

			// schreibt zeilen raus
			for (int j=Pairs[i].first; j<Pairs[i].second; j++)
			{
				if (Zeilen[j].find("VERBATIM")!=Zeilen[j].npos)
					Verba = true;

				// removes verbatim blocks
				if (Verba==false)
				{
				//std::cout << counter << std::endl;
					//deletes also every \n in blocks
					breaks = Zeilen[j].find("\n");
					while (breaks!=Zeilen[j].npos)
					{
						Zeilen[j].replace(breaks, 1, "");
						breaks = Zeilen[j].find("\n");
					}

					//deletes also every \r in blocks
					breaks = Zeilen[j].find("\r");
					while (breaks!=Zeilen[j].npos)
					{
						Zeilen[j].replace(breaks, 1, "");
						breaks = Zeilen[j].find("\r");
					}
					BlockName.push_back(Zeilen[j]);
				}
				if (Zeilen[j].find("ENDVERBATIM")!=Zeilen[j].npos)
					Verba=false;

			}

		}
	}
	//std::cout << "before block" << std::endl;
	//std::cout << BlockName.size() << std::endl;
	return BlockName;
}

string Converter::Remove_all_com(string erg)
{
	size_t HFile_del;

	HFile_del = erg.find("\n");
	while (HFile_del!=erg.npos)
	{
		erg.replace(HFile_del, 1, "");
		HFile_del = erg.find("\n");
	}

		HFile_del = erg.find("\r");
	while (HFile_del!=erg.npos)
	{
		erg.replace(HFile_del, 1, "");
		HFile_del = erg.find("\r");
	}

	HFile_del = erg.find("\t");
	while (HFile_del!=erg.npos)
	{
		erg.replace(HFile_del, 1, "");
		HFile_del = erg.find("\t");
	}

	while (erg.substr(0 ,1)==" ")
	{
		erg.replace(0, 1, "");
	}

	return erg;


}


string Converter::Remove_all(string erg)
{
	size_t HFile_del;

	HFile_del = erg.find(" ");
	while (HFile_del!=erg.npos)
	{
		erg.replace(HFile_del, 1, "");
		HFile_del = erg.find(" ");
	}
	HFile_del = erg.find("\n");
	while (HFile_del!=erg.npos)
	{
		erg.replace(HFile_del, 1, "");
		HFile_del = erg.find("\n");
	}

		HFile_del = erg.find("\r");
	while (HFile_del!=erg.npos)
	{
		erg.replace(HFile_del, 1, "");
		HFile_del = erg.find("\r");
	}

	HFile_del = erg.find("\t");
	while (HFile_del!=erg.npos)
	{
		erg.replace(HFile_del, 1, "");
		HFile_del = erg.find("\t");
	}


	return erg;


}


vector<string> Converter::Remove_all(vector<string> erg)
{
	  size_t HFile_del;
	  //delete all not needed symbols out of HFile_added_Vars
	  for (size_t i=0; i<erg.size(); i++)
	  {
		  // dels all " "
		  HFile_del = erg[i].find(" ");
		  while (HFile_del!=erg[i].npos)
		  {
			  erg[i].replace(HFile_del, 1, "");
			  HFile_del = erg[i].find(" ");
		  }

		  HFile_del = erg[i].find("\n");
		  while (HFile_del!=erg[i].npos)
		  {
			  erg[i].replace(HFile_del, 1, "");
			  HFile_del = erg[i].find("\n");
		  }

		  HFile_del = erg[i].find("\r");
		  while (HFile_del!=erg[i].npos)
		  {
			  erg[i].replace(HFile_del, 1, "");
			  HFile_del = erg[i].find("\r");
		  }

		  HFile_del = erg[i].find("\t");
		  while (HFile_del!=erg[i].npos)
		  {
			  erg[i].replace(HFile_del, 1, "");
			  HFile_del = erg[i].find("\t");
		  }


		  std::cout<< "Hfile_vars: " << erg[i] << std::endl;
	  }
	  return erg;
}


void Converter::WriteStart(string filename, std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen)
{
	  // defines filenames as const char*
	  string fnamecpp = filename + ".cpp";
	  string fnameh = filename + ".h";
	  const char* filenamecpp = fnamecpp.c_str();
	  const char* filenameh = fnameh.c_str();


	  vector<string> HFile_added_Vars;



	  m_class_name = filename;


	  // open cpp file + writing init head
	  ofstream mycppfile;
	  mycppfile.open (filenamecpp);
	  std::cout << "write begins" << std::endl;
	  mycppfile << "#include \" " + filename +".h\"	\n";
	  mycppfile << "#include \" lib_grid/lg_base.h\" \n";
	  mycppfile << "#include \" lib_disc/spatial_disc/elem_disc/elem_disc_interface.h\" \n";
	  mycppfile << "#include \" lib_disc/function_spaces/grid_function.h\" \n";
	  mycppfile << "#include \" lib_disc/function_spaces/local_transfer_interface.h\" \n";
	  mycppfile << "#include \" <cmath> \" \n";

	  mycppfile << "namespace ug { \n \n \n";

	  std::vector<string> mod_funcs_names;
	  // all function have to be read in first so later need more attention more

	  std::vector<string> FUNCTION = GetBlock(Pairs, Zeilen, "FUNCTION");

	  string without[] = {"if", "else"};
	  size_t funcS, funcE;
	  bool second_time;

	  string func, funchead;
	  if (FUNCTION.size() > 0)
	  {
		  for (size_t j=0; j<FUNCTION.size(); j++)
		  {
			  second_time = false;
			  if (FUNCTION[j].find("FUNCTION")!=FUNCTION[j].npos)
			  {
				  func = build_func_head(FUNCTION[j]);
				  funchead = func_head(FUNCTION[j]);

				  mycppfile << func +" \n";
				  mycppfile << "{ \n";

				  for (size_t i = j+1; i<FUNCTION.size(); i++)
				  {
					  if (FUNCTION[i].find("FUNCTION")!=FUNCTION[i].npos)
					  {

						  //if (i!=j)
							  i = FUNCTION.size()-1;
							  second_time = true;

					  }

					  funcS = FUNCTION[i].find(funchead);
					  if (funcS!=FUNCTION[i].npos)
					  {
						  funcE = FUNCTION[i].find(" ", funcS);
						  FUNCTION[i].replace(funcS, funcE-funcS + 2, "return ");
					  }

					  if (FUNCTION[i].find("FUNCTION")==FUNCTION[i].npos)
					  {
						  if (FUNCTION[i].find("LOCAL")!=FUNCTION[i].npos)
						  {
							  FUNCTION[i].replace(FUNCTION[i].find("LOCAL"), 6, "double ");
						  }
						  if (second_time==false)
						  {
							  if ((FUNCTION[i].find(without[0])==FUNCTION[i].npos) && (FUNCTION[i].find(without[1])==FUNCTION[i].npos)
								 && (string(FUNCTION[i] + "; \n").find("};")==string(FUNCTION[i] + "; \n").npos) && (FUNCTION[i]!=""))
							  {
								  mycppfile << FUNCTION[i] + "; \n";
							 	  //std::cout << FUNCTION[i];
							  }
							  else
							  {
							 	  //std::cout << FUNCTION[i];
							 	  mycppfile << FUNCTION[i] + "\n";
							  }

						  }
					  }
				  }
				  //mycppfile << "} \n \n";

			  }
		  }

	  	  mycppfile << "\n \n";
	  }

	  // Proceduren wie Funktionen
	  std::cout << "All Functions added" << std::endl;



	  mycppfile.close();

	  ofstream myhfile;
	  myhfile.open (filenameh);
	  myhfile << "#ifndef " + filename + "_H_\n";
	  myhfile << "#define " + filename + "_H_\n";

	  myhfile << "#include \"channel_interface.h\" \n";

	  myhfile << "#include \"lib_grid/lg_base.h\" \n";
	  myhfile << "#include \"lib_grid/grid/grid_base_objects.h\" \n";
	  myhfile << "\n";

	  myhfile << "#include \"lib_disc/spatial_disc/elem_disc/elem_disc_interface.h\" \n";
	  myhfile << "#include \"lib_disc/spatial_disc/disc_util/fv1_geom.h\" \n";
	  myhfile << "#include \"lib_disc/spatial_disc/disc_util/hfv1_geom.h\" \n";
	  myhfile << "#include \"lib_disc/spatial_disc/disc_util/geom_provider.h\" \n";
	  myhfile << "#include \"lib_disc/function_spaces/grid_function.h\" \n";
	  myhfile << "#include \"lib_disc/function_spaces/local_transfer_interface.h\" \n";
	  myhfile << "#include \"lib_disc/common/local_algebra.h\" \n";
	  myhfile << "#include \"lib_disc/function_spaces/grid_function.h\" \n";
	  myhfile << "#include \"lib_disc/function_spaces/interpolate.h\" \n";
	  myhfile << "\n";

	  myhfile << "#include \"bridge/bridge.h\" \n";
	  myhfile << "#include \"bridge/util.h\" \n";
	  myhfile << "#include \"bridge/util_domain_algebra_dependent.h\" \n";
	  myhfile << "#include \"bridge/util_domain_dependent.h\" \n";
	  myhfile << "\n";

	  myhfile << "#include \"common/util/smart_pointer.h\" \n";
	  myhfile << "#include \"common/util/vector_util.h\" \n";
	  myhfile << "\n";

	  myhfile << "#include <vector> \n";
	  myhfile << "#include <stdio.h> \n";
	  myhfile << "#include \"bindings/lua/lua_user_data.h\" \n";

	  myhfile << "namespace ug { \n";
	  myhfile << "template <typename TDomain, typename TAlgebra> \n";
	  myhfile << "class " + filename + "\n";
	  myhfile << "    : public IChannel<TDomain, TAlgebra> \n";
	  myhfile << "{ \n";
	  myhfile << "    public: \n \n \n";

	  // finde die spezifischen Parameter
	  // out of neuron only GLOBALS needed
	  std::vector<string> NEURON = GetBlock(Pairs, Zeilen, "NEURON");
	  size_t Global;
	  std::vector<string> Global_vars;
	  // Globals out of Neuron
	  for (size_t i = 0; i<NEURON.size(); i++)
	  {
		  Global = NEURON[i].find("GLOBAL ");
		  if (Global!=NEURON[i].npos)
		  {
			  Global_vars = build_proc_vars(NEURON[i].substr(Global+7, NEURON[i].npos));
		  }
	  }

	  // writting all globals
	  if (Global_vars.size() > 0)
	  {
		  for (size_t i = 0; i<Global_vars.size(); i++)
		  {
			  HFile_added_Vars.push_back(Global_vars[i]);
			  myhfile << "double " + Global_vars[i] + "; \n";
		  }

	  }
	  myhfile << "\n \n \n" ;

	  std::cout << "Globals added!" << std::endl;

	  // getting out of Params hard coded Values
	  std::vector<string> PARAMETER = GetBlock(Pairs, Zeilen, "PARAMETER");
	  size_t end;
	  string var;

	  if (PARAMETER.size() > 0)
	  {
		  for (size_t i=0; i<PARAMETER.size(); i++)
		  {
			  if (PARAMETER[i].find("=")!=PARAMETER[i].npos)
			  {
				  end = PARAMETER[i].find("(");
				  var = PARAMETER[i].substr(0, end);
				  myhfile << "const static double " + var + "; \n";
				  HFile_added_Vars.push_back(var.substr(0, var.find("=")));
			  }
		  }
	  }
	  myhfile << "std::vector<number> m_diff; \n";
	  myhfile << "\n";


	  std::cout << "All Parameters set" << std::endl;

	  size_t Ion, IonEnd;
	  string IonS;
	  std::vector<string> ListIons;

	  // adding ions from neuron block with outer conz
	  for (size_t i=0; i<NEURON.size(); i++)
	  {
		  // adding ions
		  Ion = NEURON[i].find("USEION");
		  if (Ion!=NEURON[i].npos)
		  {

			  IonEnd = NEURON[i].find(" ", Ion+7);
			  std::cout << IonEnd << ", " << Ion << std::endl;
			  IonS = NEURON[i].substr(Ion+7, (IonEnd-(Ion+7)));
			  std::cout << "Subion: "<<IonS << std::endl;
			  ListIons.push_back(IonS);
			  //writes var for outer concentrations
			  myhfile << "number " + IonS + "_out; \n";
			  HFile_added_Vars.push_back(IonS + "_out");
		  }

	  }

	  // delete all unneeded style-formats
	  Remove_all(HFile_added_Vars);


	  myhfile << "\n"; myhfile << "\n";

	  myhfile << "/// constructor \n \n";
	  myhfile << "ChannelHHNernst(const char* functions, const char* subsets) \n";
	  myhfile << ": IChannel<TDomain, TAlgebra>(functions, subsets), \n";
	  myhfile << "m_bNonRegularGrid(false), m_R(8314), m_T(279.45), m_F(96485.0) \n";
	  myhfile << "{ \n";
	  myhfile << "register_all_funcs(m_bNonRegularGrid); \n";
	  myhfile << "}; \n \n";

	  myhfile << "/// destructor \n \n";
	  myhfile << "virtual ~ChannelHHNernst() {}; \n";
	  myhfile << "// inherited from IChannel \n \n";
	  myhfile << "virtual const std::vector<number> get_diff(); \n";
	  myhfile << "virtual void set_diff(const std::vector<number>& diff); \n";
	  myhfile << "virtual void init(number time, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct); \n";
	  myhfile << "virtual void update_gating(number newTime, SmartPtr<GridFunction<TDomain, TAlgebra> > sgGridFct); \n";
	  myhfile << "virtual void ionic_current(Vertex* v, std::vector<number>& outCurrentValues); \n";
	  myhfile << "/// assembles the local right hand side \n \n";
	  myhfile << "template<typename TElem, typename TFVGeom> \n";
	  myhfile << "void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]); \n";
	  myhfile << "\n \n";
	  myhfile << "protected: \n";
	  myhfile << "bool m_bNonRegularGrid; \n \n";
	  myhfile << "void register_all_funcs(bool bHang); \n";
	  myhfile << "template <typename TElem, typename TFVGeom> \n";
	  myhfile << "void register_func(); \n \n \n";
	  myhfile << "private: \n \n";
	  // Neuron-lines with use ion

	  std::vector<string> STATE = GetBlock(Pairs, Zeilen, "STATE");





	  mycppfile.open (filenamecpp, std::ios::app);
	  mycppfile << " // Init Method for using gatings \n";
	  mycppfile << "template<typename TDomain, typename TAlgebra> \n";
	  mycppfile << "void " + filename + "<TDomain, TAlgebra>::init(number time, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) \n";
	  mycppfile << "{ \n";
	  mycppfile << "// attach attachments \n \n";



	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  IonS = ListIons[i];
		  myhfile << "ADouble " + IonS + "Gate; \n";
		  myhfile << "Grid::AttachmentAccessor<Vertex, ADouble> aa" + IonS + "Gate; \n";

		  mycppfile << "if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->" + IonS + "Gate)) \n";
		  mycppfile << "UG_THROW(\"Attachment necessary (" + IonS + "Gate) for Hodgkin and Huxley channel dynamics \"\n";
		  mycppfile << "\"could not be made, since it already exists.\"); \n";
		  mycppfile << "spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->" + IonS + "Gate); \n";
		  mycppfile << "this->aa" + IonS + "Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->" + IonS + "Gate); \n \n";
	  }
////////////////////////////////////////////////////////////////////////////////////
//// Saving all States in one Var also writting first states into cpp and h file
/////////////////////////////////////////////////////////////////////////////////////
	  /*size_t Ns_Current, Ns_CurrentEnd;
	  string Ns_CurrentS;*/

	  string addState;

	  // var in which all states and Non spec currents are saved
	  std::vector<string> State_vars;
	  size_t state, stateend;
	  size_t varSt = 0;

	  if (STATE.size()>0)
	  {
		  for (size_t i=0; i<STATE.size(); i++)
		  {
			  state = pos_letter(STATE[i]);
			  if (state!=1000)
				  if (STATE[i].find("(")==STATE[i].npos)
					  varSt = number_(state, STATE[i]);
			  std::cout << "after doppelt if" << std::endl;
			  for (size_t j=0; j<varSt; j++)
			  {
				  stateend = STATE[i].find(" ", state);
				  std::cout << "ab: " << state << " anzahl: " << stateend-state << std::endl;
				  if (STATE[i].find("}")!=STATE[i].npos)
		  			  {
					  if (stateend-state!=STATE[i].npos-state)
		  			  {
						  addState = STATE[i].substr(state, stateend-state);
		  			  }
					  else
					  {
						  std::cout << "not added!!" << std::endl;
						  addState = "";
					  }
		  			  }
				  	  else
		  			  {
		  				addState = STATE[i].substr(state, stateend-state);
		  			  }

		  			  //writting in hfile
		  			  if ((addState!="}") && (addState!="") && (addState!="\n") &&(addState!="\n}") && (addState!="}\n"))
		  			  {
		  				  State_vars.push_back(addState);
		  				  myhfile << "ADouble " + addState + "Gate; \n";
		  				  myhfile << "Grid::AttachmentAccessor<Vertex, ADouble> aa" + addState + "Gate; \n";

		  				  //writting of attachment inits for cpp file
		  				  mycppfile << "if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->" + addState + "Gate)) \n";
		  				  mycppfile << "UG_THROW(\"Attachment necessary (" + addState + "Gate) for Hodgkin and Huxley channel dynamics \"\n";
		  				  mycppfile << "\"could not be made, since it already exists.\"); \n";
		  				  mycppfile << "spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->" + addState + "Gate); \n";
		  				  mycppfile << "this->aa" + addState + "Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->" + addState + "Gate); \n \n";


		  			  }
		  			state = stateend+1;
		  		  }
		  		  //////////////////

			  varSt = 0;

		  }

	  }

	  // writting VM Parts for gatting init in cpp than close and finish h-file
	  mycppfile << "if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->vGate)) \n";
	  mycppfile << "UG_THROW(\"Attachment necessary (vGate) for Hodgkin and Huxley channel dynamics \"\n";
	  mycppfile << "\"could not be made, since it already exists.\"); \n";
	  mycppfile << "spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->vGate); \n";
	  mycppfile << "this->aavGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->vGate); \n \n";

	  // for every IonS (UseIon)+ v DoFIndex needed
	  mycppfile << "// creates multiindeces \n";
	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "std::vector<DoFIndex> multInd" + ListIons[i] +"; \n";
	  }
	  mycppfile << "std::vector<DoFIndex> multIndv; \n \n";


	  // write begin of interations over elements
	  mycppfile << "typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;  \n";
	  mycppfile << "SubsetGroup ssGrp;  \n";
	  mycppfile << "try { ssGrp = SubsetGroup(spGridFct->domain()->subset_handler(), this->m_vSubset);}  \n";
	  mycppfile << "UG_CATCH_THROW(\"Subset group creation failed.\"); \n  \n";
	  mycppfile << "for (std::size_t si = 0; si < ssGrp.size(); si++)  \n";
	  mycppfile << "{  \n";
	  mycppfile << "itType iterBegin = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);  \n";
	  mycppfile << "itType iterEnd = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); \n  \n";
	  mycppfile << "for (itType iter = iterBegin; iter != iterEnd; ++iter)  \n";
	  mycppfile << "{  \n \n";
	  mycppfile << "size_t vm = spGridFct->fct_id_by_name(\"v\");  \n";
	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "size_t " + ListIons[i] +" = spGridFct->fct_id_by_name(\"" + ListIons[i] + "\"); \n";
	  }
	  mycppfile << "\n \n";
	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "spGridFct->dof_indices(*iter, " + ListIons[i] + ", multInd" + ListIons[i] + "); \n";
	  }
	  mycppfile << "spGridFct->dof_indices(*iter, vm, multIndv);  \n \n ";
	  mycppfile << "for (size_t i=0; i<multIndv.size(); i++)  \n";
	  mycppfile << "{  \n";
	  mycppfile << "aavGate[*iter] = DoFRef(*spGridFct, multIndv[i]);  \n";
	  for (size_t i=0; i<ListIons.size(); i++)
	  {
	  mycppfile << "aa" + ListIons[i] + "Gate[*iter] = DoFRef(*spGridFct, multInd" + ListIons[i] + "[i]);  \n";
	  }
	  mycppfile << "}  \n \n \n";

	  // Reading functional info for states
	  std::vector<string> INITIAL = GetBlock(Pairs, Zeilen, "INITIAL");
	  std::vector<vector <string> > Proc_funcs = write_proc_block(Pairs, Zeilen);
	  std::vector<string> Proc_vals, locals;
	  string helper;

	  std::vector<string> new_locals, new_local_block;
	  std::vector<string> PROCEDURE = GetBlock(Pairs, Zeilen, "PROCEDURE");
	  size_t begin;
	  size_t Stats_beg = 1;


	  if (INITIAL.size() > 0)
	  {
		  begin = count_beg(INITIAL[1]);
	  }
	  std::cout << "beg: INITIAL!" << std::endl;

	  if (INITIAL.size() > 0)
	  {
		  for (size_t i=0; i<INITIAL.size(); i++)
		  {

			  for (size_t j=0; j<Proc_funcs.size(); j++)
			  {
				  begin = count_beg(INITIAL[i]);
				  if (Proc_funcs[j][0] == INITIAL[i].substr(begin, INITIAL[i].find("(")-begin))
				  {
					  Stats_beg = i;
					  // write needed values
					  for (size_t k=1; k<Proc_funcs[j].size(); k++)
					  {
						  // when only read no write read is same as use
						  if (Only_Read(Pairs, Zeilen, Proc_funcs[j][k])==true)
						  {
							  mycppfile << Write_Only_Read(Pairs, Zeilen, Proc_funcs[j][k]) + ";\n";
					  	  }
					  	  else
					  	  {
					  		  mycppfile << "double " + Proc_funcs[j][k] + " = aa"+ Proc_funcs[j][k] + "Gate[*iter]; \n \n";
					  		  std::cout << "else for: " << std::endl;
					  	  }
					  std::cout << Proc_funcs[j][k] << std::endl;
					  }
					  locals = get_local_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);

					  if (locals.size()>0)
					  {
						  for (size_t k=0; k<locals.size(); k++)
						  {
							  size_t locals_sep;
							  mycppfile << "double " + locals[k] + "; \n";
							  // removes all unneeded formation-styles
							  Remove_all(locals[k]);
							  // locals need seperation through ","
							  locals_sep = locals[k].find(",");
							  if (locals_sep!=0)
							  {
								  while (locals_sep != locals[k].npos)
								  {
									  HFile_added_Vars.push_back(locals[k].substr(0, locals_sep));
									  locals[k].replace(0, locals_sep+1, "");
									  locals_sep = locals[k].find(",");
								  }
								  HFile_added_Vars.push_back(locals[k].substr(0, locals_sep));
							  }


						  }
					  }
					  std::cout << "before proc vals" << std::endl;
					  Proc_vals = writer_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);
					  std::cout << Proc_vals.size() << std::endl;

					  bool HFile_added;
					  for (size_t k=0; k<Proc_vals.size(); k++)
					  {
						  HFile_added = false;
						  if ((Proc_vals[k]!="") && (Proc_vals[k]!="}") && (Proc_vals[k].find("LOCAL")==Proc_vals[k].npos)
							   && Proc_vals[k]!=" " && Proc_vals[k]!="\t" && Proc_vals[k]!="        ")
						  {
  							  for (size_t l=0; l<HFile_added_Vars.size(); l++)
  							  {
  								  //std::cout << "Vergleich " << Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")-1)) << " - " << Remove_all(HFile_added_Vars[l]) << std::endl;
  								  if (Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")-1))==Remove_all(HFile_added_Vars[l]))
  									  HFile_added = true;
  							  }

  							  string com_buf;

  							  // if it added in hfile write without double else with!
  							  if (HFile_added == true)
  							  {
  								  //writes comments
  								  if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
  									  Proc_vals[k].replace(Proc_vals[k].find(":"), 1, "//");
  								  std::cout << Proc_vals[k] << std::endl;
  								  mycppfile << Proc_vals[k] + "; \n";

  							  } else
  							  {
  								// writes comments, if comment for the whole than no double is allowed
  								if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
  								{
  									com_buf = Proc_vals[k];
  									com_buf = Remove_all_com(com_buf);
  									Proc_vals[k].replace(Proc_vals[k].find(":"), 1, "//");
  								}
  								std::cout << "com_buf: " << com_buf << std::endl;

  								if (com_buf.find(":")==0)
  								{
  									std::cout << Proc_vals[k] << std::endl;
  									mycppfile << Proc_vals[k] + "\n";
  								} else
  								{
  									std::cout << Proc_vals[k] << std::endl;
  									mycppfile << "double " + Proc_vals[k] + "; \n";
  								}
  							  }

						  }
					  }


				  }
			  }
		  }		  // write locals vals as double
	  }

	  std::cout<< "writing end" << std::endl;
	  for (size_t i=Stats_beg; i<INITIAL.size(); i++)
	  {
		  // writings at the end
		  if (i >= Stats_beg)
		  {
			  size_t gleich, vorgleich;
			  if ((INITIAL[i]!="") && (INITIAL[i]!="}"))
			  {
				  std::cout << "gleich Search" << std::endl;
				  vorgleich = pos_letter(INITIAL[i]);
				  gleich = INITIAL[i].find("=");

				  if (INITIAL[i].find("=")!=INITIAL[i].npos)
				  {
					  std::cout << "ab: " << vorgleich << "bis: " << gleich << std::endl;
					  mycppfile << "aa" + INITIAL[i].substr(vorgleich , gleich-vorgleich-1) + "Gate[*iter] " + INITIAL[i].substr(gleich) + "; \n";
				  } //else
			    	  //vorgleich = beg_count(INITIAL[i]);
			  }
		  }
	  }



	  std::cout << "INIT happend" << std::endl;

	  mycppfile << "}  \n";
	  mycppfile << "}  \n";
	  mycppfile << "}  \n \n \n \n";

	  mycppfile << "template<typename TDomain, typename TAlgebra> \n";
	  mycppfile << "void " + filename + "<TDomain, TAlgebra>::update_gating(number newTime, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) \n";
	  mycppfile << "{ \n";
	  mycppfile << "typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; \n";
	  mycppfile << "SubsetGroup ssGrp; \n";
	  mycppfile << "try { ssGrp = SubsetGroup(spGridFct->domain()->subset_handler(), this->m_vSubset);} \n";
	  mycppfile << "UG_CATCH_THROW(\"Subset group creation failed.\"); \n \n";
	  mycppfile << "for (std::size_t si = 0; si < ssGrp.size(); si++) \n";
	  mycppfile << "{ \n \n";
	  mycppfile << "itType iterBegin = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);  \n";
	  mycppfile << "itType iterEnd = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); \n \n";
	  mycppfile << "for (itType iter = iterBegin; iter != iterEnd; ++iter) \n";
	  mycppfile << "{ \n \n";
	  mycppfile << "Vertex* vrt = *iter; \n";
	  mycppfile << "// needed konzentration has to be set \n";

	  std::vector<string> DERIVATIVE = GetBlock(Pairs, Zeilen, "DERIVATIVE");
	  //std::vector<vector <string> > Der_funcs;


	  PROCEDURE = GetBlock(Pairs, Zeilen, "PROCEDURE");
	  Proc_funcs = write_proc_block(Pairs, Zeilen);

	  std::vector<vector< string> > Der_funcs;
	  Der_funcs = write_derivative_block(Pairs, Zeilen);


	  if (PROCEDURE.size()>0 && Proc_funcs.size()>0 && Der_funcs.size()>0)
	  {
		  for (size_t i=0; i<PROCEDURE.size(); i++)
			  {
			    std::cout << "Procedure: " << PROCEDURE[i] << std::endl;
			  }
		  for (size_t i=0; i<Proc_funcs[0].size(); i++)
			  {
			  	  for (size_t j=0; j<Proc_funcs.size(); j++)
			    std::cout << "Proc_funcs: " << Proc_funcs[j][i] << std::endl;
			  }

		  if (Der_funcs.size() >= 1)
		  {
			  for (size_t i=0; i<Der_funcs[0].size(); i++)
			  	  {
			  	  	  for (size_t j=0; j<Der_funcs.size(); j++)
			  	  		  std::cout << "Proc_funcs: " << Der_funcs[j][i] << std::endl;
			  	  }
		  }
	  }



	  mycppfile << "size_t vm = spGridFct->fct_id_by_name(\"v\");  \n";
	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "size_t " + ListIons[i] +" = spGridFct->fct_id_by_name(\"" + ListIons[i] + "\"); \n";
	  }

	  mycppfile << "\n \n";

	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "spGridFct->dof_indices(*iter, " + ListIons[i] + ", multInd" + ListIons[i] + "); \n";
	  }

	  mycppfile << "spGridFct->dof_indices(*iter, vm, multIndv);  \n \n ";
	  mycppfile << "for (size_t i=0; i<multIndv.size(); i++)  \n";
	  mycppfile << "{  \n";
	  mycppfile << "aavGate[*iter] = DoFRef(*spGridFct, multIndv[i]);  \n";

	  for (size_t i=0; i<ListIons.size(); i++)
	  {
	 	  mycppfile << "aa" + ListIons[i] + "Gate[*iter] = DoFRef(*spGridFct, multInd" + ListIons[i] + "[i]);  \n";
	  }
	  mycppfile << "}  \n \n \n";


	  // States ever needed for deriv writting out of State_vars State_vars;
	  for (size_t i=0; i<State_vars.size(); i++)
	  {
		  mycppfile << "double " + State_vars[i] + " = aa" + State_vars[i] + "Gate[*iter]; \n";
	  }



	  //Adding v because v is very often needed and always accesible
	  mycppfile << "double v  = aavGate[*iter]; \n";

	  mycppfile << "\n \n \n";

	  std::cout << "DERIV BEGINS" << std::endl;
	  if (DERIVATIVE.size() > 0)
	  	  {
	  		  begin = count_beg(DERIVATIVE[1]);
	  	  }

	  	  if (DERIVATIVE.size() > 0)
	  	  {
	  		  for (size_t i=0; i<DERIVATIVE.size(); i++)
	  		  {

	  			  for (size_t j=0; j<Proc_funcs.size(); j++)
	  			  {
	  				  begin = 0;
	  				  std::cout << "Vergleich: " << DERIVATIVE[i].substr(begin, DERIVATIVE[i].find("(")-begin)<< " - " << Proc_funcs[j][0] << std::endl;
	  				  if (Proc_funcs[j][0] == Remove_all(DERIVATIVE[i].substr(begin, DERIVATIVE[i].find("(")-begin)))
	  				  {
	  					  Stats_beg = i;
	  					  // write needed values
	  					  for (size_t k=1; k<Proc_funcs[j].size(); k++)
	  					  {
	  						  // when only read no write read is same as use
	  						  if (Only_Read(Pairs, Zeilen, Proc_funcs[j][k])==true)
	  						  {
	  							  mycppfile << Write_Only_Read(Pairs, Zeilen, Proc_funcs[j][k]) + ";\n";
	  					  	  }
	  					  	  else
	  					  	  {
	  					  		  mycppfile << "double " + Proc_funcs[j][k] + " = aa"+ Proc_funcs[j][k] + "Gate[*iter]; \n \n";
	  					  		  std::cout << "else for: " << std::endl;
	  					  	  }
	  					  std::cout << Proc_funcs[j][k] << std::endl;
	  					  }
	  					  locals = get_local_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);

	  					  if (locals.size()>0)
	  					  {
							  for (size_t k=0; k<locals.size(); k++)
							  {
								  size_t locals_sep;
								  mycppfile << "double " + locals[k] + "; \n";
								  // removes all unneeded formation-styles
								  Remove_all(locals[k]);
								  // locals need seperation through ","
								  locals_sep = locals[k].find(",");
								  if (locals_sep!=0)
								  {
									  while (locals_sep != locals[k].npos)
									  {
										  HFile_added_Vars.push_back(locals[k].substr(0, locals_sep));
										  locals[k].replace(0, locals_sep+1, "");
										  locals_sep = locals[k].find(",");
									  }
									  HFile_added_Vars.push_back(locals[k].substr(0, locals_sep));
								  }


							  }
	  					  }
	  					  std::cout << "before proc vals" << std::endl;
	  					  Proc_vals = writer_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);
	  					  std::cout << Proc_vals.size() << std::endl;


	  					  bool HFile_added;
	  					  for (size_t k=0; k<Proc_vals.size(); k++)
	  					  {
							  HFile_added = false;
							  if ((Proc_vals[k]!="") && (Proc_vals[k]!="}") && (Proc_vals[k].find("LOCAL")==Proc_vals[k].npos)
								   && Proc_vals[k]!=" " && Proc_vals[k]!="\t" && Proc_vals[k]!="        ")
							  {
	  							  for (size_t l=0; l<HFile_added_Vars.size(); l++)
	  							  {
	  								  //std::cout << "Vergleich " << Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")-1)) << " - " << Remove_all(HFile_added_Vars[l]) << std::endl;
	  								  if (Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")-1))==Remove_all(HFile_added_Vars[l]))
	  									  HFile_added = true;
	  							  }

	  							  string com_buf;

	  							  // if it added in hfile write without double else with!
	  							  if (HFile_added == true)
	  							  {
	  								  //writes comments
	  								  if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
	  									  Proc_vals[k].replace(Proc_vals[k].find(":"), 1, "//");

	  								  mycppfile << Proc_vals[k] + "; \n";

	  							  } else
	  							  {
	  								// writes comments, if comment for the whole than no double is allowed
	  								if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
	  								{
	  									com_buf = Proc_vals[k];
	  									com_buf = Remove_all_com(com_buf);
	  									Proc_vals[k].replace(Proc_vals[k].find(":"), 1, "//");
	  								}

	  								if (com_buf.find(":")==0)
	  								{
	  									mycppfile << Proc_vals[k] + "\n";
	  								} else
	  								{
	  									mycppfile << "double " + Proc_vals[k] + "; \n";
	  								}
	  							  }

							  }
	  					  }


	  				  }
	  			  }
	  		  }		  // write locals vals as double
	  	  }

	  	  std::cout<< "writing end" << std::endl;
	  	  for (size_t i=Stats_beg; i<DERIVATIVE.size(); i++)
	  	  {
	  		  // writings at the end
	  		  if (i >= Stats_beg)
	  		  {
				  // Writting derivfuncs
				  std::cout << "before deriv" << std::endl;
				  if ((i> Stats_beg) && (i < DERIVATIVE.size()-1))
				  {
					if (DERIVATIVE[i]!="")
					{
						DERIVATIVE[i].replace(DERIVATIVE[i].find("\'"), 1, "");
						DERIVATIVE[i].replace(DERIVATIVE[i].find("="), 1, "+=");
						mycppfile << DERIVATIVE[i] + "*dt; \n";
					}
				  }
	  		  }
	  	  }

	  mycppfile << "\n \n \n";



	  std::cout << "DERIV happened" << std::endl;


	  std::cout << "Adding kinetic funcs" << std::endl;

	  string needed_proc;
	  std::vector<string> KINETIC;
	  KINETIC = GetBlock(Pairs, Zeilen, "KINETIC");
	  // KINETIC used always mentioned in Breakpoint
	  std::vector<string> BREAKPOINT = GetBlock(Pairs, Zeilen, "BREAKPOINT");

	  std::vector<string> Solve_List;
	  size_t beg_sol;
	  // write solve List
	  for (size_t i=0; i< BREAKPOINT.size(); i++)
	  {
		  beg_sol = BREAKPOINT[i].find("SOLVE");
		  if (beg_sol != BREAKPOINT[i].npos)
		  {
			  needed_proc = BREAKPOINT[i-1];
			  std::cout << "first needed proc: " << needed_proc << std::endl;
			  Solve_List.push_back(BREAKPOINT[i].substr(beg_sol+6, BREAKPOINT[i].find("METHOD")-1-(beg_sol+6)));
		  }
	  }

	  size_t kin_help_beg, kin_help_end;
	  string kin_helpS;
	  vector<string> Ausgabe, Ausgabe2;


	  // if kinectic is defined
	  if (KINETIC.size() >=1)
	  {

		  if (Solve_List.size()>=1)
		  {
			  for (size_t i=0; i<Solve_List.size(); i++)
			  {
				  kin_help_beg = KINETIC[0].find("KINETIC");
				  kin_help_end = KINETIC[0].find("{");
				  kin_helpS = KINETIC[0].substr(kin_help_beg+7, kin_help_end-1-(kin_help_beg+7));
				  while (kin_helpS.find(" ")!=kin_helpS.npos)
					  kin_helpS.replace(kin_helpS.find(" "), 1, "");
				  if (kin_helpS == Solve_List[i])
				  {
					  /*std::cout << "da geht was" << std::endl;
					  std::cout << kin_helpS << std::endl;
					  std::cout << Solve_List[i] << std::endl;*/
					  for (size_t j=0; j<KINETIC.size(); j++)
					  {

						  vector<string> varsL;
						  string left, right, vars;
						  size_t komma;
						  size_t helpco, helpcc;
						  size_t c_f_o, c_l_c;
						  size_t arrows = KINETIC[j].find("<->");
						  size_t KIN_length = KINETIC[j].npos;
						  if (arrows!=KIN_length)
						  {
							  left = KINETIC[j].substr(0, arrows-1);
							  right = KINETIC[j].substr(arrows+3, KIN_length-(arrows+3));

							  c_f_o = right.find("(");
							  c_l_c = right.find_last_of(")");

							  vars = right.substr(c_f_o+1, c_l_c-(c_f_o+1));
							  right = right.substr(0, c_f_o-1);

							  std::cout << "left: " << left <<" right: "<< right << "vars: " << vars << std::endl;
							  helpco = vars.find("(");
							  helpcc = vars.find(")");
							  // del all clips
							  while ((helpco!=vars.npos) || (helpcc!=vars.npos))
							  {
								  vars.replace(helpco ,1 ,"");
								  vars.replace(helpcc ,1 ,"");
								  helpco = vars.find("(");
								  helpcc = vars.find(")");
							  }
							  // divide vars and push in var list
							  size_t beg = 0;
							  komma = vars.find(",");
							  while (komma!=vars.npos)
							  {
								  //std::cout << "var pushed_back: " << vars.substr(beg, komma-beg) << std::endl;
								  varsL.push_back(vars.substr(beg, komma-beg));
								  beg = komma+1;
								  komma = vars.find(",", beg);
							  }
							  //std::cout << "var pushed_back: " << vars.substr(beg, komma-beg) << std::endl;
							  varsL.push_back(vars.substr(beg, komma-beg));

							  Ausgabe.push_back(left + "+ =  (-" +left +"*" + varsL[0] + "+" +right +"*" + varsL[1] + ")");
							  Ausgabe.push_back(right + "+ = ("+left +"*" + varsL[0] + "+ -"+ right +"*" + varsL[1] + ")");




						  }
					  }
				  }
			  }
		  }

		  //del all " " and tabs out of name

		  size_t leer, tab;
		  leer = needed_proc.find(" ");
		  tab = needed_proc.find("\t");
		  while (leer != needed_proc.npos)
		  {
			  needed_proc.replace(leer, 1, "");
			  leer = needed_proc.find(" ");
		  }

		  while (tab != needed_proc.npos)
		  {
			  needed_proc.replace(tab, 1, "");
			  tab = needed_proc.find(" ");
		  }

		  //getting procedure for needed vars
		  size_t clip = needed_proc.find(")");
		  needed_proc.replace(clip, 1, "");
		  Ausgabe2 = GetProcEqualString(Pairs, Zeilen, needed_proc);
		  for (size_t i=0; i<Ausgabe2.size(); i++)
		  {
			  if (Ausgabe2[i].find("PROCEDURE")==Ausgabe2[i].npos)
				  mycppfile << Ausgabe2[i] + "; \n";
		  }

		  mycppfile << " \n \n \n";



		  for (size_t i=0; i<Ausgabe.size(); i++)
		  {
			  // removes all " " tabs and ~
			  leer = Ausgabe[i].find(" ");
			  while (leer!=Ausgabe[i].npos)
			  {
				  Ausgabe[i].replace(leer, 1, "");
				  leer = Ausgabe[i].find(" ");
			  }
			  leer = Ausgabe[i].find("~");
			  while (leer!=Ausgabe[i].npos)
			  {
				  Ausgabe[i].replace(leer, 1, "");
				  leer = Ausgabe[i].find("~");
			  }
			  leer = Ausgabe[i].find("\t");
			  while (leer!=Ausgabe[i].npos)
			  {
				  Ausgabe[i].replace(leer, 1, "");
				  leer = Ausgabe[i].find("\t");
			  }
			  mycppfile << Ausgabe[i] +"*dt; \n";
		  }

		  mycppfile << " \n \n \n";



	// if there is no kinetic
	  } else
	  {
		  vector<string> BProc;
		  size_t breakS;
		  string Proc;
		  if (BREAKPOINT.size()>0)
		  {
			  for (size_t i=0; i<BREAKPOINT.size(); i++)
			  {
				  breakS = BREAKPOINT[i].find("SOLVE");
				  if (breakS!=BREAKPOINT[i].npos)
				  {
					  Proc = BREAKPOINT[i].substr(breakS+6 ,BREAKPOINT[i].find(" ")-(breakS+6));
					  std::cout << "PROC: " << Proc << std::endl;
					  BProc = writer_proc_block(Pairs, Zeilen, Proc);
				  }
			  }


			  //updates States for use later
			  if (State_vars.size()>0)
			  {
			  	for (size_t i=0; i<State_vars.size(); i++)
			  	{
			  		HFile_added_Vars.push_back(State_vars[i]);
			  	}
			  	HFile_added_Vars.push_back("v");
			  }





			  bool HFile_added;
			  for (size_t i=0; i<BProc.size(); i++)
			  {
				  HFile_added = false;
				  if ((BProc[i]!="") && (BProc[i]!="}") && (BProc[i].find("LOCAL")==BProc[i].npos)
				  	   && BProc[i]!=" " && BProc[i]!="\t" && BProc[i]!="        ")
				  {
					  for (size_t l=0; l<HFile_added_Vars.size(); l++)
				  	  {
						  //std::cout << "Vergleich " << Remove_all(BProc[i].substr(0, BProc[i].find("=")-1)) << " - " << Remove_all(HFile_added_Vars[l]) << std::endl;
						  if (Remove_all(BProc[i].substr(0, BProc[i].find("=")-1))==Remove_all(HFile_added_Vars[l]))
							  HFile_added = true;
				  	  }

				  	  string com_buf;

				  	  // if it added in hfile write without double else with!
				  	  if (HFile_added == true)
				  	  {
				  		  //writes comments
				  		  if (BProc[i].find(":")!=BProc[i].npos)
				  			  BProc[i].replace(BProc[i].find(":"), 1, "//");

				  		  mycppfile << BProc[i] + "; \n";

				  	  } else
				  	  {
				  		  // writes comments, if comment for the whole than no double is allowed
				  		  if (BProc[i].find(":")!=BProc[i].npos)
				  		  {
				  			  com_buf = BProc[i];
				  			  com_buf = Remove_all_com(com_buf);
				  			  BProc[i].replace(BProc[i].find(":"), 1, "//");
				  		  }

				  		  if (com_buf.find(":")==0)
				  		  {
				  			  mycppfile << BProc[i] + "\n";
				  		  } else
				  		  {
				  			  mycppfile << "double " + BProc[i] + "; \n";
				  		  }
				  	  }

				  }
			  }
		  }

	  }








//updates Gattings
	  if (State_vars.size()>0)
	  {
		  for (size_t i=0; i<State_vars.size(); i++)
		  {
			  mycppfile << "aa" + State_vars[i] + "Gate[*iter] = " + State_vars[i] + "; \n";
		  }
	  }



	  mycppfile << " \n \n \n";
	  mycppfile << "} \n";
	  mycppfile << "} \n";
	  mycppfile << "} \n";
	  mycppfile << " \n \n \n";

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// write head of ionic flux
//////////////////////////////////////////////////////////////////////////////////////////////////////////
	  mycppfile << "template<typename TDomain, typename TAlgebra> \n";
	  mycppfile << "void " + filename + "<TDomain, TAlgebra>::ionic_current(Vertex* ver, std::vector<number>& outCurrentValues) \n";
	  mycppfile << "{ \n \n";
	  // all gates

	  if (State_vars.size()>0)
	  {
		  // States ever need new values
		  	  for (size_t i=0; i<State_vars.size(); i++)
		  	  {
		  		  mycppfile << "double " + State_vars[i] + " = aa" + State_vars[i] + "Gate[ver]; \n";
		  	  }
	  }

	  if (ListIons.size()>0)
	  {
		  // all ions
		  for (size_t i=0; i<ListIons.size(); i++)
		  {
		 	  mycppfile << "double " + ListIons[i] + " = aa" + ListIons[i] + "Gate[ver];  \n";
		  }
	  }

	  // v
	  mycppfile << "double v = aavGate[ver]; \n";
	  mycppfile << " \n";
	  mycppfile << " \n";

	  //writes Nernst Eqs
	  std::vector<string> eqs = equali(Pairs, Zeilen);

	  if (eqs.size()>0)
	  {
		  for (size_t i=0; i<eqs.size(); i++)
		  {
			  // TODO write for schleife over all ions!
			  if (eqs[i]!="1" && eqs[i]!="0" && eqs[i]!=" ca" && eqs[i]!=" na" && eqs[i]!=" k")
				  mycppfile << eqs[i] + " \n";
		  }
	  }

	  mycppfile << " \n \n";

	  bool in = true;
	  bool stand_BP = true;

	  size_t beg;
	  std::vector<string> fluxes;
	  std::vector<string> B_vars;


	  std::cout << "Before Breakpoint" << std::endl;

	  //decids which Breakpoint Method will be used
	  for (size_t i=1; i<BREAKPOINT.size()-1; i++)
	  {
		  if (BREAKPOINT[i].find("if")!=BREAKPOINT[i].npos)
		  {
			  stand_BP = false;
		  }
	  }

	  string non_spec_current;
	  size_t non_spec_beg, non_spec_end;

	  // writes Breakpoint Word for Word add out at the end
	  if (stand_BP==false)
	  {
		  //search for non spec_ionic_current
		  for (size_t i=0; i<NEURON.size(); i++)
		  {
			  non_spec_beg = NEURON[i].find("NONSPECIFIC_CURRENT");
			  if (non_spec_beg!=NEURON[i].npos)
			  {
				  non_spec_end = NEURON[i].find(" ", non_spec_beg+20);
				  non_spec_current = NEURON[i].substr(non_spec_beg+20, non_spec_end-(non_spec_beg+20));
			  }
		  }

		  string test;
		  for (size_t i=1; i<BREAKPOINT.size()-1; i++)
		  {
			  test = BREAKPOINT[i] + "; \n";
			  if ((test.find("{;")==test.npos) && (test.find("};")==test.npos))
			  {
			  	  mycppfile << BREAKPOINT[i] + "; \n";
			  }
			  else
			  {
				  mycppfile << BREAKPOINT[i] + "\n";
		 	  }
		  }

		  mycppfile << "outCurrentValues.push_back(" + non_spec_current + "); \n";

	  }






	  // writes Breakpoint if it is Standard BP
	  if (stand_BP==true)
	  {
		  for (size_t i=1; i<BREAKPOINT.size()-1; i++)
		 	  {
		 		  in = true;
		 		  if (BREAKPOINT[i].find("SOLVE")!=BREAKPOINT[i].npos)
		 			  in = false;

		 		  if (in == true)
		 		  {
		 			  std::cout << "adding flux" << std::endl;
		 			  beg = pos_letter(BREAKPOINT[i]);
		 			  if (beg!=1000)
		 				  fluxes.push_back(BREAKPOINT[i].substr(beg, BREAKPOINT[i].npos-beg));
		 		  }

		 	  }
		 	  std::cout << "fluxes size: " << fluxes.size() << std::endl;

		 	  for (size_t i=0; i<fluxes.size(); i++)
		 	  {
		 		  // make list with left handed vars
		 		  B_vars.push_back(fluxes[i].substr(0, fluxes[i].find("=")-1));
		 		  std::cout << (B_vars[i]) << std::endl;
		 	  }
		 	  // check if vars on right side an write down
		 	  for (size_t i=0; i<fluxes.size(); i++)
		 	  {
		 		  for (size_t j=0; j<B_vars.size(); j++)
		 		  {
		 			  beg = fluxes[i].find("=");
		 			  if (fluxes[i].find(B_vars[j], beg) != fluxes[i].npos)
		 			  {
		 				  size_t test = fluxes[i].find(B_vars[j], beg);
		 				  if (pos_letterb(fluxes[i].substr(test+B_vars[j].size(), 1))==false)
		 				  {
		 				  	  mycppfile << "number " + fluxes[i] + "; \n";
		 				  	  fluxes[i] = "";
		 				  }
		 			  }

		 		  }
		 	  }

		 	  vector<bool> Ion_outside;
		 	  vector<string> outs;
		 	  // first every time v
		 	  outs.push_back("");

		 	  int vm_flux_count = 0;

		 	  mycppfile << "\n \n \n";
		 	  for (size_t i=0; i<fluxes.size(); i++)
		 	  {
		 		  std::cout << "fluxes[i]: " << fluxes[i] << std::endl;

		 		  if (fluxes[i]!="")
		 		  {
		 			  if (fluxes[i].find("i")==0)
		 			  {
		 				  if (vm_flux_count>=1)
		 				  {
		 					  outs[0] += " + " + fluxes[i].substr(fluxes[i].find("=")+1, fluxes[i].npos -fluxes[i].find("=")+1);
		 					  vm_flux_count += 1;
		 				  } else {
		 					  outs[0] += fluxes[i].substr(fluxes[i].find("=")+1, fluxes[i].npos -fluxes[i].find("=")+1);
		 					  vm_flux_count += 1;
		 				  }
		 			  }
		 		  }
		 	  }

		 	  for (size_t i=0; i<fluxes.size(); i++)
		 	  {
		 		  for (size_t j=0; j<ListIons.size(); j++)
		 		  {
		 			  Ion_outside.push_back(false);
		 			  std::cout << "List Ions" << std::endl;
		 			  std::cout << "i"+ListIons[j] << std::endl;
		 			  // ToDo here we need infos out of equali
		 			  string change;
		 			  for (size_t k=0; k<eqs.size(); k++)
		 			  {
		 				  if (eqs[k]=="0")
		 				  {
		 					  size_t eqs_leer = eqs[k+1].find(" ");
		 					  while (eqs_leer!=eqs[k+1].npos)
		 					  {
		 						  eqs[k+1].replace(eqs_leer, 1, "");
		 						  eqs_leer = eqs[k+1].find(" ");
		 						  std::cout << "eqs k+1: " << eqs[k+1] << std::endl;
		 					  }
		 					  if (eqs[k+1]==ListIons[j])
		 					  {
		 						  //setting Ion_outside true to prevent output later as push_back from ionic_current
		 						  //cause ion is only working on outside
		 						  Ion_outside[j]=true;

		 						  size_t write, komma;
		 						  string left, right;
		 						  // create out of NEURON another output depending on WRITE Part!
		 						  std::cout << "changes will happen on outer concentration" << std::endl;
		 						  for (size_t l=0; l<NEURON.size(); l++)
		 						  {
		 							  if (NEURON[l].find("USEION " + eqs[k+1])!=NEURON[l].npos)
		 							  {
		 								 komma = NEURON[l].find(",");
		 								 write = NEURON[l].find("WRITE");
		 								 left = NEURON[l].substr(komma+1, write-1);
		 								 right = NEURON[l].substr(write+6, NEURON[l].npos);
		 								 if (fluxes[i].find(left)!=fluxes[i].npos)
		 								 {
		 									 mycppfile << eqs[k+1] + fluxes[i].substr(fluxes[i].find("="), fluxes[i].npos-fluxes[i].find("=")) + "; \n";
		 								 }

		 								 if (fluxes[i].find(right)!=fluxes[i].npos)
		 								 {
		 									 mycppfile << eqs[k+1] + "_out " + fluxes[i].substr(fluxes[i].find("="), fluxes[i].npos-fluxes[i].find("="))+ "; \n";
		 								 }

		 							  }
		 						  }

		 					  }
		 				  }


		 			  }
		 			  if ((fluxes[i].find("i" + ListIons[j])!=fluxes[i].npos) && (Ion_outside[j]==false))
		 				  outs.push_back(fluxes[i].substr(fluxes[i].find("=")+1, fluxes[i].npos -fluxes[i].find("=")+1));
		 		  }
		 	  }


		 	  mycppfile << "outCurrentValues.push_back(" + outs[0] + "); \n";
		 	  for (size_t i=1; i<outs.size(); i++)
		 	  {
		 		  mycppfile << "outCurrentValues.push_back(" + outs[i] + "/m_F ); \n";
		 	  }
	  }

	  mycppfile << "  \n";
	  mycppfile << "}  \n";



	  mycppfile << "  \n";
	  mycppfile << "  \n";
	  mycppfile.close();


	  myhfile << "ADouble v; \n";
	  myhfile << "Grid::AttachmentAccessor<Vertex, ADouble> aav; \n \n";
	  myhfile << "//nernst const values \n";

	  myhfile << "number m_R; \n";
	  myhfile << "number m_T; \n";
	  myhfile << "number m_F; \n \n";

	  myhfile << "/// Base type \n";
	  myhfile << "typedef IChannel<TDomain, TAlgebra> base_type; \n";
	  myhfile << "///	Own type \n";
	  myhfile << "typedef ChannelHH<TDomain, TAlgebra> this_type; \n";
	  myhfile << "/// GridFunction type \n";
	  myhfile << "typedef GridFunction<TDomain, TAlgebra> TGridFunction; \n \n";

	  myhfile << "}; \n \n";

	  myhfile << "} // namespace ug \n \n \n";

	  myhfile << "#endif // " + filename + "_H_\n";

	  myhfile.close();
}

