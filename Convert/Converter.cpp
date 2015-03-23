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
#include <unistd.h>
#include <sstream>

//Todo add setter for generator
// closer loook to unit conv
// perhaps java interface for converter

Converter::Converter() {
	// TODO Auto-generated constructor stub

}

Converter::~Converter() {
	// TODO Auto-generated destructor stub
}


// Giving back all needed Nernst-Potentials
std::vector<string> Converter::Find_all_Eqs(std::vector<string> Neuron)
{
	std::vector<string> Erg;
	std::vector<string> Read  = Find_all_read(Neuron);
	std::vector<string> Write = Find_all_write(Neuron);

	for (size_t i=0; i<Read.size(); i++)
	{
		string test_e = Remove_all(Read[i]);
		if (test_e[0]=='e')
		{
			Erg.push_back(Read[i]);
		}
	}

	for (size_t i=0; i<Write.size(); i++)
	{
		string test_e = Remove_all(Write[i]);
		if (test_e[0]=='e')
		{
			Erg.push_back(Write[i]);
		}
	}

	return Erg;
}

// Tests if a Ion name also is a gatting name and gives back GatingName of Ion
// here working with complete string out of DERIV - Block
string Converter::ProofSGatingGating(string s, std::vector<string> SGating)
{
	string Gating_Name;
	string erg = "";
	string left, right;
	size_t rBegin;
	// IF there is any = else all marked as right side
	if (s.find("=")!=s.npos)
	{
		left = s.substr(0, s.find("=")-1);
		right = s.substr(s.find("=")+1, s.npos-(s.find("=")+1));
	}
	else
	{
		left = "";
		right = s;
	}

	// left side easy to work like cause only one var
	for (size_t i=0; i<SGating.size(); i++)
	{
		if (Remove_all(left) + "S"== SGating[i])
			left = Remove_all(left) + "S";
	}

	//std::cout << "after first for" << std::endl;
	// handling right side
	string buf_right = right;
	for (size_t i=0; i<SGating.size(); i++)
	{
		Gating_Name = SGating[i].substr(0, SGating[i].size()-1);
		rBegin = right.find(Gating_Name);
		while (rBegin!=right.npos)
		{
			if (is_single_letter(right[rBegin+Gating_Name.size()]) == false && is_single_number(right[rBegin+Gating_Name.size()]) == false)
			{
				right = right.replace(rBegin, Gating_Name.size(), SGating[i]);
				rBegin = right.find(Gating_Name, rBegin+1);
			}
			else
			{
				rBegin = right.find(Gating_Name, rBegin+1);
			}
		}
	}
	erg = left + " += " + right;

	return erg;
}


// Tests if a Ion name also is a gatting name and gives Back gattingAtachment of Gating with Ion name
// Here only working with the var itself
string Converter::ProofSGatingInit(string s, std::vector<string> SGating)
{
	string erg = Remove_all(s);


	for (size_t i=0; i<SGating.size(); i++)
	{
		//std::cout << "ProofSGating: " << erg << " vs " << SGating[i] << std::endl;
		if (erg + "S"== SGating[i])
			erg = "aa" + erg + "SGate[vrt]";
	}


	return erg;
}

//giving back all vars which will be changed throug new Channel-Interface
std::vector<string> Converter::Find_all_write(std::vector<string> Neuron)
{
	std::vector<string> result;
	string all_writes;
	size_t komma, begin;

	for (size_t i=0; i<Neuron.size(); i++)
	{
		if (Neuron[i].find("USEION")!=Neuron[i].npos)
		{
			if (Neuron[i].find("WRITE")!=Neuron[i].npos)
			{
				all_writes = Neuron[i].substr(Neuron[i].find("WRITE")+6, Neuron[i].npos-(Neuron[i].find("WRITE")+6));
				komma = all_writes.find(",");
				begin = 0;
				while (komma!=all_writes.npos)
				{
					result.push_back(Remove_all(all_writes.substr(begin, komma)));
					begin = komma;
					komma = all_writes.find(",", begin+1);
				}
				result.push_back(Remove_all(all_writes.substr(begin, komma)));
			}
		}
	}

	return result;
}

// giving back all information which is read from Neuron
std::vector<string> Converter::Find_all_read(std::vector<string> Neuron)
{
	std::vector<string> result;
	string all_reads;
	size_t komma, begin;

	for (size_t i=0; i<Neuron.size(); i++)
	{
		if (Neuron[i].find("USEION")!=Neuron[i].npos)
		{
			if (Neuron[i].find("READ")!=Neuron[i].npos)
			{
				all_reads = Neuron[i].substr(Neuron[i].find("READ")+5, Neuron[i].find("WRITE")-1-(Neuron[i].find("READ")+5));
				komma = all_reads.find(",");
				begin = 0;
				while (komma!=all_reads.npos)
				{
					result.push_back(Remove_all(all_reads.substr(begin, komma)));
					// +1 neede for ","
					begin = komma+1;
					komma = all_reads.find(",", begin+1);
				}
				result.push_back(Remove_all(all_reads.substr(begin, komma)));
			}
		}
	}

	return result;
}

double Converter::Unit_Conv_All(string s)
{
	double erg = 1;
	double break_dash;
	double above, below;

	break_dash = s.find("/");
	//std::cout << "break dash: " << break_dash << std::endl;
	if (break_dash!=s.npos)
	{
		//std::cout << "above string: " << s.substr(0, break_dash) <<std::endl;
		//std::cout << "below string: " << s.substr(break_dash+1, s.npos) <<std::endl;
		above = Unit_Conv_Value(s.substr(0, break_dash));
		below = Unit_Conv_Value(s.substr(break_dash+1, s.npos));
		erg = (above/below);
	}
	else
	{
		erg = Unit_Conv_Value(s);
	}
	//std::cout << "above: " << above << std::endl;
	//std::cout << "below: " << below << std::endl;



	return erg;
}


double Converter::Unit_Conv_Value(string s)
{
	double erg = 1;
	//needed units /c /m^2 /mV /ms
	// existing units amper: pA, mA, A
	// existing units volt: pV, mV, V
	// um2 mm2 cm2 dm2 m2 km2
	// existing units time: ms, s,
	// working with Siemens!!
	// existing unigs simens: pS, mS, S
	// existing units concentrats: mM,
	// existing units temp: degC
	// (um) = (micron) ????
	// mho/cm2 mA/cm2 S/cm2
	// 2 standing for ²

	bool onlyVA=true;
	// In our model we have mA and mV
	if ((s.find("V")!=s.npos) && (s.find("A")!=s.npos))
	{
		if ((s.find("pV")!=s.npos) || (s.find("pA")!=s.npos))
		{
			if (erg==0)
				erg = 1e-9;
			else erg *= 1e-9;

			onlyVA = false;
		}

		if ((s.find("nV")!=s.npos) || (s.find("nA")!=s.npos))
		{
			if (erg==0)
				erg = 1e-6;
			else erg *= 1e-6;

			onlyVA = false;
		}
		if (onlyVA == true)
		{
			if ((s.find("V")!=s.npos) || (s.find("A")!=s.npos))
			{
				if (erg==0)
					erg = 1e3;
				else erg *= 1e3;
			}
		}
	}

	//special mho is siemens
	// In our model we need c/m^2/mV/ms
	// Siemens = 1A/V
	// coloumb = 1A * 1s
	// Siemens = 1000*1000/1000
	bool onlyS=true;
	if ((s.find("S")!=s.npos))
	{
		if ((s.find("pS")!=s.npos))
		{
			if (erg==0)
				erg = 1e15;
			else erg *= 1e15;
			onlyS = false;
		}

		if ((s.find("nS")!=s.npos))
		{
			if (erg==0)
				erg = 1e12;
			else erg *= 1e12;
			onlyS = false;
		}



		if ((s.find("mS")!=s.npos))
		{
			if (erg==0)
				erg = 1e6;
			else erg *= 1e6;
			onlyS = false;
		}

		if (onlyS == true)
		{
			if (erg==0)
				erg = 1e3;
			else erg *= 1e3;
		}
	}





	// In our model we have m/m2/m3/m4
	string dims[] = {"", "2", "3", "4"};
	int dimi[] = {1, 100, 10000, 1000000};
	double factor = 0.0;
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
			factor = i+1.0;

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


//Gives Ion back which is in kontakt with string s for example if s = ena than erg = na
string Converter::In_NeuronUse_List(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s)
{
	string erg ="";
	size_t end;

	std::vector<string> NEURON = GetBlock(Pairs, Zeilen, "NEURON");

	for (size_t i=0; i<NEURON.size(); i++)
	{
		if (NEURON[i].find("USEION")!=NEURON[i].npos)
		{
			if (NEURON[i].find(Remove_all(s))!=NEURON[i].npos)
			{
				if (NEURON[i].find("READ")!=NEURON[i].npos)
				{
					end = NEURON[i].find("READ")-1;
				}
				else
				{
					if (NEURON[i].find("WRITE")!=NEURON[i].npos)
					{
						end = NEURON[i].find("WRITE")-1;
					}
					else
					{
						end = NEURON[i].npos;
					}
				}
				// writting who it is called in our prog else ""
				erg = Remove_all(NEURON[i].substr(NEURON[i].find("USEION")+7, end-(NEURON[i].find("USEION")+7)));
			}
		}

	}

	return erg;

}

// tests if single char is a letter
bool Converter::is_single_letter(char s)
{
	// little letters
	for ( char i( 97 ); i < 122; ++i )
	{
		//97-122
		if (s==i)
			return true;
	}

	// big letters
	for ( char i( 65 ); i < 90; ++i )
	{
		//97-122
		if (s==i)
			return true;
	}
	return false;
}

// test if single char is a number
bool Converter::is_single_number(char s)
{
	bool erg = false;

	if (s=='0' || s=='1' || s=='2' || s=='3' || s=='4' || s=='5' || s=='6' || s=='7' || s=='8'  || s=='9' || s=='.')
		erg = true;

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
		////std::cout << "vergleich von " << test << "mit: " << std::endl;
		//std::cout << PROCEDURE[i] << std::endl;
		if (PROCEDURE[i].find(test)!=PROCEDURE[i].npos)
		{
			beg = i;
			//std::cout << "gleich" << std::endl;
		}
		if ((PROCEDURE[i].find("PROCEDURE")!=PROCEDURE[i].npos) && (i!=beg))
		{
			end = i;
			i = PROCEDURE.size();
		}
	}
	if (end == 0)
		end = PROCEDURE.size()-1;

	//std::cout << "begin: " << beg << " ende: " << end << std::endl;

	//if (beg != 0)
	//{
		//std::cout << "GetProcEqual does anything" << std::endl;
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
				//std::cout << "beg: " << pos+7 << " dur: " << end-pos+7 << "end: "<< end << std::endl;
				erg = erg + NEURON[i].substr(pos+7, (end-(pos+7))-1)+"Gate[*iter]";
			}
		}
	}
	return (erg + "; \n");
}


bool Converter::begG(string s)
{
	////std::cout << "fehler in beg" << std::endl;
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
	out.push_back("const number helpV = 1e3*(m_R*m_T)/m_F;");

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
				//std::cout << IonRend << ", " << Ion << std::endl;
				IonS = NEURON[i].substr(Ion+7, IonRead-1-(Ion+7));
				//std::cout << "Subion: "<< IonS << std::endl;
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
			//std::cout << "read ion: " + ListIonRead[j] + "List ion: " + ListIon[j] << std::endl;
			out.push_back("number " + ListIonRead[j] + " = helpV*(log(m_pVMDisc->" + ListIon[j] + "_out/" + ListIon[j] + "));");
			}
		}


		for (size_t j=0; j<ListIon.size(); j++)
		{
			if (version[j]==true)
			{
			//std::cout << "read ion: " + ListIonRead[j] + "List ion: " + ListIon[j] << std::endl;
			out.push_back("number " + ListIonRead[j] + " = helpV*(log(m_pVMDisc->" + ListIon[j] + "_out/" + ListIon[j] + "));");
			}
		}



	return out;
}

// making out of complete function head only function name
string Converter::func_head(string s)
{
	string fname;
	size_t startf, endf;
	std::vector<string> vars;
	startf = s.find(" ");
	endf = s.find("(");
	fname = s.substr(startf+1, endf-startf-1);
	//std::cout << "fname: " << fname << std::endl;
	return fname;
}

// Writing heads an vars of all Procedures which where opend in Procedure-Block
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

// Unit change not really needed but could be better with
// Writing heads an vars of all Procedures which where opend in Derivative-Block
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



// giving back all local vars in a Procedure-block
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


// search for different vars which not needed
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



// Very import function which is recursive
std::vector<string> Converter::writer_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name)
{
	std::vector<string> fileInfo;
	std::vector<string> vars;
	std::vector<string> fileInfo_norm;
	std::vector< vector <string> > fileInfo_rec;
	size_t beg;

	std::vector<vector <string> > Proc_list;

	std::vector<string>  PROCEDURE = GetBlock(Pairs, Zeilen, "PROCEDURE");

	// only [j][0] have functionnames
	Proc_list = write_proc_block(Pairs, Zeilen);

	for (size_t i = 0; i<PROCEDURE.size(); i++)
	{
		beg = PROCEDURE[i].find("PROCEDURE " + name);
		if (beg!=PROCEDURE[i].npos)
		{
			for (size_t j = i; j<PROCEDURE.size(); j++)
				{
					// search for first line to write
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
							////std::cout << "zeug mal wieder: " << std::endl;
							////std::cout << PROCEDURE[k].substr(PROCEDURE[k].find("(")+1, PROCEDURE[k].find(")")-PROCEDURE[k].find("(")-1) << std::endl;
							////std::cout << Proc_list[l][1] << std::endl;
							fileInfo.push_back(Proc_list[l][1] + " = " + PROCEDURE[k].substr(PROCEDURE[k].find("(")+1, PROCEDURE[k].find(")")-PROCEDURE[k].find("(")-1));
							fileInfo_rec.push_back(writer_proc_block(Pairs, Zeilen, Proc_list[l][0]));
							//fileInfo = writer_proc_block(Pairs, Zeilen, Proc_list[l][0]);
							rec = true;

						}

					}

					// writing everything that is not recursiv
					if (rec==false)
						//fileInfo.push_back(PROCEDURE[k]);
						fileInfo_norm.push_back(PROCEDURE[k]);

				}
		}
	}



	bool is_added = false;

	// adding all out of recursion in fileInfo output
	for (size_t l=0; l<fileInfo_rec.size(); l++)
	{
		for (size_t m=0; m<fileInfo_rec[l].size(); m++)
		{
			is_added = false;
			// adde es nur wenn es noch nicht vorhanden ist
			for (size_t n=0; n<fileInfo.size(); n++)
			{
				// if it is already written do not write again
				if (fileInfo[n]==fileInfo_rec[l][m])
					is_added = true;

			}
			if (is_added == false)
				fileInfo.push_back(fileInfo_rec[l][m]);

		}

	}

	// adding all of not recursiv
	for (size_t l=0; l<fileInfo_norm.size(); l++)
	{
		fileInfo.push_back(fileInfo_norm[l]);
	}

	return fileInfo;
}


//searching for all Vars in a Procedure call
//with given PROCEDURE string
std::vector<string> Converter::build_proc_vars(string s)
{
	string fname, var;
	size_t startf, endf, endvars, buf, buf1;
	std::vector<string> vars;
	startf = s.find(" ");
	endf = s.find("(");
	fname = s.substr(startf+1, endf-startf-1);
	////std::cout << fname << std::endl;
	endvars = s.find_last_of(")");
	////std::cout << "beg " << endf << " end " << endvars << std::endl;
	var = s.substr(endf+1, endvars-endf-1);
	////std::cout << var << std::endl;
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

//function Build func head out of complete string
string Converter::build_func_head(string s)
{
	string out = "double " + m_class_name + "::";
	string fname, var;
	size_t startf, endf, endvars, buf, buf1;
	std::vector<string> vars;
	startf = s.find(" ");
	endf = s.find("(");
	fname = s.substr(startf+1, endf-startf-1);
	////std::cout << fname << std::endl;
	endvars = s.find_last_of(")");
	////std::cout << "beg " << endf << " end " << endvars << std::endl;
	var = s.substr(endf+1, endvars-endf);
	// if there is also a unit del it
	//std::cout << "var: "<< var << std::endl;
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
			//std::cout << "buf2: " << buf2 << std::endl;
		}
		vars.push_back(buf2);
		buf = buf1+1;
	}
	//std::cout << "buf:  " << buf <<"bis: " << var.find("(", buf)-buf << std::endl;
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
	//if (counter > 1)
	// if any comment there we need to make on var lower
	if (s.find(":")!=s.npos)
		counter -=1;

		counter +=1;
	return counter;
}

//giving position of a letter
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

	////std::cout << "pos_letter does anything" << std::endl;
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

//checks if there is a next letter
bool Converter::pos_letterb(string s)
{
	vector<size_t> pos;
	//std::cout << "pos_letter does anything " << s << std::endl;
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


// ließt einzelne NModl-File ein
std::vector<string> Converter::Openfile(string filename)
{

	std::vector<string> zeilen;
	//try {

	const char* filenamechar = filename.c_str();
	ifstream ModFile(filenamechar, std::ios::out);

	string zeile;
		while (!ModFile.eof())          // Solange noch Daten vorliegen
	{
		////std::cout << "while" << std::endl;
	    getline(ModFile, zeile);        // Lese eine Zeile
	    /*if (zeile.find(":")!=zeile.npos)
	    {
	    	//std::cout << "replace happens" << std::endl;
	    	zeile.replace(zeile.find(":"), 1,";//");
	    }// replace : with // as comment*/
	    zeilen.push_back(zeile);    // Zeige sie auf dem Bildschirm
	}
	ModFile.close();    //}
	  //catch (std::ifstream::failure e) {
	    //std::cerr << "Exception opening/reading/closing file\n";
	  //}

	//std::cout << "einlesen erfolgreich" << std::endl;
	return zeilen;
}

//Defining all single Blocks without Block-Name Information
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

// Writing needed Block out of Nmodl-file
std::vector<string> Converter::GetBlock(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name)
{
	std::vector<string> BlockName;
	bool Verba = false;

	for (size_t i=0; i<Pairs.size(); i++)
	{
		////std::cout << "first" << Pairs.size() << std::endl;
		// first find of right name and second find that not outsurced
		if ((Zeilen[Pairs[i].first].find(name)!=Zeilen[Pairs[i].first].npos) && (Zeilen[Pairs[i].first].find(":"+name)==Zeilen[Pairs[i].first].npos))
		{
			////std::cout << Zeilen[Pairs[i].first] << std::endl;
			//überprüft wenn weitere klammern da sind, dass nur die letzte verwendet wird;

			while (Pairs[i].first == Pairs[i+1].first)
			{
				////std::cout << "in whhile by i: " << i << std::endl;
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
				////std::cout << counter << std::endl;
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

		// Part convertes all "^" to pow
					string left, right;
					string zeile, new_zeile, zeile_without;
					string orginal;
					zeile_without = Remove_all(Zeilen[j]);
					breaks = zeile_without.find("^");
					while (breaks!=zeile_without.npos)
					{
						//std::cout << "^ gefunden" << std::endl;
						// writes into left every part of one number
						left = "";

						//std::cout << "breaks: " << breaks << std::endl;
						for (size_t k=breaks-1; k>1; k--)
						{
							//sleep(10);
							//std::cout<< "k schleife with k: " << k << std::endl;
							//std::cout << zeile_without[k] << std::endl;
							if ((is_single_number(zeile_without[k])==true) || (is_single_letter(zeile_without[k])==true))
							{
								left = left + zeile_without[k];
							}
							else
								{
								//std::cout << "in else" << std::endl;
								k = 1; //break for for
								}
						}
						// left needs rewinded
						std::reverse(left.begin(), left.end());

						right = "";
						size_t clip_on, clip_off;
						clip_on = 0;
						clip_off = 0;
						bool change = false;
						bool clips = false;
						size_t length = zeile_without.size();
						for (size_t k=breaks; k<=length; k++)
						{
							if (zeile_without[breaks+1]=='(')
							{
								clips = true;
								if (zeile_without[k]=='(' || zeile_without[k]==')')
								{

									if (zeile_without[k]=='(')
										clip_on +=1;
									if (zeile_without[k]==')')
										clip_off +=1;

									// number of clips equal means end of pow second
									if (clip_on == clip_off)
									{
										right = zeile_without.substr(breaks+1, k-(breaks));
										k = zeile_without.size();
									}

								}
							}
							// if the problem have to be solved without clips
							if (clips == false)
							{
								//std::cout << "in clip false" << std::endl;
								// do only one times change length and k
								if (change == false)
								{
									k=k+1;
									length = length-1;
									change = true;
								}

								if ((is_single_number(zeile_without[k])==true) || (is_single_letter(zeile_without[k])==true))
								{
									right = right + zeile_without[k];
								}
								else
								{
									k = zeile_without.size();
								}
							}
						}



						//std::cout << right << std::endl;
						//std::cout << left << std::endl;
						////std::cout << "zeile: " << zeile << std::endl;
						new_zeile = zeile_without.replace(zeile_without.find(right), right.size(), "");
						//std::cout << "new_zeile1: " << new_zeile << std::endl;
						new_zeile = new_zeile.replace(new_zeile.find(left), left.size(), "");
						//std::cout << "new_zeile2: " << new_zeile << std::endl;
						new_zeile = new_zeile.replace(new_zeile.find("^"), 1, " pow(" + left + " , " + right + ")");
						//std::cout << "new_zeile3: " << new_zeile << std::endl;
						////std::cout << "while schleife geht net mit Zeilen[j] = " << Zeilen[j] << std::endl;

						Zeilen[j] = new_zeile;
						breaks = Zeilen[j].find("^");
						////std::cout << "breaks " << breaks << std::endl;
						//end of converting action
					}
					// writting formatted Zeile in Blockname
					BlockName.push_back(Zeilen[j]);
				}
				if (Zeilen[j].find("ENDVERBATIM")!=Zeilen[j].npos)
					Verba=false;

			}

		}
	}
	////std::cout << "before block" << std::endl;
	////std::cout << BlockName.size() << std::endl;
	return BlockName;
}

// removes all comments out of a string giving back new string
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



// Function which removes all " " out of a string giving back new string
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

// Function which removes all " " out of a vector<string> giving back new vector<string>
vector<string> Converter::Remove_all(vector<string> erg)
{
	  size_t HFile_del;
	  //delete all not needed symbols out of a list
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


		  //std::cout<< "Hfile_vars: " << erg[i] << std::endl;
	  }
	  return erg;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Function which writes H File and Cpp File
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Converter::WriteStart(string filename, std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, bool force)
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
	  //std::cout << "write begins" << std::endl;
	  mycppfile << "#include \"" + filename +".h\"	\n";
	  mycppfile << "#include \"lib_grid/lg_base.h\" \n";
	  mycppfile << "#include \"lib_disc/spatial_disc/elem_disc/elem_disc_interface.h\" \n";
	  mycppfile << "#include \"lib_disc/function_spaces/grid_function.h\" \n";
	  mycppfile << "#include \"lib_disc/function_spaces/local_transfer_interface.h\" \n";
	  mycppfile << "#include <cmath> \n";

	  mycppfile << "namespace ug { \n \n \n";

	  std::vector<string> mod_funcs_names;
	  // all function have to be read in first so later need more attention more

	  std::vector<string> FUNCTION = GetBlock(Pairs, Zeilen, "FUNCTION");

	  string without[] = {"if", "else"};
	  size_t funcS, funcE;
	  bool second_time;

	  vector<string> func_hfile;

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
				  func_hfile.push_back(func);

				  // writing needed templates in function!
				  string left = func.substr(0, func.find("::"));
				  string right = func.substr(func.find("::")+2, func.npos-(func.find("::")+2));
				  func = left + "<TDomain>::" + right;
				  mycppfile << "template<typename TDomain> \n";
				  // Todo template params in head benötigt
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
								 && (string(FUNCTION[i] + "; \n").find("};")==string(FUNCTION[i] + "; \n").npos) && (FUNCTION[i]!="")
								 && Remove_all(FUNCTION[i])!="}")
							  {
								  mycppfile << FUNCTION[i] + "; \n";
							 	  ////std::cout << FUNCTION[i];
							  }
							  else
							  {
							 	  ////std::cout << FUNCTION[i];
							 	  mycppfile << FUNCTION[i] + "\n";
							  }

						  }
					  }
				  }
				  //mycppfile << "} \n \n";

			  }
		  }

	  	  mycppfile << "\n \n \n";
	  }

	  // Proceduren wie Funktionen
	  //std::cout << "All Functions added" << std::endl;

	  mycppfile << "// adding function which always inits_attachments \n";
	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void "+filename+"<TDomain>::vm_disc_available()  \n";
	  mycppfile << "{  \n";
	  mycppfile << "	init_attachments();  \n";
	  mycppfile << "}  \n \n \n \n";

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

	  myhfile << "#include \"VM_Disc.h\" \n \n";

	  myhfile << "#include <vector> \n";
	  myhfile << "#include <stdio.h> \n";
	  myhfile << "#include \"bindings/lua/lua_user_data.h\" \n";

	  myhfile << "namespace ug { \n \n";
	  myhfile << "// forward declaration \n";
	  myhfile << "template <typename TDomain> \n";
	  myhfile << "class VMDisc; \n \n";

	  myhfile << "template <typename TDomain> \n";
	  myhfile << "class " + filename + "\n";
	  myhfile << "    : public IChannel<TDomain> \n";
	  myhfile << "{ \n";
	  myhfile << "    public: \n ";
	  myhfile << "using IChannel <TDomain>::m_pVMDisc; \n \n";


	  // finde die spezifischen Parameter
	  // out of neuron only GLOBALS needed
	  std::vector<string> NEURON = GetBlock(Pairs, Zeilen, "NEURON");
	  size_t Global;
	  std::vector<string> Global_vars;
	  // allreadings and all writings needed later!
	  std::vector<string> writes_Ions = Find_all_write(NEURON);
	  std::vector<string> read_Ions = Find_all_read(NEURON);

	  /*
	  //std::cout << "Find_all_write and Find all read testing area: " << std::endl;
	  for (size_t i=0; i<writes_Ions.size(); i++)
	  {
		  //std::cout << "writer: " << writes_Ions[i] << std::endl;
	  }

	  for (size_t i=0; i<read_Ions.size(); i++)
	  {
		  //std::cout << "reader: " << read_Ions[i] << std::endl;
	  }*/


	  // get all vars which nernst equalis is used later on
	  std::vector<string> All_Eqs = Find_all_Eqs(NEURON);

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
		  for (size_t i = 0; i<Global_vars.size()-1; i++)
		  {
			  HFile_added_Vars.push_back(Global_vars[i]);
			  //myhfile << "double " + Global_vars[i] + "; \n";
		  }

	  }
	  myhfile << "\n \n \n" ;

	  //std::cout << "Globals added!" << std::endl;



//////////////////////////////////////////////////////////////////////////////////////////////
// Writing out all hard-coded Parameters and change in needed Unit
//////////////////////////////////////////////////////////////////////////////////////////////


	  // getting out of Params hard coded Values
	  std::vector<string> PARAMETER = GetBlock(Pairs, Zeilen, "PARAMETER");
	  size_t end, endUnit;
	  string var, UnitS;
	  double Unit;

	  std::vector<string> Params_Unit;
	  // Vector with all setted Params
	  // this list is needed for example if ena has a hard coded value
	  std::vector<string> Parameter;

	  // setting all parameters and converting units
	  //std::cout << "Setting Parameters and converting in needed Units" << std::endl;
	  if (PARAMETER.size() > 0)
	  {
		  for (size_t i=0; i<PARAMETER.size(); i++)
		  {
			  std::stringstream Units;
			  Unit = 0;
			  if (PARAMETER[i].find("=")!=PARAMETER[i].npos)
			  {
				  end = PARAMETER[i].find("(");
				  endUnit = PARAMETER[i].find(")");

				  // no clips but comment add the end
				  if (endUnit==PARAMETER[i].npos && end==PARAMETER[i].npos)
				  {
					  if (PARAMETER[i].find(":")!=PARAMETER[i].npos)
					  {
						  PARAMETER[i] = PARAMETER[i].substr(0, PARAMETER[i].find(":")-1);
					  }
					  Unit = 1;
				  }

				  var = PARAMETER[i].substr(0, end);
				  if (end!= PARAMETER[i].npos && endUnit!= PARAMETER[i].npos)
				  {
					  Unit = Unit_Conv_All(PARAMETER[i].substr(end, endUnit));
				  }
				  Units << Unit;
				  Params_Unit.push_back(var + "*" + Units.str());
				  HFile_added_Vars.push_back(var.substr(0, var.find("=")));
				  Parameter.push_back(var.substr(0, var.find("=")));
			  }
		  }
	  }
	  //myhfile << "std::vector<number> m_diff; \n";
	  myhfile << "\n";


	  //std::cout << "All Parameters set" << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////
// Writing all
//////////////////////////////////////////////////////////////////////////////////////////////

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
			  //std::cout << IonEnd << ", " << Ion << std::endl;
			  IonS = NEURON[i].substr(Ion+7, (IonEnd-(Ion+7)));
			  //std::cout << "Subion: "<<IonS << std::endl;
			  ListIons.push_back(IonS);
			  //writes var for outer concentrations
			  //myhfile << "number " + IonS + "_out; \n";
			  HFile_added_Vars.push_back(IonS + "_out");
		  }

	  }


	  // delete all unneeded style-formats
	  Remove_all(HFile_added_Vars);
	  // make some copys of HFile_added_Vars for later use
	  //because in different functions different locals are used
	  vector<string> HFile_added_Vars_org = HFile_added_Vars;
	  vector<string> HFile_added_Vars_Deriv = HFile_added_Vars;
	  vector<string> HFile_added_Vars_Init = HFile_added_Vars;


	  myhfile << "\n"; myhfile << "\n";

	  myhfile << "/// @copydoc IChannel<TDomain>::IChannel(cont char*) \n";
	  myhfile << filename + "(const char* functions, const char* subsets) \n";
	  myhfile << "try : IChannel<TDomain>(functions, subsets), \n";
	  myhfile << "m_R(8.314), m_T(293.0), m_F(96485.0), \n";
	  /// adding all standard vals of Parameters
	  size_t gleich;
	  string add;
	  for (size_t i=0; i<Params_Unit.size(); i++)
	  {
		  gleich = Params_Unit[i].find("=");
		  add = Params_Unit[i].substr(0, gleich) + "(" + Params_Unit[i].substr(gleich+1) + ")";
		  if (i+1<Params_Unit.size())
			  add = add + ", \n";
		  else
			  add = add + " {} \n";
		  myhfile << add;
	  }
	  myhfile << "UG_CATCH_THROW(\"Error in "+ filename + " initializer list. \") \n \n \n";

	  myhfile << "/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) \n";
	  myhfile << filename + "(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) \n";
	  myhfile << "try : IChannel<TDomain>(functions, subsets), \n";
	  myhfile << "m_R(8.314), m_T(293.0), m_F(96485.0), \n";
	  /// adding all standard vals of Parameters
	  size_t gleich1;
	  string add1;
	  for (size_t i=0; i<Params_Unit.size(); i++)
	  {
	  	gleich1 = Params_Unit[i].find("=");
	  	add1 = Params_Unit[i].substr(0, gleich1) + "(" + Params_Unit[i].substr(gleich1+1) + ")";
	  	if (i+1<Params_Unit.size())
	  		add1 = add1 + ", \n";
	  	else
	  		add1 = add1 + " {} \n";
	  	myhfile << add1;
	  }
	  myhfile << "UG_CATCH_THROW(\"Error in "+ filename + " initializer list. \") \n";

	  myhfile << "/// destructor \n \n";
	  myhfile << "virtual ~"+ filename + "() {}; \n";

	  string name;

	  // adding function out of function list
	  for (size_t i=0; i<func_hfile.size() ; i++)
	  {
		  name = func_hfile[i].substr(func_hfile[i].find(filename), func_hfile[i].find("::")-func_hfile[i].find(filename)+2);
		  //std::cout << "name: " << name << std::endl;
		  func_hfile[i].replace(func_hfile[i].find(name), name.size(), "");
		  myhfile << func_hfile[i] + "; \n";
	  }


	  /// adding functions of IChannel
	  myhfile << "/// create attachments and accessors \n";
	  myhfile << "void init_attachments(); \n";
	  myhfile << "// inherited from IChannel \n \n";
	  myhfile << "virtual void init(const LocalVector& u, Edge* e); \n";
	  myhfile << "virtual void update_gating(number newTime, const LocalVector& u, Edge* e); \n";
	  myhfile << "virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); \n";
	  myhfile << "virtual void vm_disc_available(); \n";
	  myhfile << "\n \n";
	  myhfile << "protected: \n";

	  //myhfile << "void register_func(); \n \n \n";
	  myhfile << "private: \n \n";
	  // Neuron-lines with use ion







	  std::vector<string> STATE = GetBlock(Pairs, Zeilen, "STATE");

	  // needed var for Gating params having same value of used Ion
	  std::vector<string> SGating;


	  // creating init Attachments
	  mycppfile.open (filenamecpp, std::ios::app);
	  mycppfile << " // creating Method for attachments \n";
	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::init_attachments() \n";
	  mycppfile << "{ \n";
	  mycppfile << "// inits temperatur from kalvin to celsius and some other typical neuron values\n";
	  mycppfile << "m_pVMDisc->celsius = m_T - 273; \n \n \n";
	  mycppfile << "SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); \n";

////////////////////////////////////////////////////////////////////////////////////
//// Saving all States in one Var also writting first states into cpp and h file
/////////////////////////////////////////////////////////////////////////////////////
	  /*size_t Ns_Current, Ns_CurrentEnd;
	  string Ns_CurrentS;*/

	  string addState;

	  // var in which all states and Non spec currents are saved
	  std::vector<string> State_vars;
	  size_t state, stateend, komm_beg;
	  size_t varSt = 0;
	  bool clip = false;
	 // size_t tabs;

	  if (STATE.size()>0)
	  {
		  for (size_t i=0; i<STATE.size(); i++)
		  {
			  state = pos_letter(STATE[i]);
			  if (state!=1000)
			  {
				  if (STATE[i].find("(")==STATE[i].npos)
				  {
					  varSt = number_(state, STATE[i]);
				  }
				  else
				  {
					  varSt = 1;
					  clip=true;
				  }
			  }
			  ////std::cout << "after doppelt if" << std::endl;
			  for (size_t j=0; j<varSt; j++)
			  {
				  ////std::cout << "in state " << std::endl;
				  komm_beg = STATE[i].find(":");
				  stateend = STATE[i].find(" ", state);
				  if (clip==true)
					  stateend = STATE[i].find("(", state);
				  clip = false;
				  ////std::cout << "ab: " << state << " anzahl: " << stateend-state << std::endl;
				  if (stateend==komm_beg+1)
				  {
					  stateend -= 1;
				  }
				  if ((stateend<=komm_beg) && (state<komm_beg))
				  {
					  if (STATE[i].find("}")!=STATE[i].npos)
		  			  {

						  if (stateend-state!=STATE[i].npos-state)
						  {
							  addState = Remove_all(STATE[i].substr(state, stateend-state));
						  }
						  else
						  {
							  ////std::cout << "not added!!" << std::endl;
							  addState = "";
						  }
		  			  }
					  else
		  			  {
						  ////std::cout << "vergleich: " << stateend << " - "<<komm_beg << std::endl;
						  addState = Remove_all(STATE[i].substr(state, stateend-state));
		  			  }

		  			  //writting in hfile
		  			  if ((addState!="}") && (addState!="") && (addState!="\n") &&(addState!="\n}") && (addState!="}\n"))
		  			  {
		  				// If State is already in Ion List we need another name for State
		  				// doing this by adding "S"
		  				  for (size_t k=0; k<ListIons.size(); k++)
		  				  {
		  					  if (Remove_all(addState)==Remove_all(ListIons[k]))
		  					  {
		  						  addState = addState+"S";
		  						  SGating.push_back(addState);
		  					  }
		  				  }
		  				  State_vars.push_back(addState);
		  				  myhfile << "ADouble " + addState + "Gate; \n";
		  				  myhfile << "Grid::AttachmentAccessor<Vertex, ADouble> aa" + addState + "Gate; \n";

		  				  //writting of attachment inits for cpp file
		  				  mycppfile << "if (spGrid->has_vertex_attachment(this->" + addState + "Gate)) \n";
		  				  mycppfile << "UG_THROW(\"Attachment necessary (" + addState + "Gate) for " +filename+" channel dynamics \"\n";
		  				  mycppfile << "\"could not be made, since it already exists.\"); \n";
		  				  mycppfile << "spGrid->attach_to_vertices(this->" + addState + "Gate); \n";
		  				  mycppfile << "this->aa" + addState + "Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->" + addState + "Gate); \n \n";

		  			  }

		  		  }
				  state = stateend+1;
		  		  //////////////////
			  }
			  varSt = 0;


		  }

	  }
	  mycppfile << "} \n \n \n \n";



	  // Method for initialization of States with values read out of Initial
	  mycppfile << " // Init Method for using gatings \n";
	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::init(const LocalVector& u, Edge* edge) \n";
	  mycppfile << "{ \n";
	  mycppfile << "//get celsius \n";
	  mycppfile << "number celsius = m_pVMDisc->celsius; \n";
	  mycppfile << "// make preparing vor getting values of every edge \n";
	  mycppfile << "typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; \n";
	  mycppfile << "vrt_list vl; \n";
	  mycppfile << "m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); \n \n \n";
	  mycppfile << "//over all edges \n";
	  mycppfile << "for (size_t l = 0; l< vl.size(); l++) \n";
	  mycppfile << "{ \n";
	  mycppfile << "\t Vertex* vrt = vl[l]; \n \n \n";
	  mycppfile << "number v = u(m_pVMDisc->_v_, l); \n";

	  // TODO we only need read ions in right context
	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "number " + ListIons[i] +" = u(m_pVMDisc->_"+ ListIons[i] +"_, l); \n";
	  }
	  mycppfile << "\n \n";


	  std::cout << "Start working on Initial-Block" << std::endl;

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
	  ////std::cout << "beg: INITIAL!" << std::endl;

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
					  		  // add only if they are in gating or ion list
					  		  string GatingName = In_NeuronUse_List(Pairs, Zeilen, Proc_funcs[j][k]);
					  		  if (GatingName!="")
					  		  {
					  			  mycppfile << "double " + Remove_all(Proc_funcs[j][k]) + " = aa"+ GatingName + "Gate[vrt]; \n \n";
					  			  ////std::cout << "else for: " << std::endl;
					  		  }
					  	  }
					  ////std::cout << Proc_funcs[j][k] << std::endl;
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
									  HFile_added_Vars_Init.push_back(locals[k].substr(0, locals_sep));
									  locals[k].replace(0, locals_sep+1, "");
									  locals_sep = locals[k].find(",");
								  }
								  HFile_added_Vars_Init.push_back(Remove_all(locals[k].substr(0, locals_sep)));
							  }


						  }
					  }
					  ////std::cout << "before proc vals" << std::endl;
					  Proc_vals = writer_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);
					  ////std::cout << Proc_vals.size() << std::endl;

					  bool HFile_added;
					  for (size_t k=0; k<Proc_vals.size(); k++)
					  {
						  HFile_added = false;
						  if ((Proc_vals[k]!="") && (Proc_vals[k]!="}") && (Proc_vals[k].find("LOCAL")==Proc_vals[k].npos)
							   && Proc_vals[k]!=" " && Proc_vals[k]!="\t" && Proc_vals[k]!="        ")
						  {
  							  for (size_t l=0; l<HFile_added_Vars_Init.size(); l++)
  							  {
  								  ////std::cout << "Vergleich " << Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")-1)) << " - " << Remove_all(HFile_added_Vars[l]) << std::endl;
  								  if (Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")))==Remove_all(HFile_added_Vars_Init[l]))
  									  HFile_added = true;
  							  }

  							  string com_buf;

  							  // if it added in hfile write without double else with!
  							  if (HFile_added == true)
  							  {
  								  //writes comments
  								  if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
  									  Proc_vals[k].replace(Proc_vals[k].find(":"), 1, "//");
  								  ////std::cout << Proc_vals[k] << std::endl;
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
  								////std::cout << "com_buf: " << com_buf << std::endl;

  								if (com_buf.find(":")==0)
  								{
  									////std::cout << Proc_vals[k] << std::endl;
  									mycppfile << Proc_vals[k] + "\n";
  								} else
  								{
  									////std::cout << Proc_vals[k] << std::endl;
  									if (Remove_all(Proc_vals[k])!="")
  										mycppfile << "double " + Proc_vals[k] + "; \n";

  									// perhaps needed and also needed in deriv
  									// if var used a second time it does not need double at beginning
  									//HFile_added_Vars_Init.push_back(Remove_all(Proc_vals))

  								}
  							  }

						  }
					  }


				  }
			  }
		  }		  // write locals vals as double
	  }

	  ////std::cout<< "writing end" << std::endl;
	  for (size_t i=Stats_beg; i<INITIAL.size(); i++)
	  {
		  // writings at the end
		  if (i >= Stats_beg)
		  {
			  size_t gleich, vorgleich;
			  if ((INITIAL[i]!="") && (INITIAL[i]!="}"))
			  {
				  ////std::cout << "gleich Search" << std::endl;
				  vorgleich = pos_letter(INITIAL[i]);
				  gleich = INITIAL[i].find("=");

				  if (INITIAL[i].find("=")!=INITIAL[i].npos)
				  {
					  ////std::cout << "ab: " << vorgleich << "bis: " << gleich << std::endl;
					  // needed if they are in gatting list use as gatting else write without gatting
					  bool written = false;
					  for (size_t j=0; j< State_vars.size(); j++)
					  {
						  if (Remove_all(INITIAL[i].substr(vorgleich , gleich-vorgleich))==Remove_all(State_vars[j]))
						  {
							  mycppfile << "aa" + Remove_all(INITIAL[i].substr(vorgleich, gleich-vorgleich)) + "Gate[vrt] " + INITIAL[i].substr(gleich) + "; \n";
							  written = true;
						  }
					  }
					  {
						  if (written==false)
						  {
							  if (INITIAL[i].find("=")!=INITIAL[i].npos)
							  {
								  string left = INITIAL[i].substr(0, INITIAL[i].find("=")-1);
								  string right = INITIAL[i].substr(INITIAL[i].find("=")+1, INITIAL[i].npos-(INITIAL[i].find("=")+1));
								  mycppfile << ProofSGatingInit(left, SGating) + " = " + right + "; \n";
							  }

							  else
							  {
							  mycppfile << INITIAL[i] + "; \n";
							  }
						  }
					  }
				  } //else
			    	  //vorgleich = beg_count(INITIAL[i]);
			  }
		  }
	  }

	  mycppfile << "}  \n";
	  mycppfile << "}  \n \n \n \n";

	  std::cout << "Initial block written" << std::endl;


////////////////////////////////////////////////////////////////////////////////////////////////
// Starts writting update_gating function
///////////////////////////////////////////////////////////////////////////////////////////////



	  std::cout << "Start writting update_gating function" << std::endl;

	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) \n";
	  mycppfile << "{ \n";
	  mycppfile << "number celsius = m_pVMDisc->celsius; \n \n";

	  mycppfile << "// make preparing vor getting values of every edge \n";
	  mycppfile << "typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; \n";
	  mycppfile << "vrt_list vl; \n";
	  mycppfile << "m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); \n \n \n";
	  mycppfile << "//over all edges \n";
	  mycppfile << "for (size_t l = 0; l< vl.size(); l++) \n";
	  mycppfile << "{ \n";
	  mycppfile << "\t Vertex* vrt = vl[l]; \n \n \n";
	  mycppfile << "number dt = newTime - m_pVMDisc->m_aaTime[vrt]; \n";
	  mycppfile << "number v = u(m_pVMDisc->_v_, l); \n";

	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "number " + ListIons[i] +" = u(m_pVMDisc->_"+ ListIons[i] +"_, l); \n";
	  }
	  mycppfile << "\n \n";

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
			    //std::cout << "Procedure: " << PROCEDURE[i] << std::endl;
			  }
		  for (size_t i=0; i<Proc_funcs[0].size(); i++)
			  {
			  	  for (size_t j=0; j<Proc_funcs.size(); j++)
			  	  {
			    //std::cout << "Proc_funcs: " << Proc_funcs[j][i] << std::endl;
			  	  }
			  }

		  if (Der_funcs.size() >= 1)
		  {
			  for (size_t i=0; i<Der_funcs[0].size(); i++)
			  	  {
			  	  	  for (size_t j=0; j<Der_funcs.size(); j++)
			  	  	  {
			  	  		  //std::cout << "Proc_funcs: " << Der_funcs[j][i] << std::endl;
			  	  	  }
			  	  }
		  }
	  }




	  // States ever needed for deriv writting out of State_vars State_vars;
	  for (size_t i=0; i<State_vars.size(); i++)
	  {
		  mycppfile << "double " + State_vars[i] + " = aa" + State_vars[i] + "Gate[vrt]; \n";
		  HFile_added_Vars_Deriv.push_back(Remove_all(State_vars[i]));
	  }



	  //Adding v because v is very often needed and always accesible
	  //mycppfile << "double v  = aavGate[*iter]; \n";

	  mycppfile << "\n \n \n";

	  std::cout << "Working on Derivative-Block starts" << std::endl;
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
	  				  ////std::cout << "Vergleich: " << DERIVATIVE[i].substr(begin, DERIVATIVE[i].find("(")-begin)<< " - " << Proc_funcs[j][0] << std::endl;
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
						  		  // add only if they are in gating or ion list
						  		  string GatingName = In_NeuronUse_List(Pairs, Zeilen, Proc_funcs[j][k]);
						  		  if (GatingName!="")
						  		  {
						  			  mycppfile << "double " + Proc_funcs[j][k] + " = aa"+ GatingName + "Gate[vrt]; \n \n";
						  			  ////std::cout << "else for: " << std::endl;
						  		  }
	  					  	  }
	  					  ////std::cout << Proc_funcs[j][k] << std::endl;
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
										  HFile_added_Vars_Deriv.push_back(Remove_all(locals[k].substr(0, locals_sep)));
										  locals[k].replace(0, locals_sep+1, "");
										  locals_sep = locals[k].find(",");
									  }
									  HFile_added_Vars_Deriv.push_back(Remove_all(locals[k].substr(0, locals_sep)));
								  }


							  }
	  					  }
	  					  //std::cout << "before proc vals" << std::endl;
	  					  Proc_vals = writer_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);
	  					  //std::cout << Proc_vals.size() << std::endl;


	  					  bool HFile_added;
	  					  for (size_t k=0; k<Proc_vals.size(); k++)
	  					  {
							  HFile_added = false;
							  if ((Proc_vals[k]!="") && (Proc_vals[k]!="}") && (Proc_vals[k].find("LOCAL")==Proc_vals[k].npos)
								   && Proc_vals[k]!=" " && Proc_vals[k]!="\t" && Proc_vals[k]!="        ")
							  {
	  							  for (size_t l=0; l<HFile_added_Vars_Deriv.size(); l++)
	  							  {
	  								  //std::cout << "Vergleich " << Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")-1)) << " - " << Remove_all(HFile_added_Vars_Deriv[l]) << std::endl;
	  								  if (Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")))==Remove_all(HFile_added_Vars_Deriv[l]))
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
	  									//std::cout<< Proc_vals[k] << std::endl;
	  									if (Remove_all(Proc_vals[k])!="")
	  										mycppfile << "double " + Proc_vals[k] + "; \n";
	  								}
	  							  }

							  }
	  					  }


	  				  }
	  			  }
	  		  }		  // write locals vals as double
	  	  }

	  	  //std::cout<< "writing end" << std::endl;
	  	  for (size_t i=Stats_beg; i<DERIVATIVE.size(); i++)
	  	  {
	  		  // writings at the end
	  		  if (i >= Stats_beg)
	  		  {
				  // Writting derivfuncs
				  //std::cout << "before deriv" << std::endl;
				  if ((i> Stats_beg) && (i < DERIVATIVE.size()-1))
				  {
					if (DERIVATIVE[i]!="")
					{
						if ((DERIVATIVE[i].find("\'")!=DERIVATIVE[i].npos) && (DERIVATIVE[i].find("=")!=DERIVATIVE[i].npos))
						{
							DERIVATIVE[i].replace(DERIVATIVE[i].find("\'"), 1, "");
							DERIVATIVE[i].replace(DERIVATIVE[i].find("="), 1, "+=");
							DERIVATIVE[i] = DERIVATIVE[i] + "*dt";
							DERIVATIVE[i] = ProofSGatingGating(DERIVATIVE[i], SGating);
							// if that var is a Gating and Ion at same time
						}

						bool used = false;
						//std::cout << "before added Vars" << std::endl;
						for (size_t j = 0; j<HFile_added_Vars_Deriv.size(); j++)
						{
							if (DERIVATIVE[i].find("=")!=DERIVATIVE[i].npos)
							{
								//std::cout << "in derivv find =" << std::endl;
								if (Remove_all(DERIVATIVE[i].substr(0, DERIVATIVE[i].find("=")-1))==Remove_all(HFile_added_Vars_Deriv[j]))
									used = true;
							}
						}
						if (DERIVATIVE[i].find("if")!=DERIVATIVE[i].npos)
							used = true;


						if (used==false)
						{
							//std::cout << "in used false" << std::endl;
							if (Remove_all(DERIVATIVE[i])!="")
								mycppfile << "double " + DERIVATIVE[i] + "; \n";
						}
						else
						{
							mycppfile << DERIVATIVE[i] + "; \n";
						}
					}
				  }
	  		  }
	  	  }

	  mycppfile << "\n \n \n";



	  std::cout << "Derivate-Block is written" << std::endl;


	  std::cout << "Start working on kinetic functions" << std::endl;

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
			  //std::cout << "first needed proc: " << needed_proc << std::endl;
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

							  //std::cout << "left: " << left <<" right: "<< right << "vars: " << vars << std::endl;
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


		  if (needed_proc.find("BREAKPOINT")==needed_proc.npos)
		  {
			  //getting procedure for needed vars
			  size_t clip = needed_proc.find(")");
			  if (clip!=needed_proc.npos)
				  needed_proc.replace(clip, 1, "");
			  Ausgabe2 = GetProcEqualString(Pairs, Zeilen, needed_proc);
			  for (size_t i=0; i<Ausgabe2.size(); i++)
			  {
				  if (Ausgabe2[i].find("PROCEDURE")==Ausgabe2[i].npos)
					  mycppfile << Ausgabe2[i] + "; \n";
			  }
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

		  std::cout<< "Working on Kinetic block finished" << std::endl;

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
					  //std::cout << "PROC: " << Proc << std::endl;
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
			  mycppfile << "aa" + State_vars[i] + "Gate[vrt] = " + State_vars[i] + "; \n";
		  }
	  }



	  mycppfile << " \n \n \n";
	  mycppfile << "} \n";
	  mycppfile << "} \n";
	  mycppfile << " \n \n \n";



	  std::cout << "update_gating Function is written" << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// write head of ionic flux
//////////////////////////////////////////////////////////////////////////////////////////////////////////

	  std::cout << "Start writting Ionic_current function"<< std::endl;


	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) \n";
	  mycppfile << "{ \n \n";
	  // all gates

	  if (State_vars.size()>0)
	  {
		  // States ever need new values
		  	  for (size_t i=0; i<State_vars.size(); i++)
		  	  {
		  		  mycppfile << "number " + State_vars[i] + " = aa" + State_vars[i] + "Gate[ver]; \n";
		  	  }
	  }

	  if (ListIons.size()>0)
	  {
		  // all ions
		  for (size_t i=0; i<ListIons.size(); i++)
		  {
		 	  mycppfile << "number " + ListIons[i] + " = vrt_values[VMDisc<TDomain>::_" + ListIons[i] + "_]; \n";
		  }
	  }

	  // v needed every time
	  mycppfile << "number v =  vrt_values[VMDisc<TDomain>::_v_]; \n";
	  mycppfile << " \n";
	  mycppfile << " \n";

	  //writes Nernst Eqs
	  std::vector<string> eqs = equali(Pairs, Zeilen);

	  //var vor latest use
	  std::vector<string> comp_outs;
	  string comp_outs_left, comp_outs_right;

	  std::cout << "Start preparing Nernst-Equations" << std::endl;

	  if (eqs.size()>0)
	  {
		  for (size_t i=0; i<eqs.size(); i++)
		  {
			  bool is_comp_out = false;
			  // TODO write for schleife over all ions!
			  if (eqs[i]!="1" && eqs[i]!="0" && Remove_all(eqs[i])!="ca" && Remove_all(eqs[i])!="na" && Remove_all(eqs[i])!="k")
			  {
				  //std::cout << eqs[i] + " \n" << std::endl;
				  // other handle if it is ica, ik or ina which means outflux
				  if (eqs[i].find("=")!=eqs[i].npos)
				  {
					  comp_outs_right = eqs[i].substr(eqs[i].find("=")+1, eqs[i].npos - (eqs[i].find("=")+1));
					  comp_outs_left = eqs[i].substr(0, eqs[i].find("=")-1);
					  if ((comp_outs_left.find("ica")!=comp_outs_left.npos)
					  ||  (comp_outs_left.find("ina")!=comp_outs_left.npos)
					  ||  (comp_outs_left.find("ik")!= comp_outs_left.npos))
					  {
						  comp_outs_left = comp_outs_left.replace(comp_outs_left.find("number"), 6,"");
						  comp_outs_left = Remove_all(comp_outs_left);
						  comp_outs_right = ProofSGatingGating(comp_outs_right, SGating);
						  comp_outs_right = comp_outs_right.replace(comp_outs_right.find(";"), 1, "");
						  comp_outs_right = comp_outs_right.replace(comp_outs_right.find("+="), 2, "=");
						  comp_outs.push_back(comp_outs_left +  comp_outs_right);
						  is_comp_out = true;
					  }
				  }

				  bool rev_pot_is_setted = false;
				  // if values in eqs setted before in NModl file do not calculate Nernstequatation of this param
				  for (size_t j=0; j<Parameter.size(); j++)
				  {

					  string verg = eqs[i].substr(0, eqs[i].find("=")-1);
					  verg = verg.replace(verg.find("number"), 6, "");
					  //std::cout << "Vergleich: " << Remove_all(verg) << " vs " << Remove_all(Parameter[j]) << std::endl;
					  if (Remove_all(verg) == Remove_all(Parameter[j]))
						  rev_pot_is_setted = true;
				  }

				  if (is_comp_out == false && rev_pot_is_setted == false)
					  mycppfile << eqs[i] + " \n";
			  }
		  }
	  }

	  mycppfile << " \n \n";

	  std::cout << "Nernst-Equations prepared" << std::endl;


	  bool in = true;
	  bool stand_BP = true;

	  size_t beg;
	  std::vector<string> fluxes;
	  std::vector<string> B_vars;


	  std::cout << "Working on Breakpoint-Block" << std::endl;

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
		 			  //std::cout << "adding flux" << std::endl;
		 			  beg = pos_letter(BREAKPOINT[i]);
		 			  if (beg!=1000)
		 				  fluxes.push_back(BREAKPOINT[i].substr(beg, BREAKPOINT[i].npos-beg));
		 		  }

		 	  }
		 	  std::cout << "Breakpoint-Block Infos written in fluxes" << std::endl;

		 	  for (size_t i=0; i<fluxes.size(); i++)
		 	  {
		 		  // make list with left handed vars
		 		  B_vars.push_back(fluxes[i].substr(0, fluxes[i].find("=")-1));
		 		  //std::cout << (B_vars[i]) << std::endl;
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
		 					  //hier??
		 					  //std::cout << "number " + fluxes[i] + "; \n" << std::endl;
		 				  	  mycppfile << "number " + fluxes[i] + "; \n";
		 				  	  fluxes[i] = "";
		 				  }
		 			  }

		 		  }
		 	  }

		 	  vector<bool> Ion_outside;
		 	  vector<string> outs;
		 	  vector<string> outsHelp;
		 	  // first every time v
		 	  outs.push_back("");

		 	  int vm_flux_count = 0;

		 	  mycppfile << "\n \n \n";

			  //// adding writen comp_outs to file
			  for (size_t i = 0; i<comp_outs.size(); i++)
			  {
				  fluxes.push_back(comp_outs[i]);
			  }

		 	  for (size_t i=0; i<fluxes.size(); i++)
		 	  {
		 		  //std::cout << "fluxes[i]: " << fluxes[i] << std::endl;

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
		 			  //std::cout << "List Ions" << std::endl;
		 			  //std::cout << "i"+ListIons[j] << std::endl;
		 			  string change;
		 			  // Here we are ousing infos out of equality
		 			  for (size_t k=0; k<eqs.size(); k++)
		 			  {
		 				  if (eqs[k]=="0")
		 				  {
		 					  size_t eqs_leer = eqs[k+1].find(" ");
		 					  while (eqs_leer!=eqs[k+1].npos)
		 					  {
		 						  eqs[k+1].replace(eqs_leer, 1, "");
		 						  eqs_leer = eqs[k+1].find(" ");
		 						  //std::cout << "eqs k+1: " << eqs[k+1] << std::endl;
		 					  }
		 					  if (eqs[k+1]==ListIons[j])
		 					  {
		 						  //setting Ion_outside true to prevent output later as push_back from ionic_current
		 						  //cause ion is only working on outside
		 						  Ion_outside[j]=true;

		 						  size_t write, komma;
		 						  string left, right;
		 						  // create out of NEURON another output depending on WRITE Part!
		 						  ////std::cout << "changes will happen on outer concentration" << std::endl;
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
		 									 //std::cout << eqs[k+1] + fluxes[i].substr(fluxes[i].find("="), fluxes[i].npos-fluxes[i].find("=")) + "; \n" << std::endl;
		 									 mycppfile << eqs[k+1] + fluxes[i].substr(fluxes[i].find("="), fluxes[i].npos-fluxes[i].find("=")) + "; \n";
		 								 }

		 								 if (fluxes[i].find(right)!=fluxes[i].npos)
		 								 {
		 									 mycppfile << "number " + eqs[k+1] + "o = m_pVMDisc->" + eqs[k+1] + "_out; \n";
		 									 mycppfile << "\n \n \n";
		 								 }

		 							  }
		 						  }

		 					  }
		 				  }


		 			  }
		 			  // adding always in outs
		 			  if ((fluxes[i].find("i"+ListIons[j])!=fluxes[i].npos))
		 			  {
		 				  outs.push_back(fluxes[i].substr(fluxes[i].find("=")+1, fluxes[i].npos -fluxes[i].find("=")+1));
		 				  outsHelp.push_back(ListIons[j]);
		 			  }
		 		  }
		 	  }





		 	  mycppfile << "outCurrentValues.push_back(" + outs[0] + "); \n";
		 	  for (size_t i=1; i<outs.size(); i++)
		 	  {
		 		  bool optional = false;
		 		  for (size_t j=0; j<writes_Ions.size(); j++)
		 		  {
		 			 // additional functionality added
		 			 //std::cout << Remove_all(writes_Ions[j]) << " vs " << Remove_all(outsHelp[i-1]+"i") << std::endl;
		 			 if (Remove_all(writes_Ions[j])==Remove_all(outsHelp[i-1]+"i"))
		 			 	 optional = true;
		 		  }

		 		  // write m_F only if there is a equilirium used!
		 		 bool m_F_needed = false;
		 		 for (size_t j=0; j<All_Eqs.size(); j++)
		 		 {
		 			 // if Nernst Pot used we need to add /m_F
		 			 if (outs[i].find(Remove_all(All_Eqs[j]))!=outs[i].npos)
		 				 m_F_needed = true;
		 		 }
		 		 //std::cout << "whats in outs: " << outs[i] << " ahh thats in it" << std::endl;

		 		 //std::cout << "optional: " << optional << " force " << force << std::endl;
		 		 if (m_F_needed==true && (optional==true || force == true))
		 			 mycppfile << "outCurrentValues.push_back(" + outs[i] + "/m_F ); \n";

		 		 if (m_F_needed==false && (optional==true || force == true))
				 {
		 			 mycppfile << "outCurrentValues.push_back(" + outs[i] + " ); \n";
				 }
		 	  }
	  }

	  std::cout << "Ionic_current Function is written" << std::endl;

	  mycppfile << " } \n \n \n";



	  std::cout << "Writing template class instantation and global Vars" << std::endl;

	  mycppfile << "//////////////////////////////////////////////////////////////////////////////// \n";
	  mycppfile << "//	explicit template instantiations \n";
	  mycppfile << "//////////////////////////////////////////////////////////////////////////////// \n";

	  mycppfile << "#ifdef UG_DIM_1 \n";
	  mycppfile << "template class "+filename+"<Domain1d>; \n";
	  mycppfile << "#endif \n \n \n";

	  mycppfile << "#ifdef UG_DIM_2 \n";
	  mycppfile << "template class "+filename+"<Domain2d>; \n";
	  mycppfile << "#endif \n \n \n";


	  mycppfile << "#ifdef UG_DIM_3 \n";
	  mycppfile << "template class "+filename+"<Domain3d>; \n";
	  mycppfile << "#endif \n \n \n";


	  mycppfile << "}  \n";



	  mycppfile << "  \n";
	  mycppfile << "  \n";
	  mycppfile.close();






	  if (Global_vars.size() > 0)
	  {
		  // check for twice entrys
		  for (size_t i = 0; i<Global_vars.size(); i++)
		  {
			  int count = 0;
			  for (size_t j = 0; j<Global_vars.size(); j++)
			  {
				  if (Remove_all(Global_vars[i])==Remove_all(Global_vars[j]))
				  {
					  count +=1;
					  if (count >= 2)
						  Global_vars[j] = "";
				  }
			  }
		  }




		  for (size_t i = 0; i<Global_vars.size(); i++)
		  {
			  //HFile_added_Vars.push_back(Global_vars[i]);
			  if (Global_vars[i]!="")
					  myhfile << "number " + Global_vars[i] + "; \n";
		  }

	  }


	  if (Params_Unit.size() >0)
	  {
		  for (size_t i = 0; i<Params_Unit.size(); i++)
		  {
			  bool is_written = false;
			  //checking that no global which is also Parameter is written twice in h-file!
			  for (size_t j=0; j<Global_vars.size(); j++)
			  {
				  //std::cout << "write vergleich: " << Remove_all(Params_Unit[i].substr(0, Params_Unit[i].find("="))) << " vs "<<Remove_all(Global_vars[j]) << std::endl;
				  if (Remove_all(Params_Unit[i].substr(0, Params_Unit[i].find("=")))==Remove_all(Global_vars[j]))
					  is_written = true;
			  }
			  if (is_written == false)
			  {
				  myhfile << "number " + Params_Unit[i].substr(0, Params_Unit[i].find("=")) + "; \n";
			  }
		  }
	  }

	  myhfile << "number m_R, m_T, m_F; \n";


	  myhfile << "}; \n \n";

	  myhfile << "} // namespace ug \n \n \n";

	  myhfile << "#endif // " + filename + "_H_\n";

	  myhfile.close();
}






////////////////////////////////////////////////////////////////////////////////////////
// Part which is needed for including in UG4 directly
////////////////////////////////////////////////////////////////////////////////////////


std::vector<string> Converter::WriteChannelFile(string Ch_Name, string filename)
{

	std::vector<string> needed_files;

	// writting in Vector channel infos
	bool write_channel = true;
	std::vector<string> ChannelFile = Openfile(filename);

	// checks if converted channel already existing
	for (size_t i=0; i<ChannelFile.size(); i++)
	{
		size_t include = ChannelFile[i].find("#include");
		if (include!=ChannelFile[i].npos)
		{
			needed_files.push_back(ChannelFile[i].substr(0, ChannelFile[i].npos));
			if (ChannelFile[i].find(Ch_Name, include)!=ChannelFile[i].npos)
			{
				write_channel = false;
			}

		}
	}


	const char* filenamechar = filename.c_str();
	ofstream mychannelfile;
	// if not existing write channel in file
	if (write_channel == true)
	{
		mychannelfile.open (filenamechar, std::ios::app);
		mychannelfile << "#include \"" + Ch_Name + ".h\" \n \n \n";
		mychannelfile << "{ \n";
		mychannelfile << "\t typedef " + Ch_Name + "<TDomain> T; \n";
		mychannelfile << "\t typedef IChannel<TDomain> TBase; \n";
		mychannelfile << "\t string name = string(\""+ Ch_Name + "\").append(suffix); \n";
		mychannelfile << "\t reg.add_class_<T, TBase >(name, grp) \n";
		mychannelfile << "\t \t .template add_constructor<void (*)(const char*, const char*)>(\"Function(s)#Subset(s)\") \n";
		mychannelfile << "\t \t .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>(\"Function(s)#Subset(s)\") \n";
		mychannelfile << "\t \t .set_construct_as_smart_pointer(true); \n";
		mychannelfile << "\t reg.add_class_to_group(name, \"" + Ch_Name + "\", tag); \n";
		mychannelfile << "} \n \n \n";
	}
	mychannelfile.close();
	return needed_files;
}

 void Converter::WriteInclude_List(std::vector<string> Includes, string filename, string sources)
 {
	 std::vector<string> IncludeFile = Openfile(filename);

	 const char* filenamechar = filename.c_str();
	 ofstream myincludefile;
	 myincludefile.open(filenamechar, std::ios::app);

	 const char* sourceschar = sources.c_str();
	 ofstream mysourcefile;
	 mysourcefile.open(sourceschar, std::ios::app);

	 for (size_t i=0; i<Includes.size(); i++)
	 {
		 bool write = true;
		 for (size_t j=0; j<IncludeFile.size(); j++)
		 {
			 if (IncludeFile[j]==Includes[i])
				 write = false;
		 }
		 if (write==true)
		 {
			 // writting include file
			 myincludefile << Includes[i] + "\n";
			 // writting channel sources
			 size_t name_beg, name_end;
			 name_beg = Includes[i].find("#include \"") + 10;
			 name_end = Includes[i].find(".") + 1;
			 mysourcefile << Includes[i].substr(name_beg, name_end-name_beg) + ".cpp \n";
		 }
	 }
	 myincludefile.close();
	 mysourcefile.close();
 }

void Converter::WriteInPlugin(string SourceFile, string DestFile)
{
	 std::vector<string> SourceFileV = Openfile(SourceFile);

	 std::vector<string> DestFileV = Openfile(DestFile);

	 const char* DestFilechar = DestFile.c_str();
	 ofstream myDestFile;
	 myDestFile.open(DestFilechar, std::ios::trunc);

	 size_t Zeile;

	 // deletes any sources already written
	 for (size_t i=0; i<DestFileV.size(); i++)
	 {
		 for (size_t j=0; j<SourceFileV.size(); j++)
		 {
			 if (DestFileV[i].find(SourceFileV[j])!=DestFileV[i].npos)
				 SourceFileV[j] = "";
		 }

	 }


	 //finds line for beginning add sources
	 for (size_t i = 0; i<DestFileV.size(); i++)
	 {
		 if (DestFileV[i].find("set(SOURCES")!=DestFileV[i].npos)
			 Zeile = i;
	 }

	 // adding sourcefiles at right line
	 for (size_t i = 0; i<DestFileV.size(); i++)
	 {
		 if (i != Zeile)
		 {
			 myDestFile << DestFileV[i] +"\n";
		 } else
		 {
			 myDestFile << DestFileV[i] +"\n";
			 for (size_t j=0; j<SourceFileV.size(); j++)
			 {
				 if (SourceFileV[j]!="")
				 {
					 myDestFile <<  "Convert/Debug/" + SourceFileV[j] +"\n";
				 }
			 }
	 	 }
	 }
}
