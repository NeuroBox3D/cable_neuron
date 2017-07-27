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

//TODO's
// closer loook to unit conv
// perhaps java interface for converter

Converter::Converter() {
	// TODO Auto-generated constructor stub

}

Converter::~Converter() {
	// TODO Auto-generated destructor stub
}


bool Converter::Check_if(string test)
{
	test = Remove_all(test);
	if (test[0]=='i' && test[1]=='f' && test[2]=='(')
		return true;
	else
		return false;
}



// Some channels need flux informations about ions this is here setted
std::vector<string> Converter::Read_i_value(std::vector<string> Neuron)
{
	std::vector<string> erg;
	//getting all read values
	std::vector<string> All_reads = Find_all_read(Neuron);
	All_reads = Remove_all(All_reads);
	bool written = false;
	for (size_t i=0; i<All_reads.size(); i++)
	{
		if (All_reads[i]=="ica")
		{
			erg.push_back("number ica = m_pCE->flux_ca(); \n");
			written = true;
		}
		if (All_reads[i]=="ik")
		{
			erg.push_back("number ik = m_pCE->flux_k(); \n");
			written = true;
		}
		if (All_reads[i]=="ina")
		{
			erg.push_back("number ina = m_pCE->flux_na(); \n");
			written = true;
		}
		if (All_reads[i]=="i")
		{
			erg.push_back("number i = m_pCE->flux_v(); \n");
			written = true;
		}
	}
	if (written == false)
		erg.push_back("");

	return erg;
}



std::vector<string> Converter::if_handling(size_t begin, std::vector<string> Unit)
{
	std::vector<string> erg;
	size_t end =0;
	// checking if there really is any if
	if (Unit[begin].find("if")!=Unit[begin].npos)
	{
		size_t i = begin;
		// find brackets of if statement
		size_t clip_on, clip_off;
		clip_on = Unit[i].find_first_of("(");
		clip_off = Unit[i].find_last_of(")");
		string statement = Unit[i].substr(clip_on, clip_off-clip_on+1);
		// first line
		erg.push_back("if " + statement + "\n");
		// now we need end of if
		size_t of = 0;
		size_t on = 0;
		size_t beginnew;
		for (size_t j = begin; j<Unit.size(); j++)
		{
			if (Unit[j].find("{")!=Unit[j].npos)
			{
				if (on==0)
					beginnew = j;
				on += 1;
			}
			if (Unit[j].find("}")!=Unit[j].npos)
				of += 1;
			if ((on != 0 && on == of) || Unit[j].find("else")!=Unit[j].npos)
			{
				end = j;
				if (Unit[j].find("}")!=Unit[j].npos && Unit[j].find("else")!=Unit[j].npos)
					end -=1;

				break;
			}

		}
		erg.push_back("{ \n");
		bool done = false;
		for (size_t j = begin; j< end+1; j++)
		{
			// only do once!!
			if (done==false)
			{
				if (begin == beginnew)
					Unit[j] = Unit[j].substr(Unit[j].find("{"));
				if (begin == end)
					Unit[j] = Unit[j].substr(Unit[j].find("{")+1, Unit[j].find("}")-Unit[j].find("{")-1);
				done = true;
			}

			if (Remove_all(Unit[j])!="}" && Remove_all(Unit[j])!="{")
				erg.push_back(Unit[j] + "; \n");
		}
		erg.push_back("} \n");
	}

	std::cout << "end: " << end << std::endl;
	size_t start_else = 0;
	// checking for any else following
	if (Unit[end].find("else")!=Unit[end].npos)
		start_else = end;

	if (Unit[end+1].find("else")!=Unit[end+1].npos)
		start_else = end + 1;

	// if really any else exists
	if (start_else!=0)
	{
		size_t beg;
		size_t of = 0;
		size_t on = 0;
		size_t beginnew;
		erg.push_back("else \n");
		erg.push_back("{ \n");
		//search for end and beginning
		for (size_t j = start_else; j<Unit.size(); j++)
		{

			if (Unit[j].find("else")!=Unit[j].npos)
				beg = Unit[j].find("else");
			else
				beg = 0;

			std::cout << "beg: " << beg << std::endl;
			if (Unit[j].find("{", beg)!=Unit[j].npos)
			{
				if (on==0)
					beginnew = j;

				on += 1;
			}
			if (Unit[j].find("}", beg)!=Unit[j].npos)
				of += 1;

			if (on != 0 && on == of)
			{
				end = j;
				break;
			}

		}

		bool done = false;
		for (size_t j = start_else; j< end; j++)
		{
			if (done==false)
			{
				if (start_else == beginnew)
					Unit[j] = Unit[j].substr(Unit[j].find("{"));
				if (start_else == end)
					Unit[j] = Unit[j].substr(Unit[j].find("{")+1, Unit[j].find("}")-Unit[j].find("{")-1);
				done = true;
			}
			if (Remove_all(Unit[j])!="}" && Remove_all(Unit[j])!="{")
				erg.push_back(Unit[j] + "; \n");
		}


	//erg.push_back("} \n");
	}

	end += 1;
	std::stringstream ss;
	ss << end;
	erg.push_back(ss.str());


	return erg;
}


std::vector<string> Converter::buildgetter(string var, string class_name)
{
	std::vector<string> getter;


	// first all time function head for h file
	getter.push_back("double get" + var + "();");

	// second all time function head for cpp file
	getter.push_back("double "+ class_name + "<TDomain>::get" + var + "()");

	// following complete function line for line
	getter.push_back("{");
	getter.push_back("return "+ var +";");
	getter.push_back("}");
	std::cout << "before return getter;" << std::endl;

	return getter;
}

std::vector<string> Converter::buildsetter(string var, string class_name)
{
	std::vector<string> setter;


	// first all time function head for h file
	setter.push_back("void set" + var + "(double val);");

	// second all time function head for cpp file
	setter.push_back("void "+ class_name + "<TDomain>::set" + var + "(double val)");

	// following complete function line for line
	setter.push_back("{");
	setter.push_back(var + " = val" + ";");
	setter.push_back("}");

	return setter;
}



std::vector<string> Converter::GetFuncTable(std::vector<string> Zeilen)
{
	std::vector<string> erg;

	for (size_t i=0; i<Zeilen.size(); i++)
	{
		if (Zeilen[i].find("FUNCTION_TABLE")!=Zeilen[i].npos)
		{
			erg.push_back(Zeilen[i]);
		}
	}
	return erg;
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
	std::cout << erg << std::endl;

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
	/*for (size_t i=0; i<result.size() ;i++)
	{
		std::cout << "result: " <<result[i] << std::endl;
	}*/

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
	// existing units simens: pS, mS, S
	// existing units concentrats: mM,
	// existing units temp: degC
	// (um) = (micron) ????
	// mho/cm2 mA/cm2 S/cm2
	// 2 standing for

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
		// Rekursion koennte hier mehrfachaufrufe von prozeduren moeglich machen
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
	out.push_back("const number helpV = 1e3*(m_pCE->R*m_pCE->temperature())/m_pCE->F;");

	size_t Ion, IonRead, IonRend;
	string IonS;
	string PartR, PartW;

	size_t R_beg, R_end;


	std::vector<string> ListIonRead;
	std::vector<string> ListIon;
	std::vector<bool> version;

	std::vector<string> NEURON = GetBlockFunction(Pairs, Zeilen, "NEURON");

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
				if (IonRead!=NEURON[i].npos)
				{
					IonRend = NEURON[i].find(" ", IonRead+5);
					//std::cout << IonRend << ", " << Ion << std::endl;
					IonS = NEURON[i].substr(Ion+7, IonRead-1-(Ion+7));
					//std::cout << "Subion: "<< IonS << std::endl;
					ListIonRead.push_back(NEURON[i].substr(IonRead+5 , IonRend-(IonRead+5)));
					ListIon.push_back(IonS);
				}
				else
				{
					// if there is no read but a write we have also to add right way!
					IonS = NEURON[i].find("WRITE");
				}
				// Only if there is a WRITE in
				if (NEURON[i].find("WRITE")!=NEURON[i].npos)
					ListIon.push_back(IonS);

			}
		}
	}

		for (size_t j=0; j<ListIonRead.size(); j++)
		{
			std::cout << "number of lists: " << j << std::endl;
			if (version[j]==false)
			{
			std::cout << "read ion: " + ListIonRead[j] + "List ion: " + ListIon[j] << std::endl;
			out.push_back("number " + ListIonRead[j] + " = helpV*(log(m_pCE->" + ListIon[j] + "_out()/" + ListIon[j] + "));");
			}
		}


		for (size_t j=0; j<ListIonRead.size(); j++)
		{
			if (version[j]==true)
			{
			std::cout << "read ion: " + ListIonRead[j] + "List ion: " + ListIon[j] << std::endl;
			out.push_back("number " + ListIonRead[j] + " = helpV*(log(m_pCE->" + ListIon[j] + "_out()/" + ListIon[j] + "));");
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
					if (PROCEDURE[j].find("{")<PROCEDURE[j].find("LOCAL"))
						PROCEDURE[j].replace(0, PROCEDURE[j].find("{")+1, "");

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


// liest einzelne NModl-File ein
std::vector<string> Converter::Openfile(string filename)
{

	std::vector<string> zeilen;
	//try {

	const char* filenamechar = filename.c_str();
	ifstream ModFile(filenamechar, std::ios::out);


	while (!ModFile.eof())          // Solange noch Daten vorliegen
	{
		////std::cout << "while" << std::endl;
		string zeile;
	    getline(ModFile, zeile);        // Lese eine Zeile
	    if (zeile.find(":")!=zeile.npos)
	    {
	    	if (Remove_all(zeile).find(":")<=3)
	    		zeile.replace(zeile.find(":"), 1,"//--//");
	    	else
	    		zeile.replace(zeile.find(":"), 1,";//");
	    }// replace : with // as comment
	    //replace FARADAY with F
	    if (zeile.find("FARADAY")!=zeile.npos)
	    	zeile.replace(zeile.find("FARADAY"), 7, "F");

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

	bool failure;

	for (size_t i=0; i<Zeilen.size(); i++)
	{
		failure = false;

		start = Zeilen[i].find("{");
		ende = Zeilen[i].find("}");

		if (start!=Zeilen[i].npos && ende!=Zeilen[i].npos) // true, if both brackets have been found on line i
		{
			if (Zeilen[i].find("STATE")==Zeilen[i].npos)	// if STATE has not been found, ignore line i
				failure=true;
		}

		if (failure == false)
		{
			if (start!=Zeilen[i].npos)	// '{' found on line i
			{
				if (works==false)
				{
					blockH = i;		// beginning of a block
					works = true;		// ignore '{'s until a '}' has been found
				}
			}
			if (ende!=Zeilen[i].npos)
			{
				Blocks.push_back(std::make_pair(blockH, i+1));
				works = false;
			}
		}

	}
	return Blocks;

}


// Writing needed Block out of Nmodl-file
std::vector<string> Converter::GetBlockFunction(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name)
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
			//ueberprueft wenn weitere klammern da sind, dass nur die letzte verwendet wird;

			while (Pairs[i].first == Pairs[i+1].first)
			{
				std::cout << "vergleich " << Pairs[i].first << " - " << Pairs[i+1].first << std::endl;
				////std::cout << "in whhile by i: " << i << std::endl;
				i = i+1;
			}
			//ueberspringt doppelte Pairs
			//std:: cout << counter << " , " << Pairs[i+counter].second <<std::endl;

			size_t breaks;
			std::cout << "pair first: " << Pairs[i].first << " pair second: " << Pairs[i].second << std::endl;

			//if there are problems read to end of file
			if (Pairs[i].second==0)
				Pairs[i].second=Pairs[i-1].second;

			// schreibt zeilen raus
			for (int j=Pairs[i].first; j<Pairs[i].second; j++)
			{
				std::cout << name << " used pairs: " << Pairs[i].first << " - " << Pairs[i].second << std::endl;
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
						std::cout << "^ gefunden" << std::endl;
						// writes into left every part of one number
						left = "";

						int clipclose, clipopen;
						bool firstclip = false;
						// if left of ^ is a ) we have to get much more of left side until clips are closed
						if (zeile_without[breaks-1]==')')
						{
							firstclip=true;
							clipclose = 1;
							clipopen = 0;
						}



						//std::cout << "breaks: " << breaks << std::endl;
						for (size_t k=breaks-1; k>1; k--)
						{
							//sleep(10);
							//std::cout<< "k schleife with k: " << k << std::endl;
							std::cout << zeile_without[k] << std::endl;

							//normal case
							if (firstclip==false)
							{
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

							//special case we need closing brackets
							if (firstclip==true)
							{
								if (zeile_without[k]=='(')
									clipopen += 1;
								if (zeile_without[k]==')')
									clipclose += 1;
								//break condition
								if (clipopen==clipclose)
								{
									k = 1;
								}
								else
								{
									left = left + zeile_without[k];
								}


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
							// if the problem has to be solved without clips
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

								if ((is_single_number(zeile_without[k])==true) || (is_single_letter(zeile_without[k])==true)
										|| (zeile_without[k]=='(') || (zeile_without[k]==')'))
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

						//left side of = has to be same
						//Problem only right site has to replace not left!!!

						new_zeile = zeile_without.replace(zeile_without.find(right, breaks), right.size(), "");
						//std::cout << "new_zeile1: " << new_zeile << std::endl;
						new_zeile = new_zeile.replace(new_zeile.find(left, (breaks-(left.size()+1))), left.size(), "");
						//std::cout << "new_zeile2: " << new_zeile << std::endl;
						new_zeile = new_zeile.replace(new_zeile.find("^"), 1, " pow(" + left + " , " + right + ")");
						//std::cout << "new_zeile3: " << new_zeile << std::endl;
						////std::cout << "while schleife geht net mit Zeilen[j] = " << Zeilen[j] << std::endl;

						Zeilen[j] = new_zeile;
						breaks = Zeilen[j].find("^");
						////std::cout << "breaks " << breaks << std::endl;
						//end of converting action

					}



		////////////// deletes line with "at_time()"
					if (Zeilen[j].find("at_time")!=Zeilen[j].npos)
					{
						Zeilen[j] = "";
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
			//ueberprueft wenn weitere klammern da sind, dass nur die letzte verwendet wird;

			while (Pairs[i].first == Pairs[i+1].first)
			{
				////std::cout << "in whhile by i: " << i << std::endl;
				i = i+1;
			}
			//ueberspringt doppelte Pairs
			//std:: cout << counter << " , " << Pairs[i+counter].second <<std::endl;

			size_t breaks;

			//if beginning is end write one out
			if (Pairs[i].first == Pairs[i].second)
				Pairs[i].second += 1;

			// schreibt zeilen des blocks raus
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
						std::cout << "^ gefunden" << std::endl;
						// writes into left every part of one number
						left = "";

						int clipclose, clipopen;
						bool firstclip = false;
						// if left of ^ is a ) we have to get much more of left side until clips are closed
						if (zeile_without[breaks-1]==')')
						{
							firstclip=true;
							clipclose = 1;
							clipopen = 0;
						}



						//std::cout << "breaks: " << breaks << std::endl;
						for (size_t k=breaks-1; k>1; k--)
						{
							//sleep(10);
							//std::cout<< "k schleife with k: " << k << std::endl;
							std::cout << zeile_without[k] << std::endl;

							//normal case
							if (firstclip==false)
							{
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

							//special case we need closing brackets
							if (firstclip==true)
							{
								if (zeile_without[k]=='(')
									clipopen += 1;
								if (zeile_without[k]==')')
									clipclose += 1;
								//break condition
								if (clipopen==clipclose)
								{
									k = 1;
								}
								else
								{
									left = left + zeile_without[k];
								}


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
							// if the problem has to be solved without clips
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

								if ((is_single_number(zeile_without[k])==true) || (is_single_letter(zeile_without[k])==true)
										|| (zeile_without[k]=='(') || (zeile_without[k]==')'))
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

						//left side of = has to be same
						//Problem only right site has to replace not left!!!

						new_zeile = zeile_without.replace(zeile_without.find(right, breaks), right.size(), "");
						//std::cout << "new_zeile1: " << new_zeile << std::endl;
						new_zeile = new_zeile.replace(new_zeile.find(left, (breaks-(left.size()+1))), left.size(), "");
						//std::cout << "new_zeile2: " << new_zeile << std::endl;
						new_zeile = new_zeile.replace(new_zeile.find("^"), 1, " pow(" + left + " , " + right + ")");
						//std::cout << "new_zeile3: " << new_zeile << std::endl;
						////std::cout << "while schleife geht net mit Zeilen[j] = " << Zeilen[j] << std::endl;

						Zeilen[j] = new_zeile;
						breaks = Zeilen[j].find("^");
						////std::cout << "breaks " << breaks << std::endl;
						//end of converting action

					}

		// Converting all if-loops into one line
					if (Remove_all(Zeilen[j]).find("if(")!=Zeilen[j].npos)
					{
						//std::cout << "Start if reworker" << std::endl;
						bool withelse = false;
						size_t NumberOfLines = 1;
						size_t else_line = 0;

						for (size_t k = j; k < (size_t) Pairs[i].second; k++)
						{
							//std::cout << "searching else" << std::endl;
							//std::cout << "from " << k << " to " << Pairs[i].second << std::endl;
							//checking for else
							if (Zeilen[k].find("else")!=Zeilen[k].npos)
							{
								//std::cout << "finding else" << std::endl;
								withelse=true;
								NumberOfLines += 1;
								else_line = k;
								// break condition
								break;
							}
						}
						if (withelse==true)
						{
							Zeilen[else_line].replace(Zeilen[else_line].find("}"),1 ,";}");
							// checking for closing }
							size_t close = else_line + 1;
							NumberOfLines += 1;
							while (Zeilen[close].find('}')==Zeilen[close].npos)
							{
								//std::cout << Zeilen[close] << std::endl;
								//std::cout << "in while" << std::endl;
								close +=1;
								NumberOfLines +=1;
								//std::cout << "Zeilefind: " << Zeilen[close].find('}') << std::endl;
							}
						}

						size_t closewithout;

						if (withelse==false)
						{
							NumberOfLines = 1;
							closewithout = j;
							while (Zeilen[closewithout].find('}')==Zeilen[closewithout].npos)
							{
								closewithout += 1;
								NumberOfLines +=1;
							}
						}

						//std::cout << "else?? " << withelse << std::cout;
						//std::cout << "before writing linesnumbers " << NumberOfLines << std::endl;
						for (size_t k=j; k<NumberOfLines+j; k++)
						{
							if ((size_t)j != k)
							{
								Zeilen[j] = Zeilen[j] + Zeilen[k];
								Zeilen[k] = "\n";
							}
							//std::cout << "second if" << std::endl;
							if (k==(NumberOfLines+j-1))
								if (Zeilen[j].find_last_of("}") != Zeilen[j].npos)
								{
									Zeilen[j].replace(Zeilen[j].find_last_of("}"),1 ,";}");
									Zeilen[j] = Zeilen[j] + ";";
								}
						}

						//std::cout << "end-if reworker" << std::endl;
					}




		////////////// deletes line with "at_time()"
					if (Zeilen[j].find("at_time")!=Zeilen[j].npos)
					{
						Zeilen[j] = "";
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
	  std::vector<string> Zeilen2 = Zeilen;


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


	  mycppfile << "namespace ug { \n";
	  mycppfile << "namespace cable_neuron { \n \n \n";

	  std::vector<string> mod_funcs_names;
	  // all function and function_table have to be read in first so later need more attention more

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Adding of Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////


	  std::vector<string> FUNCTION = GetBlockFunction(Pairs, Zeilen2, "FUNCTION");

	  string without[] = {"if", "else"};
	  size_t funcS, funcE;
	  bool second_time;

	  std::vector<string> if_vars1;

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

				  size_t brack_open, brack_close;
				  brack_open = 0;
				  brack_close = 0;

				  // testing if brackets in func are set correct if not correcting
				  for (size_t k=0; k<func.size(); k++)
				  {
					  if (func[k]=='(')
						  brack_open += 1;

					  if (func[k]==')')
						  brack_close += 1;

				  }

				  // if there is a closed bracket to much
				  if (brack_close > brack_open)
				  {
					  size_t close_pos = func.find_last_of(")");
					  func.replace(close_pos, 1, "");
				  }

				  func_hfile.push_back(func);


				  // writing needed templates in function!
				  string left = func.substr(0, func.find("::"));
				  string right = func.substr(func.find("::")+2, func.npos-(func.find("::")+2));
				  func = left + "<TDomain>::" + right;
				  mycppfile << "template<typename TDomain> \n";

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
						  funcE = FUNCTION[i].find("=", funcS);
						  FUNCTION[i].replace(funcS, funcE-funcS + 1, "return ");
					  }

					  if (FUNCTION[i].find("FUNCTION")==FUNCTION[i].npos)
					  {
						  if (FUNCTION[i].find("LOCAL")!=FUNCTION[i].npos)
						  {
							  FUNCTION[i].replace(FUNCTION[i].find("LOCAL"), 6, "double ");
						  }
						  if (second_time==false)
						  {

							  if (FUNCTION[i].find("if")!=FUNCTION[i].npos)
							  {
								  std::cout << "first if!!" << std::endl;
								  if (i+1<=FUNCTION.size())
								  {
									  if (FUNCTION[i+1].find("=")!=FUNCTION[i+1].npos)
									  {
										  std::cout << "push_back!! " << std::endl;
										  if_vars1.push_back(Remove_all(FUNCTION[i+1].substr(0, FUNCTION[i+1].find("="))));
									  }
								  }
							  }


							  if ((FUNCTION[i].find(without[0])==FUNCTION[i].npos) && (FUNCTION[i].find(without[1])==FUNCTION[i].npos)
								 && (string(FUNCTION[i] + "; \n").find("};")==string(FUNCTION[i] + "; \n").npos) && (FUNCTION[i]!="")
								 && Remove_all(FUNCTION[i])!="}")
							  {

								  	  bool var_in_ii = false;
									  // normal vars handling
									  if (FUNCTION[i].find(" = ")!=FUNCTION[i].npos)
									  {
										  if (if_vars1.size()>0)
										  {
											  var_in_ii=false;
											  for (size_t ii=0; ii<if_vars1.size(); ii++)
											  {
												  if (FUNCTION[i].find(if_vars1[ii])!=FUNCTION[i].npos)
													  var_in_ii=true;
											  }
										  }
										  if (var_in_ii==true)
										  {
											  std::cout << FUNCTION[i] << std::endl;
											  mycppfile << FUNCTION[i] + "; \n";
										  }
										  else
										  {
											  std::cout << "double " + FUNCTION[i] << std::endl;
											  mycppfile << "double " + FUNCTION[i] + "; \n";
										  }
									  }
									  else
									  {
										  std::cout << FUNCTION[i] << std::endl;
									////std::cout << FUNCTION[i];
										  mycppfile << FUNCTION[i] + "; \n";
									  }
								  //}
							  }
							  else
							  {
								  if ((if_vars1.size()>0) and (FUNCTION[i].find(without[0])!=FUNCTION[i].npos))
								  {
									  //bool adden = false;
									  for (size_t ii=0; ii<if_vars1.size(); ii++)
									  {
										  std::cout << "laufen ii: " << ii << std::endl;
										  if (i+1<=FUNCTION.size())
										  {
											  std::cout << FUNCTION[i] << " - " << FUNCTION[i+1] << std::endl;
											  std::cout << if_vars1[ii] << " - " << Remove_all(FUNCTION[i+1].substr(0, FUNCTION[i+1].find("="))) << std::endl;
											  if (if_vars1[ii]==Remove_all(FUNCTION[i+1].substr(0, FUNCTION[i+1].find("="))))
											  {
												  std::cout << "mano " << std::endl;
												  mycppfile << "double " << if_vars1[ii] << "; \n";
											  }
										  }
									  }
									  std::cout << FUNCTION[i] << std::endl;
									  ////std::cout << FUNCTION[i];
									  mycppfile << FUNCTION[i] + "\n";
								  }
								  else
								  {
									  std::cout << FUNCTION[i] << std::endl;
									  ////std::cout << FUNCTION[i];
									  mycppfile << FUNCTION[i] + "\n";
								  }
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

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Adding of Function tables
/////////////////////////////////////////////////////////////////////////////////////////////////////
	  std::vector<string> FUNCTION_TABLE = GetFuncTable(Zeilen);


	  string Table_func_head_complet;

	  std::vector<string> all_func_vars;
	  if (FUNCTION_TABLE.size() > 0)
	  {
		  //std::cout << "Function table is happening" << std::endl;
		  for (size_t i=0; i< FUNCTION_TABLE.size(); i++)
		  {
			  // getting function head
			  if (FUNCTION_TABLE[i].find("FUNCTION_TABLE")!=FUNCTION_TABLE[i].npos)
			  {
				  std::cout << "is working" << std::endl;
				  string FT_head =  Remove_all(FUNCTION_TABLE[i].substr(FUNCTION_TABLE[i].find("FUNCTION_TABLE")+14 ,FUNCTION_TABLE[i].find("(") - (FUNCTION_TABLE[i].find("FUNCTION_TABLE")+14) ));

				  Table_func_head_complet += "double " + FT_head;



				  size_t varsbeg = FUNCTION_TABLE[i].find(FT_head) + FT_head.size();
				  int clips = 0;
				  string vars = "";
				  // counts clips
				  //std::cout << "befor clips" << std::endl;
				  for (size_t j=varsbeg; j<FUNCTION_TABLE[i].size(); j++)
				  {
					  if (FUNCTION_TABLE[i][j] == ' ')
						  j +=1;
					  if (FUNCTION_TABLE[i][j]=='(')
						  clips += 1;
					  if (FUNCTION_TABLE[i][j]==')')
						  clips -= 1;
					  std::cout << "clips: "<< clips << std::endl;
					  if (clips <= 0)
					  {
						  vars = Remove_all(FUNCTION_TABLE[i].substr(varsbeg, j-varsbeg));
						  break;
					  }
				  }
				  //std::cout << "vars: " << vars << std::endl;
				  // time vor parsing vars in single vars

				  if (vars[0]=='(')
					  vars[0] = ' ';
				  vars = Remove_all(vars);
				  // now delete all Units
				  size_t clip_on = vars.find("(");
				  size_t clip_off = vars.find(")");

				  // if there is an opend and close clip
				  while ((clip_on!=vars.npos) && (clip_off!=vars.npos))
				  {
					  //std::cout << "clip_on: " << clip_on << " clip_off: " << clip_off << std::endl;
					  vars.replace(clip_on, clip_off+1-clip_on, "");
					  clip_on = vars.find("(");
					  clip_off = vars.find(")");
					  //std::cout << "vars zwischen: " <<vars << std::endl;
				  }
				  //std::cout<< vars << std::endl;
				  size_t komma = vars.find(",");
				  size_t beg = 0;
				  while (komma!=vars.npos)
				  {
					  all_func_vars.push_back(vars.substr(beg, komma-beg));
					  beg = komma+1;
					  komma = vars.find(",", beg);
				  }
				  all_func_vars.push_back(vars.substr(beg, komma-beg));
			  }

		  }
	  }

	  // adding function with func infos into cppfile and hfile




	  if (all_func_vars.size()>0)
	  {
		  Table_func_head_complet += "(";
		  for (size_t i=0; i<all_func_vars.size(); i++)
		  {
			  if (i+1 == all_func_vars.size())
				  Table_func_head_complet += "double " + all_func_vars[i];
			  else
				  Table_func_head_complet += "double " + all_func_vars[i] + ", ";
		  }


		  Table_func_head_complet += ")";

		  mycppfile << "template<typename TDomain> \n";
		  mycppfile << "double "+filename+"<TDomain>::" + Table_func_head_complet.substr(Table_func_head_complet.find("double ")+7) + "\n";
		  mycppfile << "{  \n";
		  mycppfile << "return "+ all_func_vars[0] +"; \n";
		  mycppfile << "} \n \n";

	  }


	  mycppfile << "// adding function which always inits_attachments \n";
	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void "+filename+"<TDomain>::ce_obj_available()  \n";
	  mycppfile << "{  \n";
	  mycppfile << "	init_attachments();  \n";
	  mycppfile << " 	F = m_pCE->F; \n R = m_pCE->R; \n K = m_pCE->temperature(); \n celsius = m_pCE->temperature_celsius(); \n";
	  mycppfile << "}  \n \n \n \n";
	  mycppfile.close();

	  ofstream myhfile;
	  myhfile.open (filenameh);
	  myhfile << "#ifndef " + filename + "_H_\n";
	  myhfile << "#define " + filename + "_H_\n";

	  myhfile << "#include \"../../cable_membrane_transport_interface.h\" \n";

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

	  myhfile << "#include \"../../../cable_disc/cable_equation.h\" \n \n";

	  myhfile << "#include <vector> \n";
	  myhfile << "#include <stdio.h> \n";
	  myhfile << "#include \"bindings/lua/lua_user_data.h\" \n";

	  myhfile << "namespace ug {\n";
	  myhfile << "namespace cable_neuron {\n\n\n";
	  myhfile << "// forward declaration \n";
	  myhfile << "template <typename TDomain> \n";
	  myhfile << "class CableEquation; \n \n";

	  myhfile << "template <typename TDomain> \n";
	  myhfile << "class " + filename + "\n";
	  myhfile << "    : public ICableMembraneTransport<TDomain> \n";
	  myhfile << "{ \n";
	  myhfile << "    public: \n ";
	  myhfile << "using ICableMembraneTransport <TDomain>::m_pCE; \n \n";


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
// Also setting all clear Parameter to 0 and preparing setters and getters for every Parameter
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
					  if (PARAMETER[i].find(";//")!=PARAMETER[i].npos)
					  {
						  PARAMETER[i] = PARAMETER[i].substr(0, PARAMETER[i].find(";//")-1);
					  }
					  Unit = 1;
				  }

				  var = PARAMETER[i].substr(0, end);
				  if (end!= PARAMETER[i].npos && endUnit!= PARAMETER[i].npos)
				  {
					  Unit = Unit_Conv_All(PARAMETER[i].substr(end, endUnit));
				  }
				  Units << Unit;
				  HFile_added_Vars.push_back(var.substr(0, var.find("=")));
				  // Hard-Coded VM-Params can not be overiden
				  if ((Remove_all(var.substr(0, var.find("=")))!="celsius")
					 && (Remove_all(var.substr(0, var.find("=")))!="F")
					 && (Remove_all(var.substr(0, var.find("=")))!="R")
					 && (Remove_all(var.substr(0, var.find("=")))!="K")
					 )
				  {
					  Params_Unit.push_back(var + "*" + Units.str());
				  	  Parameter.push_back(var.substr(0, var.find("=")));
				  }
			  }
			  // there are sometimes parameters which do not get a value
			  // but they need without the converted file is not working
			  // also standard parameters like celsius farraday are used here
			  // but they get value from CableEquation-class
			  else
			  {
				  string Param;
				  if (PARAMETER[i].find("=")==PARAMETER[i].npos)
				  {
					  // if there is any comment in beginning of line kill line
					  string param_ = Remove_all(PARAMETER[i]);
					  if (param_.find(":")==0)
						  PARAMETER[i] = "";
					  // if there is any comment in line
					  if (PARAMETER[i].find(":")!=PARAMETER[i].npos)
						  Param = Remove_all(PARAMETER[i].substr(0, PARAMETER[i].find(":")-1));
					  else
						  Param = Remove_all(PARAMETER[i]);

					  // Remove Units
					  if (Param.find("(")!=Param.npos)
						  Param = Param.substr(0, Param.find("("));

					  Param = Remove_all(Param);
					  // hard coded params like celsius should not be overriden
					  if (Param!="celsius" && Param!="dt" && Param!="FARADAY" && Param!="v" && Param!="" && Param!=" "
							  && Param!="na" && Param!="ca" && (Param.find("PARAMETER")==Param.npos) && Param!="}")
					  {
						  // only var name
						  Parameter.push_back(Param);
						  HFile_added_Vars.push_back(Param);
						  // changing Param to 0 value
						  Param += " = 0";
						  // writing complete Param
						  Params_Unit.push_back(Param);
					  }

				  }
			  }
		  }
	  }
	  //myhfile << "std::vector<number> m_diff; \n";
	  myhfile << "\n";

	  // Pointers needs the same treatments as params without values
	  // so we add here the pointers out of Neuron-Block
	  if (NEURON.size()>0)
	  {
		  for (size_t i=0; i<NEURON.size(); i++)
		  {
			  if (NEURON[i].find("POINTER")!=NEURON[i].npos)
			  {
				  size_t PointEnd = NEURON[i].find(":");
				  string Point = Remove_all(NEURON[i].substr(NEURON[i].find("POINTER")+7, PointEnd));
				  Parameter.push_back(Point);
				  HFile_added_Vars.push_back(Point);
				  Point += " =0";
				  Params_Unit.push_back(Point);
			  }
		  }
	  }





//////////////////////////////////////////////////////////////////////////////////////////////
	  /// Writing all ions out
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
			  //std::cout << "Ionsout: " << IonS + "_out" << std::endl;
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


//////////////////////////////////////////////////////////////////////////////////////////////
	 /// Building State_Vars for beginning
//////////////////////////////////////////////////////////////////////////////////////////////

  	  std::vector<string> STATE = GetBlock(Pairs, Zeilen, "STATE");

  	  std::cout << "before new state" << std::endl;

  	  if (STATE.size()>0)
  	  {
  		  for (size_t i=0; i<STATE.size(); i++)
  		  {
  			  if (STATE[i].find("//")!=STATE[i].npos)
  			  {
  				  STATE[i] = STATE[i].substr(0, STATE[i].find("//")-1);
  			  }
  		  }
  	  }

	  // needed var for Gating params having same value of used Ion
	  std::vector<string> SGating;

  	  string addState;

  	  // var in which all states and Non spec currents are saved
  	  std::vector<string> State_vars;
  	  size_t state, stateend, stateend2, komm_beg;
  	  size_t varSt = 0;
  	  bool clip = false;
  	 // size_t tabs;

  	  std::cout << "before states" << std::endl;
  	  std::cout << "State_size " << STATE.size() << std::endl;
  	  //std::cout << "State0: " << STATE[0] << std::endl;

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
  			  std::cout << "number vars: " << varSt <<std::endl;
  			  if (varSt>0)
  			  {
  				  for (size_t j=0; j<varSt; j++)
  				  {
  					  ////std::cout << "in state " << std::endl;
  					  komm_beg = STATE[i].find(":");
  					  stateend = STATE[i].find(" ", state);
  					  stateend2 = STATE[i].find("}", state);
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
  						  if ((addState!="}") && (addState!=" ") && (addState!="") && (addState!="\n") &&(addState!="\n}") && (addState!="}\n")
  								  && (addState.find("}")==addState.npos) && (addState.find("{")==addState.npos))
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
  							  std::cout << "State: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  " << addState << std::endl;
  						  }
  					  }
  					  state = stateend+1;
  					  //////////////////
  				  }
  			  }
  			  else
  			  {
  				  stateend2=1;
  			  }
  			  varSt = 0;
  			  std::cout << "state really 0" << std::endl;
  			  //adding last state if there is one
  			  if (STATE[i][stateend2-1]!=' ' && stateend2!=STATE[i].npos)
  			  {
  				  std::cout << "state number: " << stateend2 << std::endl;
  				  addState = STATE[i][stateend2-1];
  				  std::cout << "komisch " << addState << std::endl;
	  			  if ((addState!="}") && (addState!="") && (addState!=" ") && (addState!="\n") &&(addState!="\n}") && (addState!="}\n")
	  			  && (addState.find("}")==addState.npos)&& (addState.find("{")==addState.npos) && pos_letter(addState)!=1000)

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
	  				  std::cout<< "state added: " << addState << std::endl;
	  				  State_vars.push_back(addState);
	  			  }
  			  }

  		  }
  	  }




	  // for all Parameters we need getter and setters (also for Pointers)
	  std::vector<vector <string> > Getters;
	  std::vector<vector <string> > Setters;

	  // Preparing for all Parameters getters and setters
	  for (size_t i = 0; i<Parameter.size(); i++)
	  {
		  if (Parameter[i].find("//--//")==Parameter[i].npos)
		  {
			  std::cout << Parameter[i] << std::endl;
			  Getters.push_back(buildgetter(Remove_all(Parameter[i]), filename));
			  Setters.push_back(buildsetter(Remove_all(Parameter[i]), filename));
		  }
	  }



	  	//Setters

	  std::cout << "All Parameters set" << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////
// Writing all interface parts
//////////////////////////////////////////////////////////////////////////////////////////////


	  myhfile << "\n"; myhfile << "\n";

	  myhfile << "/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) \n";
	  myhfile << filename + "(const char* functions, const char* subsets) \n";
	  myhfile << "try : ICableMembraneTransport<TDomain>(functions, subsets), \n";
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
		  {
			  if (State_vars.size()==0)
				  add = add + "{} \n";
			  else
				  add = add + ", \n";
		  }
		  if (add.find("//--//")==add.npos)
			  myhfile << add;
	  }
	  for (size_t i=0; i<State_vars.size(); i++)
	  {
		  if (i+1<State_vars.size())
			  myhfile << "m_log_" << State_vars[i] << "Gate(false), \n";
		  else
			  myhfile << "m_log_" << State_vars[i] << "Gate(false) {} \n";
  	  }
	  myhfile << "UG_CATCH_THROW(\"Error in "+ filename + " initializer list. \"); \n \n \n";

	  myhfile << "/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) \n";
	  myhfile << filename + "(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) \n";
	  myhfile << "try : ICableMembraneTransport<TDomain>(functions, subsets), \n";
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
	  	{
	    if (State_vars.size()==0)
	    	add1 = add1 + "{} \n";
		else
			add1 = add1 + ", \n";
	  	}
	  	if (add1.find("//--//")==add1.npos)
	  		myhfile << add1;
	  }
	  for (size_t i=0; i<State_vars.size(); i++)
	  {
		  if (i+1<State_vars.size())
			  myhfile << "m_log_" << State_vars[i] << "Gate(false), \n";
		  else
			  myhfile << "m_log_" << State_vars[i] << "Gate(false) {} \n";
  	  }
	  myhfile << "UG_CATCH_THROW(\"Error in "+ filename + " initializer list. \"); \n";

	  myhfile << "/// destructor \n \n";
	  myhfile << "virtual ~"+ filename + "() {}; \n";


	  std::cout << "bla" << std::endl;

	  string name;

	  // adding function out of function list
	  for (size_t i=0; i<func_hfile.size() ; i++)
	  {
		  name = func_hfile[i].substr(func_hfile[i].find(filename), func_hfile[i].find("::")-func_hfile[i].find(filename)+2);
		  //std::cout << "name: " << name << std::endl;
		  func_hfile[i].replace(func_hfile[i].find(name), name.size(), "");
		  myhfile << func_hfile[i] + "; \n";
	  }


	  if (FUNCTION_TABLE.size()>0)
	  {
		  // adding Table function
		  myhfile << Table_func_head_complet + ";\n";
	  }


	  /// adding functions of ICableMembraneTransport
	  myhfile << "/// create attachments and accessors \n";
	  myhfile << "void init_attachments(); \n";
	  myhfile << "// inherited from ICableMembraneTransport \n \n";
	  myhfile << "virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); \n";
	  myhfile << "virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); \n";
	  myhfile << "virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); \n";
	  myhfile << "virtual void ce_obj_available(); \n";
	  myhfile << "virtual std::vector<number> state_values(number x, number y, number z); \n";
	  myhfile << "\n \n";



	  std::cout << "before getters!" << std::endl;


	  // Writting all parameter getter and setters
	  for (size_t i = 0; i<Getters.size(); i++)
	  {
		  myhfile << Getters[i][0] + " \n";
	  }

	  std::cout << "getters are builded!" << std::endl;

	  for (size_t i = 0; i<Setters.size(); i++)
	  {
		  myhfile << Setters[i][0] + " \n";
	  }

	  std::cout << "setters are builded!" << std::endl;

	  // writing function for m_log_gates
	  for (size_t i=0; i<State_vars.size(); i++)
	  {
  		  myhfile << "void set_log_" << State_vars[i] << "Gate(bool bLog" << State_vars[i] << "Gate); \n";
  	  }

	  myhfile << "\n \n";
	  myhfile << "protected: \n";

	  //myhfile << "void register_func(); \n \n \n";
	  myhfile << "private: \n \n";
	  // Neuron-lines with use ion
	  myhfile << "virtual void specify_write_function_indices(); \n";



	  std::cout << "getters and setters are builded!" << std::endl;





	  // creating init Attachments
	  mycppfile.open (filenamecpp, std::ios::app);


//////////////////////////////////////////////////////////////////////////////
// Adding getters and setters in cpp file
//////////////////////////////////////////////////////////////////////////////

	  for (size_t i = 0; i<Getters.size(); i++)
	  {
	  // starting by one cause 0 is for h-file
		  mycppfile << "template<typename TDomain> \n";
		  for (size_t j=1; j<Getters[i].size(); j++)
	  	  {
			  std::cout << i << " - " <<j << std::endl;
			  std::cout << "write " << Getters[i][j] << std::endl;
			  mycppfile << Getters[i][j] + " \n";
	  	  }
	  }

	  for (size_t i = 0; i<Setters.size(); i++)
	  {
	  // starting by one cause 0 is for h-file
		  mycppfile << "template<typename TDomain> \n";
		  for (size_t j=1; j<Setters[i].size(); j++)
	  	  {
			  std::cout << i << " - " <<j << std::endl;
			  std::cout << "write " << Setters[i][j] << std::endl;
			  mycppfile << Setters[i][j] + " \n";
	  	  }
	  }



//////////////////////////////////////////////////////////////////////////////
// Start Writting interface-main-functions
//////////////////////////////////////////////////////////////////////////////







	  // creating init Attachments
	  mycppfile << " // creating Method for attachments \n";
	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::init_attachments() \n";
	  mycppfile << "{ \n";

	  mycppfile << "SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid(); \n";

////////////////////////////////////////////////////////////////////////////////////
//// Saving all States in one Var also writting first states into cpp and h file
/////////////////////////////////////////////////////////////////////////////////////
	  /*size_t Ns_Current, Ns_CurrentEnd;
	  string Ns_CurrentS;*/

	  varSt = 0;
	  clip = false;
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
  			  std::cout << "number vars: " << varSt <<std::endl;
  			  for (size_t j=0; j<varSt; j++)
  			  {
  				  ////std::cout << "in state " << std::endl;
  				  komm_beg = STATE[i].find(":");
  				  stateend = STATE[i].find(" ", state);
  				  stateend2 = STATE[i].find("}", state);
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
  		  			  if ((addState!="}") && (addState!=" ") && (addState!="") && (addState!="\n") &&(addState!="\n}") && (addState!="}\n")
  		  			  && (addState.find("}")==addState.npos) && (addState.find("{")==addState.npos))
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
  			  //adding last state if there is one
  			  if (STATE[i][stateend2-1]!=' ' && stateend2!=STATE[i].npos)
  			  {
  				  std::cout << "state number: " << stateend2 << std::endl;
  				  addState = STATE[i][stateend2-1];
  				  std::cout << "komisch " << addState << std::endl;
	  			  if ((addState!="}") && (addState!="") && (addState!=" ") && (addState!="\n") &&(addState!="\n}") && (addState!="}\n")
	  			  && (addState.find("}")==addState.npos)&& (addState.find("{")==addState.npos) && pos_letter(addState)!=1000)

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
	  				  std::cout<< "state added: " << addState << std::endl;
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

  		  }
  	  }
	  mycppfile << "} \n \n \n \n";


	  /*
	   * 		  				  myhfile << "ADouble " + addState + "Gate; \n";
		  				  myhfile << "Grid::AttachmentAccessor<Vertex, ADouble> aa" + addState + "Gate; \n";

		  				  //writting of attachment inits for cpp file
		  				  mycppfile << "if (spGrid->has_vertex_attachment(this->" + addState + "Gate)) \n";
		  				  mycppfile << "UG_THROW(\"Attachment necessary (" + addState + "Gate) for " +filename+" channel dynamics \"\n";
		  				  mycppfile << "\"could not be made, since it already exists.\"); \n";
		  				  mycppfile << "spGrid->attach_to_vertices(this->" + addState + "Gate); \n";
		  				  mycppfile << "this->aa" + addState + "Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->" + addState + "Gate); \n \n";
	   */




	  /// Writing function for State-Output
	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "std::vector<number> " + filename + "<TDomain>::state_values(number x, number y, number z) \n";
	  mycppfile << "{ \n";

	  mycppfile << "\t //var for output \n";
	  mycppfile << "\t std::vector<number> GatingAccesors; \n \n";

	  mycppfile << "\t typedef ug::MathVector<TDomain::dim> position_type; \n \n";
	  mycppfile << "\t position_type coord; \n \n";

	  mycppfile << "\t if (coord.size()==1) \n";
	  mycppfile << "\t \t coord[0]=x; \n";

	  mycppfile << "\t if (coord.size()==2) \n";
	  mycppfile << "\t { \n";
	  mycppfile << "\t \t coord[0] = x;\n";
	  mycppfile << "\t \t coord[1] = y;\n";
	  mycppfile << "\t } \n";

	  mycppfile << "\t if (coord.size()==3) \n";
	  mycppfile << "\t { \n";
  	  mycppfile << "\t \t coord[0] = x;\n";
  	  mycppfile << "\t \t coord[1] = y;\n";
  	  mycppfile << "\t \t coord[2] = z;\n";
  	  mycppfile << "\t } \n";

	  mycppfile << "\t //accesors \n";
	  mycppfile << "\t typedef Attachment<position_type> position_attachment_type; \n";
	  mycppfile << "\t typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accesor_type; \n \n";

	  mycppfile << "\t // Definitions for Iteration over all Elements \n";
	  mycppfile << "\t typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; \n";
	  mycppfile << "\t SubsetGroup ssGrp; \n";
	  mycppfile << "\t try { ssGrp = SubsetGroup(m_pCE->approx_space()->domain()->subset_handler(), this->m_vSubset);} \n";
	  mycppfile << "\t UG_CATCH_THROW(\"Subset group creation failed.\"); \n \n";

	  mycppfile << "\t itType iter; \n";
	  mycppfile << "\t number bestDistSq, distSq; \n";
	  mycppfile << "\t Vertex* bestVrt; \n \n";

	  mycppfile << "\t // Iterate only if there is one Gtting needed \n";
	  // Generating if line with all gatting logs
	  std::string Gatingif = "";
	  for (size_t i=0; i<State_vars.size(); i++)
	  {
		  std::cout << "state_sizes: "<< State_vars.size() << " - " << i+1 << std::endl;
		  if (i+1 != State_vars.size())
			  Gatingif = Gatingif + "m_log_" + State_vars[i] + "Gate || ";
		  else
			  Gatingif = Gatingif + "m_log_" + State_vars[i] + "Gate ";
	  }

	  if (Gatingif!="")
	  {
		  mycppfile << "\t if (" << Gatingif << ")\n";
		  mycppfile << "\t { \n";
		  mycppfile << "\t \t // iterating over all elements \n";
		  mycppfile << "\t \t for (size_t si=0; si < ssGrp.size(); si++) \n";
		  mycppfile << "\t \t { \n";
		  mycppfile << "\t \t \t itType iterBegin = m_pCE->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]); \n";
		  mycppfile << "\t \t \t itType iterEnd = m_pCE->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); \n \n";
		  mycppfile << "\t \t \t const position_accesor_type& aaPos = m_pCE->approx_space()->domain()->position_accessor(); \n";
		  mycppfile << "\t \t \t if (si==0) \n";
		  mycppfile << "\t \t \t { \n";
		  mycppfile << "\t \t \t \t bestVrt = *iterBegin; \n";
		  mycppfile << "\t \t \t \t bestDistSq = VecDistanceSq(coord, aaPos[bestVrt]); \n";
		  mycppfile << "\t \t \t } \n";
		  mycppfile << "\t \t \t iter = iterBegin; \n";
		  mycppfile << "\t \t \t iter++; \n";
		  mycppfile << "\t \t \t while(iter != iterEnd) \n";
		  mycppfile << "\t \t \t { \n";
		  mycppfile << "\t \t \t \t distSq = VecDistanceSq(coord, aaPos[*iter]); \n";
		  mycppfile << "\t \t \t \t { \n";
		  mycppfile << "\t \t \t \t \t bestDistSq = distSq; \n";
		  mycppfile << "\t \t \t \t \t bestVrt = *iter; \n";
		  mycppfile << "\t \t \t \t } \n";
		  mycppfile << "\t \t \t \t ++iter; \n";
		  mycppfile << "\t \t \t } \n";
		  mycppfile << "\t \t } \n";

		  // questioning which state should be outputted
		  for (size_t i=0; i<State_vars.size(); i++)
		  {
			  mycppfile << "\t \t if (m_log_" << State_vars[i] << "Gate == true) \n";
			  mycppfile << "\t \t \t GatingAccesors.push_back(this->aa" << State_vars[i] << "Gate[bestVrt]); \n";
		  }
		  mycppfile << "\t } \n";

		  mycppfile << "\t return GatingAccesors; \n";
		  mycppfile << "} \n \n";
	  }
	  else
	  {
		  mycppfile << "\t return GatingAccesors; \n";
		  mycppfile << "} \n \n";
	  }


	  mycppfile << "//Setters for states_outputs \n";
	  // Writing setters for Gates
	  for (size_t i=0; i<State_vars.size(); i++)
	  {
		  mycppfile << "template<typename TDomain> void "<< filename << "<TDomain>::set_log_" << State_vars[i] << "Gate" <<
				   "(bool bLog" << State_vars[i] << "Gate) { m_log_" << State_vars[i] << "Gate = " <<
				   "bLog" << State_vars[i] << "Gate; }\n";
	  }



	  //get information if anyion_flux is needed
	  std::vector<string> Ion_fluxes = Read_i_value(NEURON);


	  // Method for initialization of States with values read out of Initial
	  mycppfile << " // Init Method for using gatings \n";
	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) \n";
	  mycppfile << "{ \n";
	  // TODO rework that with right functions
	  mycppfile << "//get celsius and time\n";
	  mycppfile << "// inits temperatur from kalvin to celsius and some other typical neuron values\n";
	  mycppfile << "number m_T, m_R, m_F; \n";
	  mycppfile << "m_T = m_pCE->temperature(); \n";
	  mycppfile << "m_R = m_pCE->R; \n";
	  mycppfile << "m_F = m_pCE->F; \n \n \n";
	  mycppfile << "number celsius = m_pCE->temperature_celsius(); \n";
	  mycppfile << "number dt = m_pCE->time(); \n";

	  if (Ion_fluxes[0]!="")
	  {
		  for (size_t i = 0; i<Ion_fluxes.size(); i++)
		  {
			  mycppfile << Ion_fluxes[i];
		  }

	  }


	  mycppfile << "// make preparing vor getting values of every edge \n";

	  mycppfile << "number v = vrt_values[CableEquation<TDomain>::_v_]; \n";

	  // TODO we only need read ions in right context
	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "number " + ListIons[i] +" = vrt_values[CableEquation<TDomain>::_"+ ListIons[i] +"_]; \n";
	  }
	  mycppfile << "\n \n";


	  std::cout << "Start working on Initial-Block" << std::endl;

	  // Reading functional info for states
	  std::cout << "Before GetBlock" << std::endl;
	  std::vector<string> INITIAL = GetBlock(Pairs, Zeilen, "INITIAL");
	  std::cout << "After GetBlock" << std::endl;
	  std::cout << "writing proc-Block" << std::endl;
	  std::vector<vector <string> > Proc_funcs = write_proc_block(Pairs, Zeilen);
	  std::cout << "after writing proc-Block" << std::endl;
	  std::vector<string> Proc_vals, locals;
	  string helper;

	  std::vector<string> new_locals, new_local_block;
	  std::vector<string> PROCEDURE = GetBlock(Pairs, Zeilen, "PROCEDURE");
	  size_t begin;
	  size_t Stats_beg = 1;

	  // v always used as var
	  HFile_added_Vars_Init.push_back("v");


	  std::vector<int> stats_begins;

	  std::cout << "before initial size: " << INITIAL.size() << std::endl;

	  if (INITIAL.size() > 0)
	  {
		  begin = count_beg(INITIAL[1]); //original
	  }

	  ////std::cout << "beg: INITIAL!" << std::endl;

	  if (INITIAL.size() > 0)
	  {
		  for (size_t i=0; i<INITIAL.size(); i++)
		  {
			  //std::cout << INITIAL[i] << std::endl;
			  for (size_t j=0; j<Proc_funcs.size(); j++)
			  {
				  begin = count_beg(INITIAL[i]);
				  if (Proc_funcs[j][0] == INITIAL[i].substr(begin, INITIAL[i].find("(")-begin))
				  {
					  // sets beginning of different procedures because we need the lines beetween procedures
					  Stats_beg = i;
					  stats_begins.push_back(Stats_beg);
					  // searching for all vars in function initial i
					  string init_vars = Remove_all(INITIAL[i].substr(INITIAL[i].find("(")+1));
					  init_vars = init_vars.replace(init_vars.find(")"), 1, "");
					  std::cout << "string vars: " << init_vars << std::endl;
					  size_t begin_con = 0;
					  size_t con = init_vars.find(',');
					  std::cout << "first con " << con << std::endl;
					  vector<string> vars_con;

					  while (con != init_vars.npos)
					  {
						  std::cout << "in while" << std::endl;
						  vars_con.push_back(init_vars.substr(begin_con, con-begin_con));
						  std::cout << "substr works: " << init_vars.substr(begin_con, con-begin_con) << std::endl;
						  begin_con = con+1;
						  con = init_vars.find(',', begin_con);
						  std::cout << "end while with begin_con: " << begin_con << std::endl;
						  std::cout << "end while with con: " << con << std::endl;
					  }
					  std::cout << "last push_back" << std::endl;
					  vars_con.push_back(init_vars.substr(begin_con));

					  for (size_t k=1; k<Proc_funcs[j].size(); k++)
					  {
						  // only if its not known global or ion
						  if ((Remove_all(Proc_funcs[j][k])!="ca") && (Remove_all(Proc_funcs[j][k])!="na") && (Remove_all(Proc_funcs[j][k])!="v")
						  && (Remove_all(Proc_funcs[j][k])!="k" ))
						  {
							  bool global = false;
							  for (size_t m=0; m<HFile_added_Vars.size(); m++)
							  {
								  if (Remove_all(HFile_added_Vars[m])==Remove_all(Proc_funcs[j][k]))
									  global = true;
							  }
							  if (global == false)
								  mycppfile << "number " + Proc_funcs[j][k] + " = " + vars_con[k-1] + ";\n";
						  }

					  }

					  // write needed values
					  for (size_t k=1; k<Proc_funcs[j].size(); k++)
					  {
						  // when only read no write read is same as use
						  if (Only_Read(Pairs, Zeilen, Proc_funcs[j][k])==true)
						  {
							  //std::cout << Write_Only_Read(Pairs, Zeilen, Proc_funcs[j][k]) << std::endl;
							  mycppfile << Write_Only_Read(Pairs, Zeilen, Proc_funcs[j][k]) + ";\n";
					  	  }
					  	  else
					  	  {
					  		  // add only if they are in gating or ion list
					  		  string GatingName = In_NeuronUse_List(Pairs, Zeilen, Proc_funcs[j][k]);
					  		  if (GatingName!="")
					  		  {
					  			  //std::cout << "Vielleicht1: " + Remove_all(Proc_funcs[j][k]) + " = aa"+ GatingName + "Gate[vrt]" << std::endl;
					  			  mycppfile << "double " + Remove_all(Proc_funcs[j][k]) + " = "+ GatingName + "; \n \n";
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
							  std::cout << "locals " << locals[k] << std::endl;
							  mycppfile << "double " + locals[k] + "; \n";
							  // removes all unneeded formation-styles
							  Remove_all(locals[k]);
							  // locals need seperation through ","
							  locals_sep = locals[k].find(",");
							  if (locals_sep!=0)
							  {
								  while (locals_sep != locals[k].npos)
								  {
									  HFile_added_Vars_Init.push_back(Remove_all(locals[k].substr(0, locals_sep)));
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
  								  std::cout << "Init: " << HFile_added_Vars_Init[l] << std::endl;
  								  if (HFile_added_Vars_Init.size()>0)
  								  {
  								  	  if (Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")))==Remove_all(HFile_added_Vars_Init[l]) || Check_if(Proc_vals[k]))
  								  		  HFile_added = true;
  								  }
  							  }

  							  string com_buf;

  							  // if it added in hfile write without double else with!
  							  if (HFile_added == true)
  							  {
  								  //writes comments
  								  if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
  									  Proc_vals[k].replace(Proc_vals[k].find(":"), 1, ";//");
  								  //std::cout << Proc_vals[k] << std::endl;
  								  mycppfile << Proc_vals[k] + "; \n";

  							  }
  							  else
  							  {
									// writes comments, if comment for the whole than no double is allowed
									if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
									{
										com_buf = Proc_vals[k];
										com_buf = Remove_all_com(com_buf);
										Proc_vals[k].replace(Proc_vals[k].find(":"), 1, ";//");
									}
									////std::cout << "com_buf: " << com_buf << std::endl;

									if (com_buf.find(":")==0)
									{
										//std::cout << Proc_vals[k] << std::endl;
										mycppfile << Proc_vals[k] + "\n";
									} else
									{
										////std::cout << Proc_vals[k] << std::endl;
										if (Remove_all(Proc_vals[k])!="")
										{
											if (Check_if(Proc_vals[k])==false)
											{
												// if it is a line comment
												if (Proc_vals[k].find("//--//")==Proc_vals[k].npos)
												{
													mycppfile << "double " + Proc_vals[k] + "; \n";
													HFile_added_Vars_Init.push_back(Proc_vals[k].substr(0, Proc_vals[k].find("=")));
												}
												else
												{	// if it is a comment
													mycppfile << Proc_vals[k] + "\n";
												}
											}
										}


									}
  							  }

						  }
					  }



				  }
			  }
		  }		  // write locals vals as double

		  for (size_t i=0; i<stats_begins.size(); i++)
		  {
			  std::cout << "important stats: " << stats_begins[i] << std::endl;
		  }


		  size_t int_end;
		  size_t int_begin;

		  for (size_t m=0; m<stats_begins.size(); m++)
		  {
			  if (m+1<stats_begins.size())
			  {
				  int_begin = stats_begins[m];
				  int_end = stats_begins[m+1];
			  } else
			  {
				  int_begin = stats_begins[m];
				  int_end = INITIAL.size();
			  }

			  std::cout << "int being: " << int_begin << "int_end: " << int_end << std::endl;
			  //std::cout<< "writing end" << std::endl;
			  for (size_t i=int_begin; i<int_end; i++)
			  {
				  std::cout << "i: " << i << std::endl;
				  // writings at the end
				  if (i >= int_begin)
				  {
					  size_t gleich, vorgleich;
					  if ((INITIAL[i]!="") && (INITIAL[i]!="}"))
					  {
						  //std::cout << "gleich Search" << std::endl;
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

										  string left = INITIAL[i].substr(0, INITIAL[i].find("="));
										  //std::cout << left << std::endl;
										  string right = INITIAL[i].substr(INITIAL[i].find("=")+1, INITIAL[i].npos-(INITIAL[i].find("=")+1));
										  //std::cout << ProofSGatingInit(left, SGating) + " = " + right << std::endl;
										  // checks if var is global or needs check
										  bool with_double = true;
										  bool same_as_left = false;

										  string ProofedSGating_Left = ProofSGatingInit(Remove_all(left), SGating);

										  std::cout << "left vs ProofedSGaating_left: " + left + " vs " + ProofedSGating_Left << std::endl;
										  if (Remove_all(left) == Remove_all(ProofedSGating_Left))
										  {
											  same_as_left = true;
										  }

										  for (size_t k=0; k<HFile_added_Vars_Init.size(); k++)
										  {
											  std::cout << "Vergleich: " + Remove_all(ProofedSGating_Left) + " vs " + Remove_all(HFile_added_Vars_Init[k]) << std::endl;
											  if (Remove_all(ProofedSGating_Left) == Remove_all(HFile_added_Vars_Init[k]))
											  {
												  with_double = false;
											  }
											  // if there is a if in
											  if (Check_if(ProofedSGating_Left)==true)
												  with_double = false;


										  }

										  if (with_double == true && same_as_left == true)
										  {
											  mycppfile << "double " + ProofedSGating_Left + " = " + right + "; \n";
											  HFile_added_Vars_Init.push_back(ProofedSGating_Left);
										  }

										  if (with_double == false || same_as_left == false)
										  {
											  mycppfile <<  ProofedSGating_Left + " = " + right + "; \n";
										  }



									  }

									  else
									  {
										  //std::cout << "inittest: " <<INITIAL[i] << std::endl;
										  mycppfile << INITIAL[i] + "; \n";
									  }
								  }
							  }
						  } //else
							  //vorgleich = beg_count(INITIAL[i]);
					  }
				  }
			  }
		  }
	  }


	  std::cout << "all init worked" << std::endl;
	  mycppfile << "}  \n \n \n \n";

	  std::cout << "Initial block written" << std::endl;


////////////////////////////////////////////////////////////////////////////////////////////////
// Starts writting update_gating function
///////////////////////////////////////////////////////////////////////////////////////////////



	  std::cout << "Start writting update_gating function" << std::endl;

	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) \n";
	  mycppfile << "{ \n";
	  // TODO working with right functions from CableEquation
	  mycppfile << "// inits temperatur from kalvin to celsius and some other typical neuron values\n";
	  mycppfile << "number m_T, m_R, m_F; \n";
	  mycppfile << "m_T = m_pCE->temperature(); \n";
	  mycppfile << "m_R = m_pCE->R; \n";
	  mycppfile << "m_F = m_pCE->F; \n \n \n";
	  mycppfile << "number celsius = m_pCE->temperature_celsius(); \n ";
	  mycppfile << "number FARADAY = m_pCE->F; \n ";

	  if (Ion_fluxes[0]!="")
	  {
		  for (size_t i = 0; i<Ion_fluxes.size(); i++)
		  {
			  mycppfile << Ion_fluxes[i];
		  }

	  }

	  mycppfile << "number dt = newTime - m_pCE->time(); \n";
	  mycppfile << "number v = vrt_values[CableEquation<TDomain>::_v_]; \n";

	  for (size_t i=0; i<ListIons.size(); i++)
	  {
		  mycppfile << "number " + ListIons[i] +" = vrt_values[CableEquation<TDomain>::_"+ ListIons[i] +"_]; \n";
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
		  std::cout << "DERIV STATE Vars: double " + State_vars[i] + " = aa" + State_vars[i] + "Gate[vrt]; \n" << std::endl;
		  mycppfile << "double " + State_vars[i] + " = aa" + State_vars[i] + "Gate[vrt]; \n";
		  HFile_added_Vars_Deriv.push_back(Remove_all(State_vars[i]));
	  }

	  // v always used as var
	  HFile_added_Vars_Deriv.push_back("v");

	  //Adding v because v is very often needed and always accesible
	  //mycppfile << "double v  = aavGate[*iter]; \n";

	  mycppfile << "\n \n \n";


	  std::cout << "Working on Derivative-Block starts" << std::endl;




	  std::vector<int> Deriv_begins;


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
	  					 Deriv_begins.push_back(i);
	  					 Stats_beg = i;
	  					 string init_vars = Remove_all(DERIVATIVE[i].substr(DERIVATIVE[i].find("(")+1));
	  					 init_vars = init_vars.replace(init_vars.find(")"), 1, "");
	     				 std::cout << "string vars: " << init_vars << std::endl;
	  					 size_t begin_con = 0;
	  					 size_t con = init_vars.find(',');
	  					 std::cout << "first con " << con << std::endl;
	  					 vector<string> vars_con;

	  					 while (con != init_vars.npos)
	  					 {
	  						  std::cout << "in while" << std::endl;
	  						  vars_con.push_back(init_vars.substr(begin_con, con-begin_con));
	  						  std::cout << "substr works: " << init_vars.substr(begin_con, con-begin_con) << std::endl;
	  						  begin_con = con+1;
	  						  con = init_vars.find(',', begin_con);
	  						  std::cout << "end while with begin_con: " << begin_con << std::endl;
	  						  std::cout << "end while with con: " << con << std::endl;
	  					 }
	  					 std::cout << "last push_back" << std::endl;
	  					 vars_con.push_back(init_vars.substr(begin_con));

	  					 for (size_t k=1; k<Proc_funcs[j].size(); k++)
	  					 {
							  // only if its not known global or ion
							  if ((Remove_all(Proc_funcs[j][k])!="ca") && (Remove_all(Proc_funcs[j][k])!="na") && (Remove_all(Proc_funcs[j][k])!="v")
							  && (Remove_all(Proc_funcs[j][k])!="k" ))
							  {
								  bool global = false;
								  for (size_t m=0; m<HFile_added_Vars.size(); m++)
								  {
									  if (Remove_all(HFile_added_Vars[m])==Remove_all(Proc_funcs[j][k]))
										  global = true;
								  }
								  if (global == false)
									  mycppfile << "number " + Proc_funcs[j][k] + " = " + vars_con[k-1] + ";\n";
							  }
	  					 }
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
						  			  // vielleicht 2
						  			  //std::cout << "Treffer2: "+ Proc_funcs[j][k] + " = aa"+ GatingName + "Gate[vrt];" << std::endl;
						  			  if (Proc_vals[k].find("//--//")==Proc_vals[k].npos)
						  				  mycppfile << "double " + Proc_funcs[j][k] + " = " + GatingName + "; \n \n";
						  			  else
						  				  mycppfile << Proc_funcs[j][k] + " = " + GatingName + "; \n \n";
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
								  //std::cout << "double " + locals[k] << std::endl;
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
	  								  if (((Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")))==Remove_all(HFile_added_Vars_Deriv[l])) || (Check_if(Proc_vals[k])==true)))
	  									  HFile_added = true;
	  							  }

	  							  string com_buf;

	  							  // if it added in hfile write without double else with!
	  							  if (HFile_added == true)
	  							  {
	  								  //writes comments
	  								  if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
	  									  Proc_vals[k].replace(Proc_vals[k].find(":"), 1, ";//");

	  								  mycppfile << Proc_vals[k] + "; \n";
	  								  //std::cout << "Proc_vals[k]: " << Proc_vals[k] << std::endl;

	  							  } else
	  							  {
	  								// writes comments, if comment for the whole than no double is allowed
	  								if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
	  								{
	  									com_buf = Proc_vals[k];
	  									com_buf = Remove_all_com(com_buf);
	  									Proc_vals[k].replace(Proc_vals[k].find(":"), 1, ";//");
	  								}

	  								if (com_buf.find(":")==0)
	  								{
	  									Proc_vals[k].replace(Proc_vals[k].find(";"),1,"");
	  									mycppfile << Proc_vals[k] + "\n";
	  								} else
	  								{
	  									//std::cout<< Proc_vals[k] << std::endl;
	  									if (Remove_all(Proc_vals[k])!="")
	  									{
	  										//std::cout << "double " + Proc_vals[k] << std::endl;
	  										//checking if it is a comment
	  										if (Proc_vals[k].find("//--//")==Proc_vals[k].npos)
	  										{
	  											//std::cout << "write deriv: " << "double " << Proc_vals[k] << std::endl;
	  											mycppfile << "double " + Proc_vals[k] + "; \n";
	  										}
	  										else
	  										{
	  											//std::cout << "write deriv: " << Proc_vals[k] << std::endl;
	  											mycppfile << Proc_vals[k] + "; \n";
	  										}

	  										if (Proc_vals[k].find("=")!=Proc_vals[k].npos)
	  											HFile_added_Vars_Deriv.push_back(Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")-1)));
	  										else
	  											HFile_added_Vars_Deriv.push_back(Remove_all(Proc_vals[k]));
	  									}
	  								}
	  							  }

							  }
	  					  }


	  				  }
	  			  }
	  		  }		  // write locals vals as double
	  	  }

		  size_t deriv_end;
		  size_t deriv_begin;

		  for (size_t m=0; m<Deriv_begins.size(); m++)
		  {
			  if (m+1<Deriv_begins.size())
			  {
				  deriv_begin = Deriv_begins[m];
				  deriv_end = Deriv_begins[m+1];
			  } else
			  {
				  deriv_begin = Deriv_begins[m];
				  deriv_end = DERIVATIVE.size();
			  }

			  std::cout << "deriv begin: " << deriv_begin << "deriv_end: " << deriv_end << std::endl;

			  //std::cout<< "writing end" << std::endl;
			  for (size_t i=deriv_begin; i<deriv_end; i++)
			  {
				  // writings at the end
				  if (i >= (size_t) Deriv_begins[m])
				  {
					  // Writting derivfuncs
					  //std::cout << "before deriv" << std::endl;
					  if ((i > (size_t) Deriv_begins[m]) && (i < DERIVATIVE.size()-1))
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
								{
									//std::cout << "double "+ DERIVATIVE[i] << std::endl;
									if (DERIVATIVE[i].find("//--//")==DERIVATIVE[i].npos)
										mycppfile << "double " + DERIVATIVE[i] + "; \n";
									else
										mycppfile << DERIVATIVE[i] + " \n";

									std::cout << "print first: " << DERIVATIVE[i] << std::endl;

									if (DERIVATIVE[i].find("=")!=DERIVATIVE[i].npos)
										HFile_added_Vars_Deriv.push_back(Remove_all(DERIVATIVE[i].substr(0, DERIVATIVE[i].find("=")-1)));
									else
										HFile_added_Vars_Deriv.push_back(Remove_all(DERIVATIVE[i]));
								}
							}
							else
							{
								size_t oje;
								// working with data out of if_handler
								std::vector<string> if_file = if_handling(i, DERIVATIVE);
								if (if_file.size()>1)
								{
									// only if really some if statements are used
									oje = if_file.size()-1;
									istringstream f(if_file[if_file.size()-1]);
									f >> i;
								}
								else
								{
									// else like every time with known vars
									oje = 0;
									mycppfile << DERIVATIVE[i] + "; \n";
								}

								//writting of files out of if TODO add checking vars used if not declare before
								//if -loop for later use
								std::vector<string> setted_in_if;
								bool same = false;
								for (size_t n = 0; n<oje; n++)
								{
									if ((if_file[n].find("=")!=if_file[n].npos) && (if_file[n].find("if")==if_file[n].npos)
										  && (if_file[n].find("else")==if_file[n].npos))
									{
										for (size_t o=0; o<HFile_added_Vars_Deriv.size(); o++)
										{
											std::cout << "Vergleich: " << Remove_all(if_file[n].substr(0, if_file[n].find("=")-1)) << " - "<< Remove_all(HFile_added_Vars_Deriv[o]) << std::endl;
											if (Remove_all(if_file[n].substr(0, if_file[n].find("=")-1))==Remove_all(HFile_added_Vars_Deriv[o]))
												same = true;
										}

										// if there is not such a var then we need to add
										if (same == false)
										{
											setted_in_if.push_back(Remove_all(if_file[n].substr(0, if_file[n].find("=")-1)));
											HFile_added_Vars_Deriv.push_back(Remove_all(if_file[n].substr(0, if_file[n].find("=")-1)));
										}
									}
								}
								// writting needed double's
								if (setted_in_if.size()>0)
								{
									for (size_t n=0; n<setted_in_if.size(); n++)
									{
										if (Remove_all(setted_in_if[n])!="")
										{
											if (n == 0)
											{
												if (setted_in_if[m].find("//--//")==setted_in_if[m].npos)
													mycppfile << "double " + setted_in_if[m];
												else
													mycppfile << setted_in_if[m];
												std::cout<< "print: " << setted_in_if[m] << std::endl;
											}
											else
											{
												mycppfile << ", " + setted_in_if[m];
											}
										}
									}
								}

								// new zeile
								mycppfile << "; \n \n";

								// writing if-loop
								for (size_t m = 0; m<oje; m++)
								{
									mycppfile << if_file[m];
								}




							}
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
	  std::vector<string> BREAKPOINT = GetBlockFunction(Pairs, Zeilen2, "BREAKPOINT");

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
		  // testing if any Procedure is used in Kinetic block than write it down
		  for (size_t i=0; i<KINETIC.size(); i++)
		  {
			  for (size_t j=0; j<Proc_funcs.size(); j++)
			  {
				  begin = 0;
				  std::cout << "Vergleich: " << KINETIC[i].substr(begin, KINETIC[i].find("("))<< " - " << Proc_funcs[j][0] << std::endl;
				  if (Remove_all(KINETIC[i].substr(begin, KINETIC[i].find("("))) == Proc_funcs[j][0])
				  {
  					  locals = get_local_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);

  					  if (locals.size()>0)
  					  {
						  for (size_t k=0; k<locals.size(); k++)
						  {
							  size_t locals_sep;

							  mycppfile << "double " + locals[k] + "; \n";
							  //std::cout << "double " + locals[k] << std::endl;
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
  									  Proc_vals[k].replace(Proc_vals[k].find(":"), 1, ";//");

  								  mycppfile << Proc_vals[k] + "; \n";
  								  //std::cout << Proc_vals[k] << std::endl;

  							  } else
  							  {
  								// writes comments, if comment for the whole than no double is allowed
  								if (Proc_vals[k].find(":")!=Proc_vals[k].npos)
  								{
  									com_buf = Proc_vals[k];
  									com_buf = Remove_all_com(com_buf);
  									Proc_vals[k].replace(Proc_vals[k].find(":"), 1, ";//");
  								}

  								if (com_buf.find(":")==0)
  								{
  									Proc_vals[k].replace(Proc_vals[k].find(";"),1,"");
  									mycppfile << Proc_vals[k] + "\n";
  								} else
  								{
  									//std::cout<< Proc_vals[k] << std::endl;
  									if (Remove_all(Proc_vals[k])!="")
  									{
  										std::cout << "double " + Proc_vals[k] << std::endl;
  										if (Proc_vals[k].find("//--//")==Proc_vals[k].npos)
  											mycppfile << "double " + Proc_vals[k] + "; \n";
  										else
  											mycppfile << Proc_vals[k] + "; \n";
  										if (Proc_vals[k].find("=")!=Proc_vals[k].npos)
  											HFile_added_Vars_Deriv.push_back(Remove_all(Proc_vals[k].substr(0, Proc_vals[k].find("=")-1)));
  										else
  											HFile_added_Vars_Deriv.push_back(Remove_all(Proc_vals[k]));
  									}
  								}
  							 }

						  }
  					  }

				  }
			  }
		  }

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
						  // if there is no line comment
						  if (KINETIC[j].find("//--//")==KINETIC[j].npos)
						  {
							  vector<string> varsL;
							  string left, right, vars;
							  size_t komma;
							  size_t helpco, helpcc;
							  size_t c_f_o, c_l_c;
							  size_t equa = KINETIC[j].find("=");
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





								  //Fehler
								  //std::cout << left << std::endl;
								  Ausgabe.push_back(left + "+ =  (-" +left +"*" + varsL[0] + "+" +right +"*" + varsL[1] + ")");
								  Ausgabe.push_back(right + "+ = ("+left +"*" + varsL[0] + "+ -"+ right +"*" + varsL[1] + ")");




							  } else
							// if for example some vars getting set before using kinectics
							  {
								  if (equa!=KIN_length)
								  {
									  if (KINETIC[j].find("CONSERVE")==KINETIC[j].npos)
										  mycppfile << "double " + KINETIC[j] + ";\n";
								  }

							  }
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

		  std::vector<string> KineticBreakVars;
		  string varName;
		  bool outputNumber;

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
				  {
					  outputNumber = true;
					  // only need var name not whole term
					  if (Ausgabe2[i].find("=")!=Ausgabe2[i].npos)
						  varName = Remove_all(Ausgabe2[i].substr(0, Ausgabe2[i].find("=")));

					  // Checking if var is global
					  for (size_t j = 0; j<HFile_added_Vars.size(); j++)
					  {
						  std::cout << "global vergleich: " << varName << " - " << Remove_all(HFile_added_Vars[j]) << std::endl;
						  if (varName==Remove_all(HFile_added_Vars[j]))
							  outputNumber=false;
					  }

					  // Checking if var was used before
					  if (KineticBreakVars.size()>0)
					  {

						  for (size_t j = 0; j<KineticBreakVars.size(); j++)
						  {
							  if (varName==KineticBreakVars[j])
								  outputNumber=false;
						  }
					  }

					  // Write in cpp file with or without type declaration depending on outputNumber
					  if (outputNumber==true)
					  {
						  mycppfile << "number " + Ausgabe2[i] + "; \n";
						  KineticBreakVars.push_back(varName);
					  }

					  if (outputNumber==false)
						  mycppfile << Ausgabe2[i] + "; \n";


				  }
			  }
		  }

		  mycppfile << " \n \n \n";


		  // writting Outputs
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
					  for (size_t l=0; l<HFile_added_Vars_Deriv.size(); l++)
				  	  {
						  //std::cout << "Vergleich " << Remove_all(BProc[i].substr(0, BProc[i].find("=")-1)) << " - " << Remove_all(HFile_added_Vars_Deriv[l]) << std::endl;
						  if (Remove_all(BProc[i].substr(0, BProc[i].find("=")-1))==Remove_all(HFile_added_Vars_Deriv[l]))
							  HFile_added = true;
				  	  }

				  	  string com_buf;

				  	  // if it added in hfile write without double else with!
				  	  if (HFile_added == true)
				  	  {
				  		  //writes comments
				  		  if (BProc[i].find(":")!=BProc[i].npos)
				  			  BProc[i].replace(BProc[i].find(":"), 1, ";//");

				  		  mycppfile << BProc[i] + "; \n";

				  	  } else
				  	  {
				  		  // writes comments, if comment for the whole than no double is allowed
				  		  if (BProc[i].find(":")!=BProc[i].npos)
				  		  {
				  			  com_buf = BProc[i];
				  			  com_buf = Remove_all_com(com_buf);
				  			  BProc[i].replace(BProc[i].find(":"), 1, ";//");
				  		  }

				  		  if (com_buf.find(":")==0)
				  		  {
				  			  BProc[i].replace(BProc[i].find(";"),1,"");
				  			  mycppfile << BProc[i] + "\n";
				  		  } else
				  		  {
				  			  //checking if BProc has not been written before
				  			  //std::cout << "double " + BProc[i] << std::endl;
				  			  if (BProc[i].find("=")!=BProc[i].npos)
				  			  {
				  				  HFile_added_Vars_Deriv.push_back(BProc[i].substr(0, BProc[i].find("=")-1));
				  				  if (BProc[i].find("//--//")==BProc[i].npos)
				  					  mycppfile << "double " + BProc[i] + "; \n";
				  				  else
				  					  mycppfile << BProc[i] + "; \n";
				  			  }
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
	  mycppfile << " \n \n \n";



	  std::cout << "update_gating Function is written" << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// write head of ionic flux
//////////////////////////////////////////////////////////////////////////////////////////////////////////

	  std::cout << "Start writting Ionic_current function"<< std::endl;


	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) \n";
	  mycppfile << "{ \n \n";

	  mycppfile << "// inits temperatur from kalvin to celsius and some other typical neuron values\n";
	  mycppfile << "number m_T, m_R, m_F; \n";
	  mycppfile << "m_T = m_pCE->temperature(); \n";
	  mycppfile << "m_R = m_pCE->R; \n";
	  mycppfile << "m_F = m_pCE->F; \n \n \n";

	  // writing needed fluxes if any needed
	  if (Ion_fluxes[0]!="")
	  {
		  for (size_t i = 0; i<Ion_fluxes.size(); i++)
		  {
			  mycppfile << Ion_fluxes[i];
		  }
	  }


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
		 	  mycppfile << "number " + ListIons[i] + " = vrt_values[m_pCE->_" + ListIons[i] + "_]; \n";
		  }
	  }



	  // v needed every time
	  mycppfile << "number v =  vrt_values[m_pCE->_v_]; \n";
	  mycppfile << " \n";
	  mycppfile << " \n";
	  // adding current time as t
	  mycppfile << "number t = m_pCE->time(); \n \n \n";


	  std::cout << "in Ions! " << std::endl;

	  std::vector<string> IonsWrittenCurrent;

	  // Adding all needed inside concentrations if some needed ki/nai/cai
	  for (size_t i = 0; i<NEURON.size(); i++)
	  {
		  size_t Useion = NEURON[i].find("USEION");
		  size_t Read = NEURON[i].find("READ");
		  if (Useion!=NEURON[i].npos)
		  {
			  if (Read!=NEURON[i].npos)
			  {
				  string read = NEURON[i].substr(Read+4);
				  // if there is any , or WRITE make string shorter
				  if (read.find(",")!=read.npos)
					  read = read.substr(0, read.find(","));
				  if (read.find("WRITE")!=read.npos)
					  read = read.substr(0, read.find("WRITE"));

				  // for being ki/nai/cai setting to right concenktration k/na/ca
				  read = Remove_all(read);
				  if (read.find("ki")!=read.npos)
				  {
					  mycppfile << "number " << read << " = " << read.substr(0, read.find("i")) << ";\n";
					  IonsWrittenCurrent.push_back(read);
				  }
				  if (read.find("nai")!=read.npos)
				  {
					  mycppfile << "number " << read << " = " << read.substr(0, read.find("i")) << ";\n";
					  IonsWrittenCurrent.push_back(read);
				  }
				  if (read.find("cai")!=read.npos)
				  {
					  mycppfile << "number " << read << " = " << read.substr(0, read.find("i")) << ";\n";
					  IonsWrittenCurrent.push_back(read);
				  }

			  }
		  }
		  std::cout << "in Neuronline: " << i << std::endl;

	  }

	  mycppfile << "\n \n \n";

	  std::cout << "after Neuronline" << std::endl;





	  //writes Nernst Eqs
	  std::vector<string> eqs = equali(Pairs, Zeilen);

	  //removing wrong ION
	  if (eqs.size()>0)
	  {
		  for (size_t i=0; i<eqs.size(); i++)
		  {
			  if (eqs[i].find("ION")!=eqs[i].npos)
				  eqs[i].replace(eqs[i].find("ION"),3, "");
		  }
	  }


	  //var vor latest use
	  std::vector<string> comp_outs;
	  string comp_outs_left, comp_outs_right;

	  std::cout << "Start preparing Nernst-Equations" << std::endl;

	  if (eqs.size()>0)
	  {
		  for (size_t i=0; i<eqs.size(); i++)
		  {
			  bool writtenIonBefore = false;
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
					  // not only there have to be ica/ina/ik also there has to be next no letter
					  if (((comp_outs_left.find("ica")!=comp_outs_left.npos) && ((is_single_letter(comp_outs_left[(comp_outs_left.find("ica"))+3]))==false))
					  || ((comp_outs_left.find("ina")!=comp_outs_left.npos) && ((is_single_letter(comp_outs_left[(comp_outs_left.find("ina"))+3]))==false))
					  || ((comp_outs_left.find("ik")!= comp_outs_left.npos) && ((is_single_letter(comp_outs_left[(comp_outs_left.find("ik"))+2]))==false)))
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
					  std::cout << "eqs: " << eqs[i] << std::endl;
					  std::cout << "verg2: " << verg << std::endl;
					  if (verg.find("number")!=verg.npos)
						  verg = verg.replace(verg.find("number"), 6, "");
					  //std::cout << "Vergleich: " << Remove_all(verg) << " vs " << Remove_all(Parameter[j]) << std::endl;
					  if (Remove_all(verg) == Remove_all(Parameter[j]))
						  rev_pot_is_setted = true;
				  }

				  if (is_comp_out == false && rev_pot_is_setted == false)
				  {
					  string test = eqs[i];
					  if (test.find("number") != test.npos)
						  test = test.replace(test.find("number"), 6, "");
					  //writting in if-schleife if there is hard coded ena, eca or ek given use this instead of calculation
					  if (i==0)
					  {
						  mycppfile << eqs[i] + " \n";
					  }
					  // if really only eca ek or ena is meant
					  else if (Remove_all(test)[0]!='e')
					  {
						  if (IonsWrittenCurrent.size()>0)
						  {
							  //only if not added before

							  for (size_t l=0; l<IonsWrittenCurrent.size(); l++)
							  {
								  //std::cout << "neuvergleich: " << Remove_all(test.substr(0, test.find("="))) << " - " << (IonsWrittenCurrent[l]) << std::endl;
								  if (Remove_all(test.substr(0, test.find("=")))==(IonsWrittenCurrent[l]))
									  writtenIonBefore = true;
							  }

						  }
						  std::cout << "bool value: " << writtenIonBefore << std::endl;

						  if (writtenIonBefore==false)
							  mycppfile << eqs[i] + " \n";
					  }
					  else
					  {
						  string e_name2;
						  string e_name = eqs[i].substr(0, eqs[i].find("=")-1);
						  mycppfile << e_name + "; \n";
						  eqs[i] = eqs[i].replace(eqs[i].find("number"), 6, "");
						  e_name2 = e_name.replace(e_name.find("number"), 6, "");
						  //std::cout << "equis: " << eqs[i] << std::endl;
						  mycppfile << "if (m_pCE->" + Remove_all(e_name2) + "() == 0) \n";
						  mycppfile << "{ \n";
						  mycppfile << "\t " + eqs[i] + " \n";
						  mycppfile << "} \n";
						  mycppfile << "else \n";
						  mycppfile << "{ \n";
						  mycppfile << "\t " + e_name + " = m_pCE->" + Remove_all(e_name2) + "(); \n";
						  mycppfile << "} \n";
					  }
				  }
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

	  if (BREAKPOINT.size()>0)
	  {
		  //remove local vals
		  for (size_t i=0; i<BREAKPOINT.size(); i++)
		  {
			  if (BREAKPOINT[i].find("LOCAL")!=BREAKPOINT[i].npos)
				  BREAKPOINT[i]="";
		  }

		  //decids which Breakpoint Method will be used
		  for (size_t i=1; i<BREAKPOINT.size()-1; i++)
		  {
			  std::cout << "Breakpoint number: " << i << std::endl;
			  if (BREAKPOINT[i].find("if")!=BREAKPOINT[i].npos)
			  {
				  stand_BP = false;
			  }
		  }

		  std::cout << "end error" << std::endl;

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
					  std::cout << "NONSPEC_CURRENT: " << non_spec_current;
				  }
			  }

			  string test;
			  for (size_t i=1; i<BREAKPOINT.size()-1; i++)
			  {

				  if (BREAKPOINT[i].find("if")!=BREAKPOINT[i].npos)
				  {
					  std::vector<string> break_if = if_handling(i, BREAKPOINT);

					  //writting of files out of if TODO add checking vars used if not declare before
					  //if -loop for later use
					  std::vector<string> setted_in_if;
					  bool same = false;
					  for (size_t m = 0; m<break_if.size()-1; m++)
					  {
						  if ((break_if[m].find("=")!=break_if[m].npos) && (break_if[m].find("if")==break_if[m].npos)
								  && (break_if[m].find("else")==break_if[m].npos))
						  {
							  for (size_t n=0; n<HFile_added_Vars_org.size(); n++)
							  {
								  if (Remove_all(break_if[m].substr(0, break_if[m].find("=")-1))==Remove_all(HFile_added_Vars_org[n]))
									  same = true;
							  }

							  // if there is not such a var then we need to add
							  if (same == false)
							  {
								  setted_in_if.push_back(Remove_all(break_if[m].substr(0, break_if[m].find("=")-1)));
								  HFile_added_Vars_org.push_back(Remove_all(break_if[m].substr(0, break_if[m].find("=")-1)));
							  }
						  }
					  }
					  // writting needed double's
					  if (setted_in_if.size()>0)
					  {
						  for (size_t m=0; m<setted_in_if.size(); m++)
						  {
							  if (Remove_all(setted_in_if[m])!="")
							  {
								  if (m == 0)
									  mycppfile << "double " + setted_in_if[m];
								  else
									  mycppfile << ", " + setted_in_if[m];
							  }
						  }
					  }

					  // new zeile
					  mycppfile << "; \n \n";

					  // writing if-loop
					  for (size_t m = 0; m<break_if.size()-1; m++)
					  {
						  mycppfile << break_if[m];
					  }





					  // writing position to work on
					  istringstream f(break_if[break_if.size()-1]);
					  f >> i;
				  }

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
			  if (Remove_all(non_spec_current)=="")
				  non_spec_current="0";

			  mycppfile << "outCurrentValues.push_back(" + non_spec_current + "); \n";

		  }



		  vector<string> outsHelp;
		  vector<bool> RealWrittenOuts;
		  vector<string> Breakpoint_vars;

		  vector<string> BreakpointOutput;

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
					  //std::cout << "B_vars: " << (B_vars[i]) << std::endl;
				  }
				  // check if vars on right side an write down
				  for (size_t i=0; i<fluxes.size(); i++)
				  {
					  for (size_t j=0; j<B_vars.size(); j++)
					  {
						  beg = fluxes[i].find("=");
						  std::cout << "fluxi: " << fluxes[i] << " = bei: "<< beg <<std::endl;
						  std::cout << "B_vars: " << (B_vars[j]) << std::endl;
						  if (fluxes[i].find(B_vars[j], beg) != fluxes[i].npos)
						  {
							  std::cout << "var is on right side" << std::endl;
							  size_t test = fluxes[i].find(B_vars[j], beg);
							  if (pos_letterb(fluxes[i].substr(test+B_vars[j].size(), 1))==true)
							  {
								  std::cout << "not only part" << std::endl;
								  //here we need part of B_vars as flux ranking are the same so we can use j flux instead of i
								  //std::cout << "numberfluxes[j] " + fluxes[i] + "; \n" << std::endl;
									 string test = Remove_all(fluxes[j]);
									 if (test[0]!='i' && test[1]!='f' && test[2]!=' ')
								  {
									  BreakpointOutput.push_back("number " + fluxes[j] + "; \n");
									  fluxes[j] = "";
									  B_vars[j] = "ddddddddddddd";
								  }
								  else
								  {
									  std::cout << "else fluxes: " << fluxes[j] << std::endl;
									  // checking if var is in global list or Breakpoint list
									  bool global = false;
									  for (size_t m=0; m<Global_vars.size(); m++)
									  {
										  if (fluxes[j].find("=")!=fluxes[j].npos)
										  {
											  if (fluxes[j].substr(fluxes[j].find("="))==Global_vars[m])
											  {
												  global = true;
											  }
										  }
									  }
									  if (global == true)
									  {
										  BreakpointOutput.push_back(fluxes[j] + "; \n");
									  }
									  else
									  {
										  BreakpointOutput.push_back("number " + fluxes[j] + "; \n");
									  }
									  fluxes[j] = "";
									  B_vars[j] = "ddddddddddddd";
								  }
							  }
						  }

					  }
				  }

				  vector<bool> Ion_outside;
				  vector<string> outs;
				  // first every time v
				  outs.push_back("");

				  int vm_flux_count = 0;

				  BreakpointOutput.push_back("\n \n \n");

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
						  else
						  {
							 // writes all parts out
							 string test = Remove_all(fluxes[i]);
							 if (test[0]!='i' && test[1]!='f' && test[2]!=' ')
								 BreakpointOutput.push_back("number " + fluxes[i] + "; \n");
							 else
								 BreakpointOutput.push_back(fluxes[i] + "; \n");
							 //std::cout << "fluxesout[i] " << fluxes[i] << std::endl;
						  }
					  }
				  }




				  for (size_t i=0; i<fluxes.size(); i++)
				  {
					  for (size_t j=0; j<ListIons.size(); j++)
					  {
						  std::cout << "ListIons: " << ListIons[j] << std::endl;
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

									  //setting Ion_outside true to prevent output later as push_back from current
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
												 //std::cout << "komplikated: " <<eqs[k+1] + fluxes[i].substr(fluxes[i].find("="), fluxes[i].npos-fluxes[i].find("=")) + "; \n" << std::endl;
												 mycppfile << eqs[k+1] + fluxes[i].substr(fluxes[i].find("="), fluxes[i].npos-fluxes[i].find("=")) + "; \n";
											 }

											 if (fluxes[i].find(right)!=fluxes[i].npos)
											 {
												 mycppfile << "number " + eqs[k+1] + "o = m_pCE->" + eqs[k+1] + "_out(); \n";
												 mycppfile << "\n \n \n";
											 }

										  }
									  }

								  }
							  }


						  }
						  std::cout << "fluxes!!! " << fluxes[i] << "!!!!!!!!!!!" <<std::endl;
						  // adding always in outs but only if there is afterwards no other letter
						  if ((fluxes[i].find("i"+ListIons[j])!=fluxes[i].npos)
						  && (is_single_letter(fluxes[i][fluxes[i].find("i"+ListIons[j]) + ListIons[j].size()+1])==false))
						  {
							  outs.push_back(fluxes[i].substr(fluxes[i].find("=")+1, fluxes[i].npos -fluxes[i].find("=")+1));
							  outsHelp.push_back(ListIons[j]);
						  }
					  }
				  }

				  ///////////writing outputs
				  for (size_t t=0; t<BreakpointOutput.size(); t++)
				  {
					  mycppfile << BreakpointOutput[t];
				  }



			//search for non spec_ionic_current and set it if used
			for (size_t i=0; i<NEURON.size(); i++)
			{
			non_spec_beg = NEURON[i].find("NONSPECIFIC_CURRENT");
			if (non_spec_beg!=NEURON[i].npos)
				{
					non_spec_end = NEURON[i].find(" ", non_spec_beg+20);
					non_spec_current = NEURON[i].substr(non_spec_beg+20, non_spec_end-(non_spec_beg+20));
					std::cout << "NONSPEC_CURRENT: " << non_spec_current;
				}
			}



		 	  // if output is empty write 0 defekt
		 	  if ((Remove_all(outs[0])=="") || (outs[0]==" "))
		 		  outs[0]="0";

		 	  if (non_spec_current!="")
		 		  mycppfile << "outCurrentValues.push_back(" + outs[0] + " + " + non_spec_current + "); \n";
		 	  else
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
		 		 bool m_F_needed_2 = false;
		 		 for (size_t j=0; j<All_Eqs.size(); j++)
		 		 {
		 			 // if Nernst Pot used we need to add /m_F
		 			 if (outs[i].find(Remove_all(All_Eqs[j]))!=outs[i].npos)
		 			 {
		 				 //calcium hast valence of two so we need 2*m_F
		 				 m_F_needed = true;

		 				 //if there is ca
		 				 if (Remove_all(All_Eqs[j]).find("ca")!= Remove_all(All_Eqs[j]).npos)
		 						 m_F_needed_2 = true;
		 			 }
		 		 }
		 		 //std::cout << "whats in outs: " << outs[i] << " ahh thats in it" << std::endl;

		 		 //std::cout << "optional: " << optional << " force " << force << std::endl;
		 		 if (m_F_needed==true && (optional==true || force == true))
		 		 {
		 			 if (m_F_needed_2 == true)
		 			 {
		 				 mycppfile << "outCurrentValues.push_back(" + outs[i] + "/2 * m_F ); \n";
		 				 RealWrittenOuts.push_back(true);
		 			 }
		 			 else
		 			 {
		 				 mycppfile << "outCurrentValues.push_back(" + outs[i] + "/m_F ); \n";
		 				 RealWrittenOuts.push_back(true);
		 			 }
		 		 }

		 		 if (m_F_needed==false && (optional==true || force == true))
				 {
		 			 mycppfile << "outCurrentValues.push_back(" + outs[i] + " ); \n";
		 			 RealWrittenOuts.push_back(true);
				 }
		 	  }

	  }
	  mycppfile << "} \n \n \n";

	  std::cout << "Ionic_current Function is written" << std::endl;

	  std::cout << "writting specify_write_functions_indices()..." << std::endl;

	  mycppfile << "template<typename TDomain> \n";
	  mycppfile << "void " + filename + "<TDomain>::specify_write_function_indices() \n";
	  mycppfile << "{ \n \n";
	  // First value is always for v reserved
	  mycppfile << "this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_); \n";
	  // All other outs have to be checked
 	  for (size_t i=0; i<writes_Ions.size(); i++)
 	  {
 		 for (size_t j=0; j<outsHelp.size(); j++)
 		 {
 			 if (Remove_all(writes_Ions[i])==("i"+Remove_all(outsHelp[j])))
 			 {
 				if (RealWrittenOuts.size()>0)
 				{
 					if (RealWrittenOuts[i]==true)
 						mycppfile << "this->m_vWFctInd.push_back(CableEquation<TDomain>::_"<< writes_Ions[i].replace(writes_Ions[i].find("i"),1,"") <<"_); \n";
 				}
 			 }
 		 }
 	  }
	  mycppfile << "} \n \n \n";

	  std::cout << "specify_write_functions_indices() is written!" << std::endl;


	  }





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


	  mycppfile << "} // namespace cable_neuron\n";
	  mycppfile << "} // namespace ug\n\n\n";

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

		  // Before writting all Global vars check for ,
		  for (size_t i = 0; i<Global_vars.size(); i++)
		  {
			 size_t komma = Global_vars[i].find(",");
			 while (komma!=Global_vars[i].npos)
			 {
				 Global_vars.push_back(Global_vars[i].substr(0, komma));
				 Global_vars[i] = Global_vars[i].replace(0, komma+1, "");
				 komma = Global_vars[i].find(",");
				 //std::cout << komma << std::endl;
			 }
		  }
	  }


	  std::vector<string> already_added;


	  if (Params_Unit.size() >0)
	  {
		  for (size_t i = 0; i<Params_Unit.size(); i++)
		  {
			  if (Params_Unit[i].find("//--//")==Params_Unit[i].npos)
			  {
				  myhfile << "number " + Params_Unit[i].substr(0, Params_Unit[i].find("=")) + "; \n";
				  already_added.push_back(Remove_all(Params_Unit[i].substr(0, Params_Unit[i].find("="))));
			  }
		  }
	  }

	  if (Global_vars.size() > 0)
	  {
		  for (size_t i = 0; i<Global_vars.size(); i++)
		  {
			  bool adding = true;
			  for (size_t j =0; j<already_added.size(); j++)
			  {
				  if (Remove_all(already_added[j]) == Remove_all(Global_vars[i]))
					  adding = false;
			  }

			  if (Global_vars[i]!="" && adding==true)
			  {
				  if (Global_vars[i].find("//--//")==Global_vars[i].npos)
				  {
					  std::cout << "Global vars: " << Global_vars[i] << std::endl;
					  myhfile << "number " + Remove_all(Global_vars[i]) + "; \n";
					  already_added.push_back(Global_vars[i]);
				  }
			  }
		  }
	  }

	  for (size_t i=0; i<State_vars.size(); i++)
	  {
		  myhfile << "bool m_log_" << State_vars[i] << "Gate; \n";
	  }

	  //Writing Neuron constant files out of const block
	  std::vector<string> CONSTANT = GetBlock(Pairs, Zeilen, "CONSTANT");
	  if (CONSTANT.size() >= 1)
	  {
		  for (size_t i = 0; i<CONSTANT.size(); i++)
		  {
			  size_t klammer = CONSTANT[i].find("(");
			  size_t eckklammer = CONSTANT[i].find("<");
			  if (klammer==CONSTANT[i].npos and eckklammer==CONSTANT[i].npos and CONSTANT[i].find("{")==CONSTANT[i].npos
			    and CONSTANT[i].find("}")==CONSTANT[i].npos and Remove_all(CONSTANT[i])!="")
				  myhfile << "const number " << CONSTANT[i] << "; \n";
			  else
			  {
				  if (klammer < eckklammer)
				  {
					  if (CONSTANT[i].find("{")==CONSTANT[i].npos and CONSTANT[i].find("}")==CONSTANT[i].npos and Remove_all(CONSTANT[i])!="")
						  myhfile << "const number " << CONSTANT[i].substr(0, klammer) << "; \n";
				  }
				  else
				  {
					  if (CONSTANT[i].find("{")==CONSTANT[i].npos and CONSTANT[i].find("}")==CONSTANT[i].npos and Remove_all(CONSTANT[i])!="")
						  myhfile << "const number " << CONSTANT[i].substr(0, eckklammer) << "; \n";
				  }
			  }
		  }
	  }


	  //Standard Neuron params
	  myhfile << "// Standard-NModl-File-Params \n";
	  myhfile << "number F, R, K, celsius; \n";
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  // Getting some initial values out of Neuron as const numbers
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  std::vector<bool> is_state;

	  std::cout << "before first" << std::endl;

	  std::cout << "State_vars: " << State_vars.size() << " INITIAL: " << INITIAL.size() << std::endl;

	  // initialising is_state vars
	  for (size_t j=0; j<INITIAL.size(); j++)
	  {
		  is_state.push_back(false);
	  }

	  //making bool list with only states numbers are true
	  for (size_t i=0; i<State_vars.size(); i++)
	  {
		  for (size_t j=0; j<INITIAL.size(); j++)
		  {
			  size_t gleich = INITIAL[j].find("=");
			  if (gleich != INITIAL[j].npos)
			  {
				  if (State_vars[i]==Remove_all(INITIAL[j].substr(0, gleich)))
					  is_state[j]=true;
			  }
		  }
	  }

	  std::cout << "first worked" << std::endl;

	  // if hase not been added already
	  for (size_t i=0; i<INITIAL.size(); i++)
	  {
		  size_t gleich = INITIAL[i].find("=");
		  if (gleich!=INITIAL[i].npos)
		  {
			  std::cout << "wokring... " << std::endl;
			  if (is_state[i]==false)
			  {
				  bool test = false;

				  for (size_t j=0; j<already_added.size(); j++)
				  {
					  if (Remove_all(INITIAL[i].substr(0, gleich))==already_added[j])
						  test=true;
				  }
				  if (test == false)
				  {
					  if (INITIAL[i].find("if")==INITIAL[i].npos)
						  myhfile << "const double " << INITIAL[i] << "; \n";
				  }
			  }
		  }
	  }






	  myhfile << "}; \n \n";

	  myhfile << "} // namespace cable_neuron\n";
	  myhfile << "} // namespace ug\n\n\n";

	  myhfile << "#endif // " + filename + "_H_\n";

	  myhfile.close();
}


////////////////////////////////////////////////////////////////////////////////////////
// Part which is needed for including in UG4 directly
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
// Returns all needed Sets for one function
////////////////////////////////////////////////////////////////////////////////////////

std::vector<string> Converter::Get_Set_Registry(string filename)
{
	std::vector<string> erg;
	std::vector<string> lines;

	lines = Openfile(filename);

	for (size_t i=0; i<lines.size(); i++)
	{
		size_t sets = lines[i].find("sets");
		size_t set = lines[i].find("set");
		if (sets==lines[i].npos)
		{
			//.add_method("set_synapse_provider_factory", &T::set_synapse_provider_factory)
			if (set!=lines[i].npos)
			{
				string func_name = lines[i].substr(set, lines[i].find("(")-set);
				erg.push_back(".add_method(\"" + func_name + "\" , &T::" + func_name + ")" );
			}
		}

	}



	return erg;
}


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


	std::vector<string> setters;

	const char* filenamechar = filename.c_str();
	ofstream mychannelfile;
	// if not existing write channel in file
	if (write_channel == true)
	{
		mychannelfile.open (filenamechar, std::ios::app);
		mychannelfile << "#include \"" + Ch_Name + ".h\" \n \n \n";
		mychannelfile << "{ \n";
		mychannelfile << "\t typedef " + Ch_Name + "<TDomain> T; \n";
		mychannelfile << "\t typedef ICableMembraneTransport<TDomain> TBase; \n";
		mychannelfile << "\t string name = string(\""+ Ch_Name + "\").append(suffix); \n";
		mychannelfile << "\t reg.add_class_<T, TBase >(name, grp) \n";
		mychannelfile << "\t \t .template add_constructor<void (*)(const char*, const char*)>(\"Function(s)#Subset(s)\") \n";
		mychannelfile << "\t \t .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>(\"Function(s)#Subset(s)\") \n";

		//ToDO write log_states methodes


		// Adding getters and setter
		setters = Get_Set_Registry(Ch_Name+".h");
		for (size_t i=0; i<setters.size(); i++)
		{
			mychannelfile << setters[i] + "\n";
		}



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
			 name_end = Includes[i].find(".");
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

	 // del all double ..
	 for (size_t i=0; i<SourceFileV.size(); i++)
	 {
		 size_t pointpoint = SourceFileV[i].find("..");
	 	 if (pointpoint!=SourceFileV[i].npos)
	 		 SourceFileV[i].replace(pointpoint, 1, "");
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
					 myDestFile <<  "membrane_transport/nmodl_converter/converted/" + SourceFileV[j] +"\n";
				 }
			 }
	 	 }
	 }
}
