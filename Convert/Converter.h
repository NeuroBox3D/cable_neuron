/*
 * Converter.h
 *
 *  Created on: 04.02.2015
 *      Author: unheiliger
 */


#include <vector>
#include <string>


#ifndef CONVERTER_H_
#define CONVERTER_H_

using namespace std;

class Converter {
public:
	Converter();
	virtual ~Converter();

	// Unit Convert functions
	double Unit_Conv_Value(string s);
	double Unit_Conv_All(string s);

	// Proofs if string is Gating string and give Gating-name for that var
	// Different versions for init and update_gating are needed
	string ProofSGatingInit(string s, std::vector<string> SGating);
	string ProofSGatingGating(string s, std::vector<string> SGating);

	// finding READ and WRITE values behind USEION
	std::vector<string> Find_all_read(std::vector<string> Neuron);
	std::vector<string> Find_all_write(std::vector<string> Neuron);
	std::vector<string> Find_all_Eqs(std::vector<string> Neuron);
	std::vector<string> Read_i_value(std::vector<string> Neuron);

	// functions for removing all unused chars and comments
	vector<string> Remove_all(vector<string> erg);
	string Remove_all(string erg);
	string Remove_all_com(string erg);

	// gives ion name out
	string In_NeuronUse_List(std::vector<pair<int, int> > Pairs, std::vector<string>, string s);


/////////////////////////////////////////////////////////////////////////////////
	// Different Search/Find/Decide Functions
/////////////////////////////////////////////////////////////////////////////////
	// functions testing if char is letter/number
	bool is_single_letter(char s);
	bool is_single_number(char s);

	// testing functions
	bool pos_letterb(string s);
	bool begG(string s);

	// checking for if-loop
	bool Check_if(string test);

	// some different Search Beginnings
	size_t find_beg(string beg);
	size_t find_begg(string beg);
	size_t number_(size_t pos, string s);
	size_t pos_letter(string s);
	size_t count_beg(string beg);
	string writing_starts(string s);

/////////////////////////////////////////////////////////////////////////////////
	// Opensfile-Functions
/////////////////////////////////////////////////////////////////////////////////

	// Open mod file
	std::vector<string> Openfile(string filename);

	std::vector<string> GetFuncTable(std::vector<string> Zeilen);

	// writes file for later including in plugin main
	std::vector<string> WriteChannelFile(string Ch_Name, string filename);
	void WriteInclude_List(std::vector<string> Includes, string filename, string sources);
	void WriteInPlugin(string SourceFile, string DestFile);
	std::vector<string> Get_Set_Registry(string filename);

	// Find all Blocks of NModl-File
	std::vector<std::pair<int, int> > FindBlocks(std::vector<string> Zeilen);

	// Get Blocks like "NEURON"-Block form openend NModl-File
	std::vector<string> GetBlock(std::vector<pair<int, int> >, std::vector<string> Zeilen, string name);
	std::vector<string> GetBlockFunction(std::vector<pair<int, int> >, std::vector<string> Zeilen, string name);


	//testing if on ion is only read needed
	bool Only_Read(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s);
	//writing the only_read Part of ion
	string Write_Only_Read(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s);





	//building Function heads
	string build_func_head(string s);
	string func_head(string s);


	//functions for writting setters and getters for one var
	// Usage first line setter for h-file
	// second and following lines for cpp-file line for line
	std::vector<string> buildgetter(string var, string class_name);
	std::vector<string> buildsetter(string var, string class_name);


	std::vector<string> if_handling(size_t begin, std::vector<string> Unit);


	// Search some needed Procedure out of "PROCEDURE" Block
	string Search_for_Proc(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, std::vector<string> name);

	// Building Vars out of Procfunctionhead
	std::vector<string> build_proc_vars(string s);
	std::vector<vector<string> > write_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);

	std::vector<string> writer_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name);
	std::vector<string> get_local_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name);


	std::vector<vector<string> > write_derivative_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);



	std::vector<string> equali(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);
	std::vector<string> GetProcEqualString(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s);


/// Main Function using all stuff
	void WriteStart(string filename, std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, bool force);



private:

	string m_class_name;

};

#endif /* CONVERTER_H_ */
