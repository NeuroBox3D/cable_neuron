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

	double Unit_Conv_Value(string s);
	double Unit_Conv_All(string s);


	// functions for removing all unused chars and comments
	vector<string> Remove_all(vector<string> erg);
	string Remove_all(string erg);
	string Remove_all_com(string erg);

	// gives ion name out
	string In_NeuronUse_List(std::vector<pair<int, int> > Pairs, std::vector<string>, string s);

	string Convert_Line(std::vector<string> Known_vars, string s);

	// functions testing if char is letter/number
	bool is_single_letter(char s);
	bool is_single_number(char s);

	std::vector<string> Openfile(string filename);

	std::vector<std::pair<int, int> > FindBlocks(std::vector<string> Zeilen);

	bool Only_Read(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s);
	string Write_Only_Read(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s);


	std::vector<string> GetBlock(std::vector<pair<int, int> >, std::vector<string> Zeilen, string name);
	size_t find_beg(string beg);
	size_t find_begg(string beg);
	size_t number_(size_t pos, string s);
	size_t pos_letter(string s);
	string build_func_head(string s);
	string func_head(string s);

	bool pos_letterb(string s);

	bool begG(string s);

	string writing_starts(string s);

	string Search_for_Proc(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, std::vector<string> name);

	std::vector<string> build_proc_vars(string s);
	std::vector<vector<string> > write_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);
	size_t count_beg(string beg);
	std::vector<string> writer_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name);

	std::vector<string> get_local_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name);
	std::vector<vector<string> > write_derivative_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);

	std::vector<string> equali(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);

	std::vector<string> GetProcEqualString(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string s);

	void WriteStart(string filename, std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);



private:

	string m_class_name;

};

#endif /* CONVERTER_H_ */
