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

	double Unit_Conv(string s);

	std::vector<string> Openfile(string filename);

	std::vector<std::pair<int, int> > FindBlocks(std::vector<string> Zeilen);

	std::vector<string> GetBlock(std::vector<pair<int, int> >, std::vector<string> Zeilen, string name);
	size_t find_beg(string beg);

	size_t number_(size_t pos, string s);
	size_t pos_letter(string s);
	string build_func_head(string s);
	string func_head(string s);

	bool pos_letterb(string s);

	std::vector<string> build_proc_vars(string s);
	std::vector<vector<string> > write_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);
	size_t count_beg(string beg);
	std::vector<string> writer_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name);

	std::vector<string> get_local_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name);
	std::vector<vector<string> > write_derivative_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);

	std::vector<string> equali(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);

	void WriteStart(string filename, std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen);



private:

	string m_class_name;

};

#endif /* CONVERTER_H_ */
