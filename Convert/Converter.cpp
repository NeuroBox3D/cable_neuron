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




Converter::Converter() {
	// TODO Auto-generated constructor stub

}

Converter::~Converter() {
	// TODO Auto-generated destructor stub
}

double Converter::Unit_Conv(string s)
{
	// all units have to be in ms m mV
}


std::vector<string> Converter::equali(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen)
{
	std::vector<string> out;
	out.push_back("const number helpV = (m_R*m_T)/m_F;");

	size_t Ion, IonRead, IonRend;
	string IonS;
	std::vector<string> ListIonRead;
	std::vector<string> ListIon;

	std::vector<string> NEURON = GetBlock(Pairs, Zeilen, "NEURON");

	for (size_t i=0; i<NEURON.size(); i++)
	{
		Ion = NEURON[i].find("USEION");
		if (Ion!=NEURON[i].npos)
		{
			IonRead = NEURON[i].find("READ", Ion+7);
			IonRend = NEURON[i].find(" ", IonRead+5);
			std::cout << IonRend << ", " << Ion << std::endl;
			IonS = NEURON[i].substr(IonRead+5, (IonRend-(IonRead+5)));
			std::cout << "Subion: "<< IonS << std::endl;
			ListIonRead.push_back(IonS);
			ListIon.push_back(IonS.substr(1, IonS.size()-1));
		}
	}

	for (size_t i=0; i<ListIon.size(); i++)
	{
		out.push_back("number " + ListIonRead[i] + " = helpV*(log(" + ListIon[i] + "_out/" + ListIon[i] + "));");
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







std::vector<string> Converter::writer_proc_block(std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen, string name)
{
	std::vector<string> fileInfo;
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
				if (PROCEDURE[j].find("UNITSOFF")!=PROCEDURE[j].npos)
				{
					PROCEDURE[j].replace(PROCEDURE[j].find("UNITSOFF"), 8, "");
					for (size_t k = j; k<PROCEDURE.size(); k++)
					{
						if ((PROCEDURE[k].find(":")==PROCEDURE[k].npos) && (PROCEDURE[k].find("PROCEDURE")==PROCEDURE[i].npos))
						{
							fileInfo.push_back(PROCEDURE[k]);
						}
					}
				}
			}
		}
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
	endvars = s.find(")");
	//std::cout << "beg " << endf << " end " << endvars << std::endl;
	var = s.substr(endf+1, endvars-endf-1);
	//std::cout << var << std::endl;
	buf = 0;
	// writs all vars out
	while (var.find(",", buf)!=var.npos)
	{
		buf1 = var.find(",", buf);
		vars.push_back(var.substr(buf, buf1-buf));
		buf = buf1+1;
	}
	vars.push_back(var.substr(buf, var.npos));

	out += fname + "(";

	for (size_t i=0; i<vars.size()-1; i++)
	{
		out += "double " + vars[i] + ", ";
	}
	out += "double " + vars[vars.size()-1] + ")";
    return out;
}


size_t Converter::count_beg(string beg)
{
	size_t begi = 0;
	if (beg.find("\t")!=beg.npos)
	{
		begi = 1;
	}

	while (beg.find(" ", begi)!=beg.npos)
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


size_t Converter::number_(size_t pos, string s)
{
	size_t counter = 0;
	size_t begin;
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
	//std::cout << "pos_letter does anything" << std::endl;
	for ( char i( 97 ); i < 122; ++i )
	{
		//97-122
		if (s.find(i)!=s.npos)
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
	for ( char i( 97 ); i < 122; ++i )
	{
		//97-122
		if (s.find(i)!=s.npos)
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
	for (size_t i=0; i<Pairs.size(); i++)
	{
		//std::cout << "first" << Pairs.size() << std::endl;
		if (Zeilen[Pairs[i].first].find(name)!=Zeilen[Pairs[i].first].npos)
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

			// schreibt zeilen raus
			for (int j=Pairs[i].first; j<Pairs[i].second; j++)
			{
				//std::cout << counter << std::endl;
				BlockName.push_back(Zeilen[j]);
			}

		}
	}
	std::cout << "before block" << std::endl;
	std::cout << BlockName.size() << std::endl;
	return BlockName;
}


void Converter::WriteStart(string filename, std::vector<pair<int, int> > Pairs, std::vector<string> Zeilen)
{
	  // defines filenames as const char*
	  string fnamecpp = filename + ".cpp";
	  string fnameh = filename + ".h";
	  const char* filenamecpp = fnamecpp.c_str();
	  const char* filenameh = fnameh.c_str();



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
	  mycppfile << "namespace ug { \n \n \n";

	  std::vector<string> mod_funcs_names;
	  // all function have to be read in first so later need more attention more

	  std::vector<string> FUNCTION = GetBlock(Pairs, Zeilen, "FUNCTION");

	  string func, funchead;
	  if (FUNCTION.size() > 0)
	  {
		  func = build_func_head(FUNCTION[0]);
	  	  funchead = func_head(FUNCTION[0]);
	  }
	  size_t funcS, funcE;
	  mycppfile << func +" \n";
	  mycppfile << "{ \n";
	  for (size_t i = 1; i<FUNCTION.size(); i++)
	  {
		  funcS = FUNCTION[i].find(funchead);
		  if (funcS!=FUNCTION[i].npos)
		  {
			 funcE = FUNCTION[i].find(" ", funcS);
			 FUNCTION[i].replace(funcS, funcE-funcS + 2, "return ");
		  }
		  mycppfile << FUNCTION[i] + "\n";
	  }

	  mycppfile << "\n \n";
	  mycppfile.close();

	  ofstream myhfile;
	  myhfile.open (filenameh);
	  myhfile << "#ifndef " + filename + "_H_\n";
	  myhfile << "#define " + filename + "_H_\n";
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
	  for (size_t i = 0; i<Global_vars.size(); i++)
	  {
		  myhfile << "double " + Global_vars[i] + "; \n";
	  }
	  myhfile << "\n \n \n" ;



	  // getting out of Params hard coded Values
	  std::vector<string> PARAMETER = GetBlock(Pairs, Zeilen, "PARAMETER");
	  size_t end;
	  string var;


	  for (size_t i=0; i<PARAMETER.size(); i++)
	  {
		  if (PARAMETER[i].find("=")!=PARAMETER[i].npos)
		  {
			   end = PARAMETER[i].find("(");
			   var = PARAMETER[i].substr(0, end);
			   myhfile << "const static " + var + "; \n";
		  }
	  }
	  myhfile << "std::vector<number> m_diff; \n";
	  myhfile << "\n";


	  size_t Ion, IonEnd;
	  string IonS;
	  std::vector<string> ListIons;

	  for (size_t i=0; i<NEURON.size(); i++)
	  {
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

		  }
	  }


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
	  size_t state, stateend;
	  size_t varSt = 0;
	  string addState;




	  std::cout << STATE.size() << std::endl;
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
			  addState = STATE[i].substr(state, stateend-state);
			  //writting in hfile
			  myhfile << "ADouble " + addState + "Gate; \n";
			  myhfile << "Grid::AttachmentAccessor<Vertex, ADouble> aa" + addState + "Gate; \n";

/////////////////////////////////////////////////////////////////////////
//writting of attachment inits for cpp file
/////////////////////////////////////////////////////////////////////////

			  mycppfile << "if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->" + addState + "Gate)) \n";
			  mycppfile << "UG_THROW(\"Attachment necessary (" + addState + "Gate) for Hodgkin and Huxley channel dynamics \"\n";
			  mycppfile << "\"could not be made, since it already exists.\"); \n";
			  mycppfile << "spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->" + addState + "Gate); \n";
			  mycppfile << "this->aa" + addState + "Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->" + addState + "Gate); \n \n";

			  state = stateend+1;
		  }
		  varSt = 0;
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

	  size_t begin;
	  size_t Stats_beg = 1;
	  if (INITIAL.size() > 1)
	  {
		  begin = count_beg(INITIAL[1]);
	  }
	  //std::cout << "beg: " << begin << std::endl;

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
					  mycppfile << "double " + Proc_funcs[j][k] + " = aa"+ Proc_funcs[j][k] + "Gate[*iter]; \n \n";
				  }

				  // write locals vals as double
				  locals = get_local_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);

				  for (size_t k=0; k<locals.size(); k++)
				  {
					  mycppfile << "double " + locals[k] + "; \n";
				  }

				  Proc_vals = writer_proc_block(Pairs, Zeilen, Proc_funcs[j][0]);
				  std::cout << Proc_vals.size() << std::endl;


				  for (size_t k=0; k<Proc_vals.size(); k++)
				  {
					  if ((Proc_vals[k]!="") && (Proc_vals[k]!="}"))
						  mycppfile << Proc_vals[k] + "; \n";

				  }
			  }
		  }

		  if (i > Stats_beg)
		  {
			  size_t gleich, vorgleich;
			  if ((INITIAL[i]!="") && (INITIAL[i]!="}"))
			  {
				  gleich = INITIAL[i].find("=");
			      vorgleich = count_beg(INITIAL[i])-4;
				  mycppfile << "aa" + INITIAL[i].substr(vorgleich , gleich-vorgleich-1) + "Gate[*iter] " + INITIAL[i].substr(gleich) + "; \n";
			  }
		  }



	  }

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
	  std::vector<vector <string> > Der_funcs;

	  Der_funcs = write_derivative_block(Pairs, Zeilen);

	  if (DERIVATIVE.size() > 1)
	  {
		  begin = count_beg(DERIVATIVE[1]);
	  }
	  //std::cout << "beg: " << begin << std::endl;
	  // writting ions

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


	  // States ever needed for deriv
	  for (size_t i=0; i<STATE.size(); i++)
	  {
		  state = pos_letter(STATE[i]);
		  if (state!=1000)
			  if (STATE[i].find("(")==STATE[i].npos)
				  varSt = number_(state, STATE[i]);
		  //std::cout << "after doppelt if" << std::endl;
		  for (size_t j=0; j<varSt; j++)
		  {
			  stateend = STATE[i].find(" ", state);
			  addState = STATE[i].substr(state, stateend-state);
			  mycppfile << "double " + addState + " = aa" + addState + "Gate[*iter]; \n";

			  state = stateend+1;
		  }
		  varSt = 0;
	  }

	  mycppfile << "\n \n \n";

	  // writting needed for DERIVATIVE
	  for (size_t i=0; i<DERIVATIVE.size(); i++)
	  {

		  for (size_t j=0; j<Der_funcs.size(); j++)
		  {
			  begin = count_beg(DERIVATIVE[i]);
			  if (Der_funcs[j][0] == DERIVATIVE[i].substr(begin, DERIVATIVE[i].find("(")-begin))
			  {
				  Stats_beg = i;
				  // write needed values
				  for (size_t k=1; k<Der_funcs[j].size(); k++)
				  {
					  mycppfile << "double " + Der_funcs[j][k] + " = aa"+ Der_funcs[j][k] + "Gate[*iter]; \n \n";
				  }
			  }
		  }
		  // Writting derivfuncs
		  if ((i>1) && (i < DERIVATIVE.size()-1))
		  {
			DERIVATIVE[i].replace(DERIVATIVE[i].find("\'"), 1, "");
			DERIVATIVE[i].replace(DERIVATIVE[i].find("="), 1, "+=");
			mycppfile << DERIVATIVE[i] + "*dt; \n";
		  }


	  }
	  mycppfile << "\n \n \n";


	  // States ever need new values
		  for (size_t i=0; i<STATE.size(); i++)
		  {
			  state = pos_letter(STATE[i]);
			  if (state!=1000)
				  if (STATE[i].find("(")==STATE[i].npos)
					  varSt = number_(state, STATE[i]);
			  //std::cout << "after doppelt if" << std::endl;
			  for (size_t j=0; j<varSt; j++)
			  {
				  stateend = STATE[i].find(" ", state);
				  addState = STATE[i].substr(state, stateend-state);
				  mycppfile << "aa" + addState + "Gate[*iter] = " + addState + "; \n";

				  state = stateend+1;
			  }
			  varSt = 0;
		  }

	  mycppfile << " \n \n \n";
	  mycppfile << "} \n";
	  mycppfile << "} \n";
	  mycppfile << "} \n";
	  mycppfile << " \n \n \n";

	  // write head of ionic flux
	  mycppfile << "template<typename TDomain, typename TAlgebra> \n";
	  mycppfile << "void " + filename + "<TDomain, TAlgebra>::ionic_current(Vertex* ver, std::vector<number>& outCurrentValues) \n";
	  mycppfile << "{ \n \n";
	  // all gates

	  // States ever need new values
	  for (size_t i=0; i<STATE.size(); i++)
	  {
		  state = pos_letter(STATE[i]);
		  if (state!=1000)
			  if (STATE[i].find("(")==STATE[i].npos)
				  varSt = number_(state, STATE[i]);
		  //std::cout << "after doppelt if" << std::endl;
		  for (size_t j=0; j<varSt; j++)
		  {
			  stateend = STATE[i].find(" ", state);
			  addState = STATE[i].substr(state, stateend-state);
			  mycppfile << "double " + addState + " = aa" + addState + "[ver]; \n";
			  state = stateend+1;
		  }
		  varSt = 0;
	  }
	  // all ions
	  for (size_t i=0; i<ListIons.size(); i++)
	  {
	 	  mycppfile << "double " + ListIons[i] + " = aa" + ListIons[i] + "[ver];  \n";
	  }
	  // v
	  mycppfile << "double v = m_aav[ver]; \n";
	  mycppfile << " \n";
	  mycppfile << " \n";

	  //writes Nernst Eqs
	  std::vector<string> eqs = equali(Pairs, Zeilen);
	  for (size_t i=0; i<eqs.size(); i++)
	  {
		  mycppfile << eqs[i] + " \n";
	  }
	  mycppfile << " \n \n";

	  bool in = true;
	  size_t beg;
	  std::vector<string> fluxes;
	  std::vector<string> B_vars;
	  std::vector<string> BREAKPOINT = GetBlock(Pairs, Zeilen, "BREAKPOINT");

	  for (size_t i=1; i<BREAKPOINT.size()-1; i++)
	  {
		  in = true;
		  if (BREAKPOINT[i].find("SOLVE")!=BREAKPOINT[i].npos)
			  in = false;

		  if (in == true)
		  {
			  beg = pos_letter(BREAKPOINT[i]);
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

	  vector<string> outs;
	  // first every time v
	  outs.push_back("");

	  mycppfile << "\n \n \n";
	  for (size_t i=0; i<fluxes.size(); i++)
	  {
		  std::cout << "fluxes[i]: " << fluxes[i] << std::endl;
		  if (fluxes[i]!="")
		  {
			  if (fluxes[i].find("i")==0)
				  outs[0] += " + " + fluxes[i].substr(fluxes[i].find("=")+1, fluxes[i].npos -fluxes[i].find("=")+1);
		  }
	  }
	  for (size_t i=0; i<fluxes.size(); i++)
	  {
		  for (size_t j=0; j<ListIons.size(); j++)
		  {
			  std::cout << "List Ions" << std::endl;
			  std::cout << "i"+ListIons[j] << std::endl;
			  if (fluxes[i].find("i" + ListIons[j])!=fluxes[i].npos)
				  outs.push_back(fluxes[i].substr(fluxes[i].find("=")+1, fluxes[i].npos -fluxes[i].find("=")+1));
		  }
	  }
	  for (size_t i=0; i<outs.size(); i++)
	  {
		  mycppfile << "outCurrentValues.push_back(" + outs[i] + "/m_F ) \n";
	  }
	  mycppfile << "  \n";
	  mycppfile << "}  \n";



	  mycppfile << "  \n";
	  mycppfile << "  \n";
	  mycppfile.close();


	  myhfile << "ADouble v; \n";
	  myhfile << "Grid::AttachmentAccessor<Vertex, ADouble> aav; \n \n";
	  myhfile << "//nernst const values \n";

	  myhfile << "number m_R \n";
	  myhfile << "number m_T \n";
	  myhfile << "number m_F \n \n";

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

