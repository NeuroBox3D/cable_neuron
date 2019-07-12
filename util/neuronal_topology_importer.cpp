/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Original author: Stephan Grein
 * Creation date: 2014-10-14
 * Restructuring and augmentation: Markus Breit (2019-04)
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "neuronal_topology_importer.h"

#include "cn_config.h"  // for project-specific definitions (NETI_ALL_FEATURES)
#include "../synapse_handling/synapses/onset_pre_synapse.h"
#include "../synapse_handling/synapses/threshold_pre_synapse.h"
#include "../synapse_handling/synapses/alpha_post_synapse.h"
#include "../synapse_handling/synapses/exp2_post_synapse.h"

// ug includes
#include "common/debug_id.h"  // DebugID
#include "common/log.h"  // UG_LOGN
#include "common/math/misc/eigenvalues.h"  // CalculateEigenvalues
#include "common/parser/rapidxml/rapidxml.hpp"
#include "common/util/file_util.h"  // FindFileInStandardPaths
#include "common/util/string_util.h"  // TokenizeString etc.
#include "lib_grid/algorithms/element_side_util.h"  // GetOpposingSide
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"  // MergeVertices
#include "lib_grid/algorithms/subset_util.h"  // AssignSubsetColors, EraseEmptySubsets
#include "lib_grid/file_io/file_io_ugx.h"  // GridWriterUGX

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <algorithm>  // std::lower_bound
#include <fstream>
#include <iterator>  // std::distance
#include <limits>
#if !defined(UG_FOR_VRL) && defined(NETI_ALL_FEATURES)
#include <regex>
#endif
#include <sstream>
#include <utility>  // std::pair


namespace ug {
namespace cable_neuron {
namespace neuronal_topology_importer {
	
const number NeuronalTopologyImporter::REMOVE_DOUBLES_THRESHOLD = 1e-8;

NeuronalTopologyImporter::NeuronalTopologyImporter()
: m_scaleLength(1e-6),  // um -> m
  m_scaleTime(1e-3),  // ms -> s
  m_scalePotential(1e-3),  // mV -> V
  m_scaleConductance(1e-6),  // uS -> S
  m_grid(SPNULL),
  m_sh(SPNULL),
  m_networkID(0),
  m_numCellTypes(5),
  m_bWithCellType(false)
{}


NeuronalTopologyImporter::~NeuronalTopologyImporter()
{}



void NeuronalTopologyImporter::set_scaling_factors(number length, number time, number potential, number conductance)
{
	m_scaleLength = length;
	m_scaleTime = time;
	m_scalePotential = potential;
	m_scaleConductance = conductance;
}

void NeuronalTopologyImporter::set_length_scaling(number lengthScale)
{
	m_scaleLength = lengthScale;
}

void NeuronalTopologyImporter::set_time_scaling(number timeScale)
{
	m_scaleTime = timeScale;
}

void NeuronalTopologyImporter::set_potential_scaling(number potentialScale)
{
	m_scalePotential = potentialScale;
}

void NeuronalTopologyImporter::set_conductance_scaling(number conductanceScale)
{
	m_scaleConductance = conductanceScale;
}


bool NeuronalTopologyImporter::import_and_generate_grid(const std::string& file)
{
	return import_and_generate_grid(FilenameAndPathWithoutExtension(file), GetFilenameExtension(file));
}



void NeuronalTopologyImporter::add_joining_criterion(const char* str)
{
	this->m_joiningCriteria.push_back(std::string(str));
}


void NeuronalTopologyImporter::rm_joining_criterion(const char* str)
{
	std::vector<std::string>::iterator it =
		std::find(m_joiningCriteria.begin(), m_joiningCriteria.end(), std::string(str));
	if (it != m_joiningCriteria.end()) m_joiningCriteria.erase(it);
}




bool NeuronalTopologyImporter::import_and_generate_grid(const std::string& file, const std::string& type)
{
	if (type == "txt")
		m_filename = FindFileInStandardPaths((file + "_secs." + type).c_str());
	else
		m_filename = FindFileInStandardPaths((file + "." + type).c_str());
	UG_COND_THROW(m_filename == "", "File '" << file + "." + type << "' not found in standard paths.");
	m_filename = FilenameAndPathWithoutExtension(m_filename);
	if (type == "txt")
		m_filename = m_filename.substr(0, m_filename.size() - 5);
	m_extension = type;

	// import file
	if (!import_geometry(file))
		return false;

	// create grid
	if (!generate_grid())
		return false;

	// handle subsets
	subset_processing();

	// save grid and indicate success
	return SaveGridToFile(*m_grid, *m_sh, (m_filename + ".ugx").c_str());
}






bool NeuronalTopologyImporter::import_geometry(const std::string& file)
{
	if (m_extension == "hoc")
		return import_geometry_from_HOC();

	if (m_extension == "swc")
		return import_geometry_from_SWC();

	if (m_extension == "ngx")
		return import_geometry_from_NGX();

	if (m_extension == "txt")
		return import_geometry_from_TXT();

	if (m_extension == "xml")
		return import_geometry_from_NML();

	UG_DLOGN(NETIGeometry, 1, "UNKNOWN filetype '" << m_extension << "'.");
	return false;
}



bool NeuronalTopologyImporter::import_geometry_from_SWC()
{
	// 1: open file
	std::ifstream inFile((m_filename + "." + m_extension).c_str());
	std::string line;

	if (!inFile)
	{
		UG_DLOGN(NETIGeometry, 1, "File **" << m_filename << "." << m_extension << " ** does not exist.");
		return false;
	}

	// 2: read line by line and split
	m_vPtInfo.clear();
	std::map<size_t, size_t> fileInd2ourInd;
	size_t lineCnt = 0;
	while (std::getline(inFile, line))
	{
		++lineCnt;

		// 2.1: ignore lines starting with #
		if (line.size() && line[0] == '#')
			continue;

		// 3: split the line into tokens
		std::istringstream buf(line);
		std::istream_iterator<std::string> beg(buf), end;
		std::vector<std::string> strs(beg, end);

		// 4: assert if tokens are the right amount
		UG_ASSERT(strs.size() == 7, "SWC not in standardized format, "
			"i. e., columns do not match the format specification.");

		// 5: populate m_compartments
		m_vPtInfo.resize(m_vPtInfo.size() + 1);
		PointInfo& swc = m_vPtInfo.back();

		// 5.1: get index and save in tmp index conversion map
		const size_t index = boost::lexical_cast<size_t>(strs[0]);
		std::map<size_t, size_t>::iterator secIt = fileInd2ourInd.lower_bound(index);
		UG_COND_THROW(secIt != fileInd2ourInd.end() && !fileInd2ourInd.key_comp()(index, secIt->first),
			"Double index " << index << " on line " << lineCnt << ".");
		secIt = fileInd2ourInd.insert(secIt, std::map<size_t, size_t>::value_type(index, m_vPtInfo.size() - 1));

		// 5.2 get type
		swc.index = boost::lexical_cast<size_t>(strs[1]);

		// 5.3: set coordinates
		swc.coordinates[0] = boost::lexical_cast<number>(strs[2]) * m_scaleLength;
		swc.coordinates[1] = boost::lexical_cast<number>(strs[3]) * m_scaleLength;
		swc.coordinates[2] = boost::lexical_cast<number>(strs[4]) * m_scaleLength;

		// 5.3: set diameter (SWC format uses radius -> 2 * radius = diameter)
		swc.diameter = boost::lexical_cast<number>(strs[5]) * m_scaleLength * 2;

		// 5.4: set connectivity
		const int connectsToInd = boost::lexical_cast<int>(strs[6]);
		if (connectsToInd < 0)
			swc.connectsTo = (size_t) -1;
		else
		{
			if ((secIt = fileInd2ourInd.find((size_t) connectsToInd)) == fileInd2ourInd.end())
				UG_THROW("Connection to unknown index " << connectsToInd << " on line " << lineCnt << ".");
			swc.connectsTo = secIt->second;
		}
	}
	return true;
}


#if !defined(UG_FOR_VRL) && defined(NETI_ALL_FEATURES)
/** \brief slurps whole file into string
 * \param[in] filename the file
 * \param[in] string the whole file in one string
 * \return \c success or not
**/
static bool slurp(const std::string& filename, std::string& string)
{
	std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
	if (in)
	{
		std::string contents;
		in.seekg(0, std::ios::end);
		contents.resize(in.tellg());
		in.seekg(0, std::ios::beg);
		in.read(&contents[0], contents.size());
		in.close();
		string = contents;
		return true;
	}
	return false;
}
#endif



bool NeuronalTopologyImporter::import_geometry_from_HOC()
{
#if !defined(UG_FOR_VRL) && defined(NETI_ALL_FEATURES)
	// open file and read whole file into string
	std::string text;
	if (!slurp(m_filename + "." + m_extension, text))
	{
		UG_LOGN("File '" << m_filename << "." << m_extension << "' could not be read.");
		return false;
	}

	// read neurites from file
	try {extract_points_and_sections_from_hoc(text);}
	UG_CATCH_THROW("Points and section extraction from hoc file failed.");

	// read connection information from file
	try {extract_connections_from_hoc(text);}
	UG_CATCH_THROW("Connection extraction from hoc file failed.");

	// read synapses from file
	try {extract_synapses_from_hoc(text);}
	UG_CATCH_THROW("Synapse extraction from hoc file failed.");

	return true;

#endif

	UG_LOGN("Either your compiler does not support C++11 or this feature "
		"has been explicitly switched off at compile-time.");

	return false;
}



bool NeuronalTopologyImporter::import_geometry_from_NGX()
{
	// open file
	std::ifstream geometryFile((m_filename + "." + m_extension).c_str());
	UG_COND_THROW(!geometryFile, "File could not be opened: " << m_filename << "." << m_extension);

	// read xml structure into tree
	boost::property_tree::ptree tree;
	read_xml(geometryFile, tree);

	// read neurites from file
	try {extract_points_and_sections_from_ngx(tree);}
	UG_CATCH_THROW("Extraction of points and sections from NGX file failed.");

	// read connection information from file
	try {extract_connections_from_ngx(tree);}
	UG_CATCH_THROW("Extraction of connections from NGX file failed.");

	// read synapses from file
	try {extract_synapses_from_ngx(tree);}
	UG_CATCH_THROW("Extraction of synapses from NGX file failed.");

	return true;
}



bool NeuronalTopologyImporter::import_geometry_from_TXT()
{
	const std::string fileNameID = m_filename + "_identifier." + m_extension;
	const std::string fileNameSecs = m_filename + "_secs." + m_extension;
	const std::string fileNameConn = m_filename + "_connex." + m_extension;
	const std::string fileNameSyn = m_filename + "_synapses." + m_extension;

	// read identifier file
	try {read_identifier_file_txt(fileNameID);}
	UG_CATCH_THROW("Failed reading identifier file.");

	m_numCellTypes = m_networkID ? 6 : 5;

	// read sections file
	//clock_t begin = clock();
	try {extract_points_and_sections_from_txt(fileNameSecs);}
	UG_CATCH_THROW("Failed extracting points and sections from txt file.");
	//clock_t end = clock();
	//UG_LOGN("Sections read from file (#" << m_mSecInfo.size() << "). Elapsed time: "
	//	<< double (end - begin) / CLOCKS_PER_SEC << "s.");

	// read connections file
	//begin = clock();
	try {extract_connections_from_txt(fileNameConn);}
	UG_CATCH_THROW("Failed extracting connections from txt file.");
	//end = clock();
	//UG_LOGN("Connections read from file. Time elapsed: "
	//		<< double (end - begin) / CLOCKS_PER_SEC << "s.");

	// read synapses file
	//begin = clock();
	try {extract_synapses_from_txt(fileNameSyn);}
	UG_CATCH_THROW("Failed extracting synapses from txt file.");
	//end = clock();
	//UG_LOGN("Synapses read from file. Time elapsed: "
	//	<< double (end - begin) / CLOCKS_PER_SEC << "s.");

	return true;
}


static void read_soma_points_nml
(
	std::vector<vector3>& vSomaPtsOut,
	const rapidxml::xml_document<char>& doc
)
{
	rapidxml::xml_node<>* mbf = doc.first_node("mbf");
	UG_COND_THROW(!mbf, "No 'mbf' tag contained in file.")
	rapidxml::xml_node<>* curNode = mbf->first_node("contour");
	UG_COND_THROW(!curNode, "No soma contained in file.")

	// read soma points
	vSomaPtsOut.clear();
	rapidxml::xml_node<>* curPtNode = curNode->first_node("point");
	UG_COND_THROW(!curNode, "No points for soma.")
	while (curPtNode)
	{
		vSomaPtsOut.resize(vSomaPtsOut.size()+1);
		vector3& coords = vSomaPtsOut.back();
		rapidxml::xml_attribute<>* attrib = curPtNode->first_attribute("x");
		UG_COND_THROW(!attrib, "No x coordinate given for point.");
		coords[0] = strtod(attrib->value(), NULL);
		attrib = curPtNode->first_attribute("y");
		UG_COND_THROW(!attrib, "No y coordinate given for point.");
		coords[1] = strtod(attrib->value(), NULL);
		attrib = curPtNode->first_attribute("z");
		UG_COND_THROW(!attrib, "No z coordinate given for point.");
		coords[2] = strtod(attrib->value(), NULL);

		curPtNode = curPtNode->next_sibling("point");
	}
}



static void soma_contour_to_segments_nml
(
	std::vector<PointInfo>& vPtInfoOut,
	const std::vector<vector3>& contourPts,
	number scaleLength
)
{
	// Conversion is done by identifying the principal components
	// and then connecting the outermost points w.r.t. the biggest component
	// with 10 cylinder segments with radius and center chosen such that
	// they approximate the given contour as closely as possible.
	// This will only work properly for convex contours.
	size_t sz = contourPts.size();
	UG_COND_THROW(!sz, "Cannot convert soma contour to segments if the contour has no points.");

	// center
	vector3 center(0.0);
	for (size_t i = 0; i < sz; ++i)
		VecAdd(center, center, contourPts[i]);
	VecScale(center, center, 1.0/sz);

	// "covariance" matrix
	vector3 c0(0.0); // matrix columns
	vector3 c1(0.0);
	vector3 c2(0.0);
	for (size_t i = 0; i < sz; ++i)
	{
		vector3 help;
		VecSubtract(help, contourPts[i], center);
		VecScaleAdd(c0, 1.0, c0, help[0], help);
		VecScaleAdd(c1, 1.0, c1, help[1], help);
		VecScaleAdd(c2, 1.0, c2, help[2], help);
	}
	matrix33 mat;
	mat(0,0) = c0[0]; mat(0,1) = c1[0]; mat(0,2) = c2[0];
	mat(1,0) = c0[1]; mat(1,1) = c1[1]; mat(1,2) = c2[1];
	mat(2,0) = c0[2]; mat(2,1) = c1[2]; mat(2,2) = c2[2];

	// principal components
	number l0, l1, l2;
	vector3 e0, e1, e2;
	bool success = CalculateEigenvalues(mat, l0, l1, l2, e0, e1, e2);
	UG_COND_THROW(!success, "Failed to determine principal components of soma.");
	VecNormalize(e2, e2);
	VecNormalize(e1, e1);

	// min and max points
	number min = std::numeric_limits<number>::max();
	number max = std::numeric_limits<number>::min();
	size_t minInd = 0;
	size_t maxInd = 0;
	for (size_t i = 0; i < sz; ++i)
	{
		vector3 help;
		VecSubtract(help, contourPts[i], center);
		number p = help * e2;
		if (p < min)
		{
			minInd = i;
			min = p;
		}
		if (p > max)
		{
			maxInd = i;
			max = p;
		}
	}
	vector3 minPt = contourPts[minInd];
	vector3 maxPt = contourPts[maxInd];

	vector3 line;
	VecSubtract(line, maxPt, minPt);

	// create segments
	vPtInfoOut.resize(10);
	size_t curDown = minInd;
	size_t curUp = minInd;
	for (size_t k = 0; k < 10; ++k)
	{
		// segment pos along line
		vector3& segPos = vPtInfoOut[k].coordinates;
		VecScaleAdd(segPos, 1.0, minPt, 0.1*k+0.05, line);

		// advance up / down until above / below current segment
		vector3 help;
		VecSubtract(help, contourPts[curDown], segPos);
		while (VecProd(help, e2) < 0)
		{
			curDown = (curDown + sz - 1) % sz;  // prevent underflow
			VecSubtract(help, contourPts[curDown], segPos);
		}
		VecSubtract(help, contourPts[curUp], segPos);
		while (VecProd(help, e2) < 0)
		{
			curUp = (curUp + 1) % sz;
			VecSubtract(help, contourPts[curUp], segPos);
		}

		// position segment in the middle w.r.t. e1
		VecScaleAdd(help, 0.5, contourPts[curDown], 0.5, contourPts[(curDown+1)%sz], -1.0, segPos);
		number distDown = VecProd(help, e1);
		VecScaleAdd(help, 0.5, contourPts[curUp], 0.5, contourPts[(curUp-1)%sz], -1.0, segPos);
		number distUp = VecProd(help, e1);
		VecScaleAdd(segPos, 1.0, segPos, (distUp+distDown)/2.0, e1);
		VecScale(segPos, segPos, scaleLength);

		vPtInfoOut[k].diameter = scaleLength * (fabs(distUp) + fabs(distDown));

		// point type is soma
		vPtInfoOut[k].index = 0;

		// connection
		vPtInfoOut[k].connectsTo = k - 1;
		vPtInfoOut[k].isDuplicate = false;
	}
}



static void createNeuroMLBranch
(
	std::vector<PointInfo>& vPtInfoOut,
	int si,
	size_t connectsTo,
	rapidxml::xml_node<char>* baseNode,
	number scaleLength
)
{
	// recursively create
	rapidxml::xml_node<char>* ptOrBranch = ((rapidxml::xml_node<>*)baseNode)->first_node();
	while (ptOrBranch)
	{
		const char* name = ptOrBranch->name();

		// regular point
		if (!strcmp(name, "point"))
		{
			const size_t nextPtIndex = vPtInfoOut.size();
			vPtInfoOut.resize(nextPtIndex+1);
			PointInfo& ptInfo = vPtInfoOut.back();

			// point type
			ptInfo.index = si;

			// coords
			vector3& coords = ptInfo.coordinates;
			rapidxml::xml_attribute<>* attrib = ptOrBranch->first_attribute("x");
			UG_COND_THROW(!attrib, "No x coordinate given for point.");
			coords[0] = scaleLength * strtod(attrib->value(), NULL);
			attrib = ptOrBranch->first_attribute("y");
			UG_COND_THROW(!attrib, "No y coordinate given for point.");
			coords[1] = scaleLength * strtod(attrib->value(), NULL);
			attrib = ptOrBranch->first_attribute("z");
			UG_COND_THROW(!attrib, "No z coordinate given for point.");
			coords[2] = scaleLength * strtod(attrib->value(), NULL);

			// diameter
			attrib = ptOrBranch->first_attribute("d");
			UG_COND_THROW(!attrib, "No diameter info for point.");
			ptInfo.diameter = scaleLength * strtod(attrib->value(), NULL);

			// connection
			ptInfo.connectsTo = connectsTo;
			connectsTo = nextPtIndex;
			ptInfo.isDuplicate = false;
		}

		// new branch
		else if (!strcmp(name, "branch"))
		{
			createNeuroMLBranch(vPtInfoOut, si, connectsTo, ptOrBranch, scaleLength);
		}

		ptOrBranch = ptOrBranch->next_sibling();
	}
}



static void add_neurite_points_nml
(
	std::vector<PointInfo>& vPtInfoOut,
	const rapidxml::xml_document<char>& doc,
	number scaleLength
)
{
	rapidxml::xml_node<char>* mbf = doc.first_node("mbf");
	UG_COND_THROW(!mbf, "No 'mbf' tag contained in file.")

	rapidxml::xml_node<char>* curTree = mbf->first_node("tree");

	const size_t nSomaSeg = vPtInfoOut.size();
	while (curTree)
	{
		// get type of neurite
		int si;
		rapidxml::xml_attribute<>* attrib = curTree->first_attribute("type");
		UG_COND_THROW(!attrib, "No type attribute given for tree.");
		const std::string type(attrib->value());

		if (type == std::string("Axon"))
			si = 1;
		else if (type == std::string("Dendrite"))
			si = 2;
		else if (type == std::string("Apical"))
			si = 3;
		else
			UG_THROW("Unknown neurite type '" << type << "' in tree.");

		// find connecting soma vertex (nearest one)
		rapidxml::xml_node<char>* ptOrBranch = curTree->first_node();
		UG_COND_THROW(!ptOrBranch, "Empty tree.");
		const char* name = ptOrBranch->name();
		UG_COND_THROW(strcmp(name, "point"),
			"Trees cannot start with anything other than a point.");

		// get first point coords
		vector3 coords;
		attrib = ptOrBranch->first_attribute("x");
		UG_COND_THROW(!attrib, "No x coordinate given for point.");
		coords[0] = scaleLength * strtod(attrib->value(), NULL);
		attrib = ptOrBranch->first_attribute("y");
		UG_COND_THROW(!attrib, "No y coordinate given for point.");
		coords[1] = scaleLength * strtod(attrib->value(), NULL);
		attrib = ptOrBranch->first_attribute("z");
		UG_COND_THROW(!attrib, "No z coordinate given for point.");
		coords[2] = scaleLength * strtod(attrib->value(), NULL);

		size_t connPt = 0;
		number minDist = std::numeric_limits<number>::max();
		for (size_t i = 0; i < nSomaSeg; ++i)
		{
			number dist = VecDistanceSq(coords, vPtInfoOut[i].coordinates);
			if (dist < minDist)
			{
				minDist = dist;
				connPt = i;
			}
		}
		UG_COND_THROW(minDist == std::numeric_limits<number>::max(),
			"No connecting soma vertex found for tree root vertex.");

		// create the whole tree
		createNeuroMLBranch(vPtInfoOut, si, connPt, curTree, scaleLength);

		// get the next tree
		curTree = curTree->next_sibling("tree");
	}
}



bool NeuronalTopologyImporter::import_geometry_from_NML()
{
	const std::string fileNameIn = m_filename + "." + m_extension;
	std::ifstream in(fileNameIn.c_str(), std::ios::binary);
	UG_COND_THROW(!in, "Could not open file '" << fileNameIn << "' for conversion.");

	rapidxml::xml_document<> doc;

	// get the length of the file
	std::streampos posStart = in.tellg();
	in.seekg(0, std::ios_base::end);
	std::streampos posEnd = in.tellg();
	std::streamsize size = posEnd - posStart;

	// go back to the start of the file
	in.seekg(posStart);

	// read the whole file en-block and terminate it with 0
	char* fileContent = doc.allocate_string(0, size + 1);
	in.read(fileContent, size);
	fileContent[size] = 0;
	in.close();

	// parse the xml-data
	doc.parse<0>(fileContent);

	// create soma segments
	try
	{
		std::vector<vector3> vSomaContourPts;
		read_soma_points_nml(vSomaContourPts, doc);

		// convert soma contour to segments
		m_vPtInfo.clear();
		soma_contour_to_segments_nml(m_vPtInfo, vSomaContourPts, m_scaleLength);
	}
	UG_CATCH_THROW("Failed reading soma points.");

	// create dendrites and axons
	try {add_neurite_points_nml(m_vPtInfo, doc, m_scaleLength);}
	UG_CATCH_THROW("Failed reading neurite points.");



	return true;
}







static void get_synapse_edge_and_local_coords
(
	Edge** e,
	number& locCoord,
	const std::string& secName,
	number secLoc,
	const std::vector<PointInfo>& vPts,
	const std::map<std::string, SectionInfo>& mSecs,
	const std::vector<Vertex*>& vVrts,
	MultiGrid& mg
)
{
	std::map<std::string, SectionInfo>::const_iterator itSec = mSecs.find(secName);
	UG_COND_THROW(itSec == mSecs.end(), "Section '" << secName << "' requested, but unknown.");
	const SectionInfo& secInfo = itSec->second;

	// find edge's vertex indices and local coordinate
	size_t vi0 = 0;
	size_t vi1 = 1;
	locCoord = -1.0;
	if (secLoc <= 0.0)
	{
		vi0 = secInfo.start;
		vi1 = vi0 + 1;
		locCoord = 0.0;
	}
	else if (secLoc >= 1.0)
	{
		vi1 = secInfo.end;
		vi0 = vi1 - 1;
		locCoord = 1.0;
	}
	else
	{
		const size_t nPts = secInfo.end - secInfo.start + 1;
		std::vector<number> vCumulLength(nPts, 0.0);
		for (size_t i = 1; i < nPts; ++i)
			vCumulLength[i] = vCumulLength[i-1] + VecDistance(vPts[secInfo.start+i-1].coordinates,
															  vPts[secInfo.start+i].coordinates);

		const number searchLength = secLoc * vCumulLength[nPts-1];
		std::vector<number>::iterator it =
			std::lower_bound(vCumulLength.begin(), vCumulLength.end(), searchLength);

		const size_t offset = std::distance(vCumulLength.begin(), it);
		vi1 = (size_t) secInfo.start + offset;
		vi0 = vi1 - 1;
		locCoord = (searchLength - vCumulLength[offset-1]) / (vCumulLength[offset] - vCumulLength[offset-1]);
	}

	// find edge
	*e = NULL;
	Vertex* vrt0 = vVrts[vi0];
	Vertex* vrt1 = vVrts[vi1];
	Grid::traits<Edge>::secure_container el;
	mg.associated_elements(el, vrt0);
	const size_t elSz = el.size();
	for (size_t i = 0; i < elSz; ++i)
	{
		if (GetOpposingSide(mg, el[i], vrt0) == vrt1)
		{
			*e = el[i];
			break;
		}
	}

	UG_COND_THROW(!*e, "Edge could not be located.");
}



void NeuronalTopologyImporter::read_identifier_file_txt(const std::string& fileName)
{
	// default values
	m_networkID = 0;
	m_bWithCellType = false;

	//UG_LOGN("Reading identifier file '" << fileName << "'.");
	std::ifstream infile(fileName.c_str());
	if (!infile.is_open())
	{
		// if no network or cell type is specified we rely on defaults which should work in any case
		UG_WARNING("Identifier file could not be opened! Using defaults.");
		return;
	}

	// read network type (neocortex or hippocampus)
	std::string networkName;
	infile >> networkName;
	if (networkName == "Neocortex")
		m_networkID = 0;
	else if (networkName == "Hippocampus")
		m_networkID = 1;
	else
		UG_THROW("Unknown network type '" << networkName << "'.");

	// read whether subsets names are to reflect cell type
	std::string line;
	infile >> line;
	m_bWithCellType = line == "true";

	infile.close();
}



void NeuronalTopologyImporter::extract_points_and_sections_from_hoc(const std::string& text)
{
#if !defined(UG_FOR_VRL) && defined(NETI_ALL_FEATURES)
	// todo: using regex here is quite slow and demands a c++11-compatible compiler; rethink!
	const std::regex pt3dadd_match_pattern("\\s*\\(([-+]?\\d+\\.?\\d*)\\s*,\\s*([-+]?\\d+\\.?\\d*)\\s*,\\s*([-+]?\\d+\\.?\\d*)\\s*,\\s*(\\d+\\.?\\d*)\\)");
	const std::regex section_match_pattern("(pt3dadd\\(.*\\))");
	const std::regex pattern("([^\\r\\n(){}]+)\\s+\\{(?:\\s*pt3dclear\\(\\)|\\s*)([^]*?)\\}");

	size_t nSecs = 0;
	std::string sectionName;
	std::vector<std::string> temp;
	m_vPtInfo.clear();
	m_mSecInfo.clear();
	const std::sregex_token_iterator end;
	for (std::sregex_token_iterator i(text.cbegin(), text.cend(), pattern, 0); i != end; ++i)
	{
		const std::string& section_match_text = *i;
		TokenizeString(section_match_text, temp, '{');
		sectionName = TrimString(temp[0]);

		// filter out some undesired matches
		if (sectionName == "forall" || (sectionName.size() >= 7 && sectionName.substr(0, 7) == "forsec "))
			continue;
		UG_COND_THROW(temp[0].empty(), "Empty section name!");

		// enter sectionName in sectionName map
		std::map<std::string, SectionInfo>::iterator secIt = m_mSecInfo.lower_bound(sectionName);
		UG_COND_THROW(secIt != m_mSecInfo.end() && !m_mSecInfo.key_comp()(sectionName, secIt->first),
			"Section '" << sectionName << "' already exists, but is about to be added again. "
			"This is not allowed.");

		secIt = m_mSecInfo.insert(secIt, std::map<std::string, SectionInfo>::
			value_type(sectionName, SectionInfo()));
		SectionInfo& secInfo = secIt->second;
		secInfo.index = nSecs;
		secInfo.start = m_vPtInfo.size();
		++nSecs;


		// find all 3d point information
		bool firstPt = true;
		for (std::sregex_token_iterator j(section_match_text.cbegin(), section_match_text.cend(), section_match_pattern, 0); j != end; ++j)
		{
			const std::string& pt3dadd_match_text = *j;

			vector3 pos(0.0);
			number diam = 0.0;
			int gotAll = 0;
			for (std::sregex_token_iterator k(pt3dadd_match_text.cbegin(), pt3dadd_match_text.cend(), pt3dadd_match_pattern, 1); k != end; ++k)
			{
				pos[0] = m_scaleLength * boost::lexical_cast<number>(*k);
				gotAll |= 1;
			}
			for (std::sregex_token_iterator k(pt3dadd_match_text.cbegin(), pt3dadd_match_text.cend(), pt3dadd_match_pattern, 2); k != end; ++k)
			{
				pos[1] = m_scaleLength * boost::lexical_cast<number>(*k);
				gotAll |= 2;
			}
			for (std::sregex_token_iterator k(pt3dadd_match_text.cbegin(), pt3dadd_match_text.cend(), pt3dadd_match_pattern, 3); k != end; ++k)
			{
				pos[2] = m_scaleLength * boost::lexical_cast<number>(*k);
				gotAll |= 4;
			}
			for (std::sregex_token_iterator k(pt3dadd_match_text.cbegin(), pt3dadd_match_text.cend(), pt3dadd_match_pattern, 4); k != end; ++k)
			{
				diam = m_scaleLength * boost::lexical_cast<number>(*k);
				gotAll |= 8;
			}

			UG_COND_THROW(gotAll != 15, "Could not get all spatial coordinates or diameter.");

			// enter points in points map
			const size_t nPts = m_vPtInfo.size();
			m_vPtInfo.resize(nPts+1);
			PointInfo& ptInfo = m_vPtInfo.back();
			ptInfo.coordinates = pos;
			ptInfo.diameter = diam;
			ptInfo.index = nSecs - 1;
			if (!firstPt)
				ptInfo.connectsTo = nPts - 1;
			else
			{
				ptInfo.connectsTo = (size_t) -1;
				firstPt = false;
			}
		}

		secInfo.end = m_vPtInfo.size() - 1;
	}

	return;
#endif

	UG_THROW("Either your compiler does not support C++11 or this feature "
		"has been explicitly switched off at compile-time.");
}



void NeuronalTopologyImporter::extract_points_and_sections_from_ngx(const boost::property_tree::ptree& tree)
{
	using boost::property_tree::ptree;

	const ptree& sectionsTree = tree.get_child("ngx.sections");

	m_vPtInfo.clear();
	m_mSecInfo.clear();
	size_t nSecs = 0;
	BOOST_FOREACH(const ptree::value_type& v, sectionsTree)
	{
		if (v.first == "soma" || v.first == "axon" || v.first == "dendrite")
		{
			const std::string& sectionName = v.second.get<std::string>("name");

			std::map<std::string, SectionInfo>::iterator secIt = m_mSecInfo.lower_bound(sectionName);
			UG_COND_THROW(secIt != m_mSecInfo.end() && !m_mSecInfo.key_comp()(sectionName, secIt->first),
				"Section '" << sectionName << "' already exists, but is about to be added again. "
				"This is not allowed.");

			secIt = m_mSecInfo.insert(secIt, std::map<std::string, SectionInfo>::
				value_type(sectionName, SectionInfo()));
			SectionInfo& secInfo = secIt->second;
			secInfo.index = nSecs;
			secInfo.start = m_vPtInfo.size();
			++nSecs;

			//f.id = v.second.get<int>("id");
			bool firstPt = true;
			const ptree& subtree = v.second;
			BOOST_FOREACH(const ptree::value_type& vv, subtree)
			{
				if (vv.first == "coordinates")
				{
					const ptree& subtree2 = vv.second;
					BOOST_FOREACH(const ptree::value_type& vvv, subtree2)
					{
						if (vvv.first == "pt3d")
						{
							const size_t nPts = m_vPtInfo.size();
							m_vPtInfo.resize(nPts + 1);
							PointInfo& ptInfo = m_vPtInfo.back();
							ptInfo.index = nSecs - 1;
							vector3& coords = ptInfo.coordinates;
							coords[0] = m_scaleLength * vvv.second.get<number>("x");
							coords[1] = m_scaleLength * vvv.second.get<number>("y");
							coords[2] = m_scaleLength * vvv.second.get<number>("z");
							ptInfo.diameter = m_scaleLength * vvv.second.get<number>("w");
							if (!firstPt)
								ptInfo.connectsTo = nPts - 1;
							else
							{
								ptInfo.connectsTo = (size_t) -1;
								firstPt = false;
							}
						}
					}
				}
			}

			secInfo.end = m_vPtInfo.size() - 1;
		}
	}
}


void NeuronalTopologyImporter::extract_points_and_sections_from_txt(const std::string& fileName)
{
	//UG_LOGN("Reading sections from file '" << fileName << "'.");
	std::ifstream secFile(fileName.c_str());
	UG_COND_THROW(!secFile.is_open(), "Sections file could not be opened!");

	std::string sectionName;
	size_t nItems;
	size_t sectionType;
	size_t cellType;
	size_t nSecs = 0;
	while (true)
	{
		if (!(secFile >> sectionName))
			break; // regular end of file

		UG_COND_THROW(!(secFile >> nItems), "Could not read number of items.");
		UG_COND_THROW(!(secFile >> sectionType), "Could not read section type.");
		if (m_bWithCellType)
		{
			UG_COND_THROW(!(secFile >> cellType), "Could not read cell type.");
			sectionType = (sectionType * m_numCellTypes) + cellType - 1;
		}
		// check that section name does not exist already
		std::map<std::string, SectionInfo>::iterator secIt = m_mSecInfo.lower_bound(sectionName);
		UG_COND_THROW(secIt != m_mSecInfo.end() && !m_mSecInfo.key_comp()(sectionName, secIt->first),
			"Section '" << sectionName << "' already exists, but is about to be added again. "
			"This is not allowed.");

		// cerate new entry in section map
		secIt = m_mSecInfo.insert(secIt, std::map<std::string, SectionInfo>::
			value_type(sectionName, SectionInfo()));
		SectionInfo& secInfo = secIt->second;
		secInfo.index = nSecs;
		secInfo.start = m_vPtInfo.size();
		++nSecs;

		// read points
		bool firstPt = true;
		for (size_t i = 0; i < nItems; ++i)
		{
			const size_t nPts = m_vPtInfo.size();
			m_vPtInfo.resize(nPts + 1);
			PointInfo& ptInfo = m_vPtInfo.back();
			ptInfo.index = sectionType;
			vector3& coords = ptInfo.coordinates;
			UG_COND_THROW(!(secFile >> coords[0]), "Failed reading x coordinate for point.");
			UG_COND_THROW(!(secFile >> coords[1]), "Failed reading y coordinate for point.");
			UG_COND_THROW(!(secFile >> coords[2]), "Failed reading z coordinate for point.");
			VecScale(coords, coords, m_scaleLength);

			UG_COND_THROW(!(secFile >> ptInfo.diameter), "Failed reading diameter information for point.");
			ptInfo.diameter *= m_scaleLength;

			if (!firstPt)
				ptInfo.connectsTo = nPts - 1;
			else
			{
				ptInfo.connectsTo = (size_t) -1;
				firstPt = false;
			}
		}

		secInfo.end = m_vPtInfo.size() - 1;
	}
	secFile.close();
}



static size_t closest_pt_to_section_position
(
	const SectionInfo& sec,
	number locPos,
	const std::vector<PointInfo> vPts
)
{
	if (locPos <= 0.0)
		return sec.start;

	if (locPos >= 1.0)
		return sec.end;

	const size_t nPts = sec.end - sec.start + 1;
	std::vector<number> vCumulLength(nPts, 0.0);
	for (size_t i = 1; i < nPts; ++i)
		vCumulLength[i] = vCumulLength[i-1] + VecDistance(vPts[sec.start+i-1].coordinates, vPts[sec.start+i].coordinates);

	const number searchLength = locPos * vCumulLength[nPts-1];
	std::vector<number>::iterator it =
		std::lower_bound(vCumulLength.begin(), vCumulLength.end(), searchLength);

	number ub = *it;
	number lb = *--it;

	size_t offset = std::distance(vCumulLength.begin(), it);
	if (searchLength - lb < ub - searchLength)
		return sec.start + offset;

	return sec.start + offset + 1;
}


void NeuronalTopologyImporter::extract_connections_from_hoc(const std::string& text)
{
#if !defined(UG_FOR_VRL) && defined(NETI_ALL_FEATURES)
	// todo: using regex here is quite slow and demands a c++11-compatible compiler; rethink!
	const std::regex connect_match_pattern(".*connect.*");
	const std::regex connect_simple_groups_match_pattern("^connect\\s([^,]+),([^,]+)");
	const std::regex extract_data_match_pattern("(\\w+)\\[?(\\d*|\\w*[-+]?\\d*)\\]?\\((\\d+\\.?\\d*)\\),\\s*(\\w+)\\[?(\\d*|\\w*[-+]?\\d*)\\]?\\((\\d+\\.?\\d*)\\)");
	const std::regex connect_for_match_pattern("for(.*)");
	const std::regex connect_for_extract_data_pattern("for\\s*(\\w+)\\s*=\\s*(\\d+)\\s*,\\s*(\\d+)\\sconnect\\s*(\\w+)\\[?(.*?)\\]?\\s*\\((.*?)\\)\\s*,\\s*(\\w+)\\[?(.*?)\\]?\\s*\\((.*?)\\)\\s*.*");

	const size_t matches_complex = 9;
	const size_t matches_simple = 6;

	// find any line containing "connect"
	const std::sregex_token_iterator end;
	for (std::sregex_token_iterator i(text.cbegin(), text.cend(), connect_match_pattern, 0); i != end; ++i) {
		// trim white space from line
		std::string connect_line = TrimString(*i);

		// find simple connections (i.e., connections like "connect dend[27](0), dend[25](1)", without for loop)
		for (std::sregex_token_iterator l(connect_line.cbegin(), connect_line.cend(), connect_simple_groups_match_pattern, 0); l != end; l++) {
			// extract data from these connection declarations
			const std::string& simple_connect_line = *l;
			std::vector<std::string> connData;
			for (size_t i = 1; i <= matches_simple; i++)
				for (std::sregex_token_iterator m(simple_connect_line.cbegin(), simple_connect_line.cend(), extract_data_match_pattern, i); m != end; m++)
					connData.push_back(*m);

			// name data parts
			std::string& secName1 = connData[0];
			const std::string& secNumber1 = connData[1];
			const number secLoc1 = boost::lexical_cast<number>(connData[2]);
			std::string& secName2 = connData[3];
			const std::string& secNumber2 = connData[4];
			const number secLoc2 = boost::lexical_cast<number>(connData[5]);

			// save connections
			if (!secNumber1.empty())
				secName1 += "[" + secNumber1 + "]";
			std::map<std::string, SectionInfo>::const_iterator secIt1 = m_mSecInfo.find(secName1);
			UG_COND_THROW(secIt1 == m_mSecInfo.end(), "Section '" << secName1 << "' in connection declarations "
				"is not a known section name.");

			if (!secNumber2.empty())
				secName2 += "[" + secNumber2 + "]";
			std::map<std::string, SectionInfo>::const_iterator secIt2 = m_mSecInfo.find(secName2);
			UG_COND_THROW(secIt2 == m_mSecInfo.end(), "Section '" << secName2 << "' in connection declarations "
				"is not a known section name.");

			UG_COND_THROW(secLoc1 != 0.0, "Connection local coordinate for first point is not zero."
				"This case is not implemented.");

			size_t connectingIndex = secIt1->second.start;
			size_t parentIndex = closest_pt_to_section_position(secIt2->second, secLoc2, m_vPtInfo);
			m_vPtInfo[connectingIndex].connectsTo = parentIndex;
		}

		// connections in for statement (i.e., connections like "for i = 29, 30 connect dend[i](0), dend[28](1)")
		for (std::sregex_token_iterator j(connect_line.cbegin(), connect_line.cend(), connect_for_match_pattern); j != end; j++) {
			// extract data from these connection declarations
			std::string for_connect_line = *j;
			std::vector<std::string> connData;
			for (size_t i = 1; i <= matches_complex; i++)
				for (std::sregex_token_iterator z(for_connect_line.cbegin(), for_connect_line.cend(), connect_for_extract_data_pattern, i); z != end; z++)
					connData.push_back(*z);
			UG_COND_THROW(connData.size() != matches_complex, "Connection regex did not match the for loop pattern!\n"
				"Check if filetype specification of hoc-files changed!");

			// name data parts
			const std::string& loopVar = connData[0];
			const int loopStart = boost::lexical_cast<int>(connData[1]);
			const int loopEnd = boost::lexical_cast<int>(connData[2]);
			const std::string& secName1 = connData[3];
			const std::string& secNumber1 = connData[4];
			const number secLoc1 = boost::lexical_cast<number>(connData[5]);
			const std::string& secName2 = connData[6];
			const std::string& secNumber2 = connData[7];
			const number secLoc2 = boost::lexical_cast<number>(connData[8]);

			// treat variables (and +/- operations) in section numbers
			std::vector<std::string> ops1Minus = TokenizeString(secNumber1, '-');
			std::vector<std::string> ops2Minus = TokenizeString(secNumber2, '-');
			std::vector<std::string> ops1Plus = TokenizeString(secNumber1, '+');
			std::vector<std::string> ops2Plus = TokenizeString(secNumber2, '+');

			std::stringstream ss;

			// for loop sections
			for (int index = loopStart; index < loopEnd+1; ++index)
			{
				// replace possibly variable arguments
				int newSecNum1 = 0;
				bool replaceSecNumber1 = false;
				if (secNumber1 == loopVar) // arg1 is the variable
				{
					newSecNum1 += index;
					replaceSecNumber1 = true;
				}
				else if (ops1Minus.size() == 2)  // arg1 is loopVar - x
				{
					if (ops1Minus[0] == loopVar)
						newSecNum1 += index;
					else
						newSecNum1 += boost::lexical_cast<int>(ops1Minus[0]);

					if (ops1Minus[1] == loopVar)
						newSecNum1 -= index;
					else
						newSecNum1 -= boost::lexical_cast<int>(ops1Minus[1]);

					replaceSecNumber1 = true;
				}
				else if (ops1Plus.size() == 2)  // arg1 is var + x
				{
					if (ops1Plus[0] == loopVar)
						newSecNum1 += index;
					else
						newSecNum1 += boost::lexical_cast<int>(ops1Plus[0]);

					if (ops1Plus[1] == loopVar)
						newSecNum1 += index;
					else
						newSecNum1 += boost::lexical_cast<int>(ops1Plus[1]);

					replaceSecNumber1 = true;
				}

				ss << secName1;
				if (replaceSecNumber1)
				{
					if (newSecNum1)
						ss << "[" << newSecNum1 << "]";
				}
				else if (!secNumber1.empty())
					ss << "[" << secNumber1 << "]";
				const std::string newSecName1 = ss.str();

				ss.str("");
				ss.clear();


				int newSecNum2 = 0;
				bool replaceSecNumber2 = false;
				if (secNumber2 == loopVar) // arg2 is the variable
				{
					newSecNum2 += index;
					replaceSecNumber2 = true;
				}
				else if (ops2Minus.size() == 2)  // arg2 is loopVar - x
				{
					if (ops2Minus[0] == loopVar)
						newSecNum2 += index;
					else
						newSecNum2 += boost::lexical_cast<int>(ops2Minus[0]);

					if (ops2Minus[1] == loopVar)
						newSecNum2 -= index;
					else
						newSecNum2 -= boost::lexical_cast<int>(ops2Minus[1]);

					replaceSecNumber2 = true;
				}
				else if (ops2Plus.size() == 2)  // arg2 is loopVar + x
				{
					if (ops2Plus[0] == loopVar)
						newSecNum2 += index;
					else
						newSecNum2 += boost::lexical_cast<int>(ops2Plus[0]);

					if (ops2Plus[1] == loopVar)
						newSecNum2 += index;
					else
						newSecNum2 += boost::lexical_cast<int>(ops2Plus[1]);

					replaceSecNumber2 = true;
				}

				ss << secName2;
				if (replaceSecNumber2)
				{
					if (newSecNum2)
						ss << "[" << newSecNum2 << "]";
				}
				else if (!secNumber2.empty())
					ss << "[" << secNumber2 << "]";
				const std::string newSecName2 = ss.str();

				ss.str("");
				ss.clear();

				// save connections
				std::map<std::string, SectionInfo>::const_iterator secIt1 = m_mSecInfo.find(newSecName1);
				UG_COND_THROW(secIt1 == m_mSecInfo.end(), "Section '" << newSecName1 << "' in connection declarations "
					"is not a known section name.");

				std::map<std::string, SectionInfo>::const_iterator secIt2 = m_mSecInfo.find(newSecName2);
				UG_COND_THROW(secIt2 == m_mSecInfo.end(), "Section '" << newSecName2 << "' in connection declarations "
					"is not a known section name.");

				UG_COND_THROW(secLoc1 != 0.0, "Connection local coordinate for first point is not zero."
					"This case is not implemented.");

				size_t connectingIndex = secIt1->second.start;
				size_t parentIndex = closest_pt_to_section_position(secIt2->second, secLoc2, m_vPtInfo);
				m_vPtInfo[connectingIndex].connectsTo = parentIndex;
			}
		}
	}

	return;
#endif

	UG_THROW("Either your compiler does not support C++11 or this feature "
		"has been explicitly switched off at compile-time.");
}


void NeuronalTopologyImporter::extract_connections_from_ngx(const boost::property_tree::ptree& tree)
{
	const boost::property_tree::ptree& connecTree = tree.get_child("ngx.connections");

	BOOST_FOREACH(const boost::property_tree::ptree::value_type& v, connecTree)
	{
		if (v.first == "connection")
		{
			const std::string& secName1 = v.second.get<std::string>("from");
			std::map<std::string, SectionInfo>::const_iterator secIt1 = m_mSecInfo.find(secName1);
			UG_COND_THROW(secIt1 == m_mSecInfo.end(), "Section '" << secName1 << "' in connection declarations "
				"is not a known section name.");

			const std::string& secName2 = v.second.get<std::string>("to");
			std::map<std::string, SectionInfo>::const_iterator secIt2 = m_mSecInfo.find(secName2);
			UG_COND_THROW(secIt2 == m_mSecInfo.end(), "Section '" << secName2 << "' in connection declarations "
				"is not a known section name.");

			const number secLoc1 = v.second.get<number>("from__loc");
			const number secLoc2 = v.second.get<number>("to__loc");

			UG_COND_THROW(secLoc1 != 0.0, "Connection local coordinate for first point is not zero."
				"This case is not implemented.");

			const size_t connectingIndex = secIt1->second.start;
			const size_t parentIndex = closest_pt_to_section_position(secIt2->second, secLoc2, m_vPtInfo);
			m_vPtInfo[connectingIndex].connectsTo = parentIndex;
		}
	}
}


void NeuronalTopologyImporter::extract_connections_from_txt(const std::string& fileName)
{
	//UG_LOGN("Reading connections from file '" << fconnections_n << "'.");
	std::ifstream connFile(fileName.c_str());
	UG_COND_THROW(!connFile.is_open(), "Connections file could not be opened!");

	std::string secName1;
	std::string secName2;
	number secLoc1;
	number secLoc2;
	while (true)
	{
		if (!(connFile >> secName1))
			break;  // regular end of file

		UG_COND_THROW(!(connFile >> secName2), "Could not read second connecting section name.");
		UG_COND_THROW(!(connFile >> secLoc1), "Could not read location of connection on child branch.");
		UG_COND_THROW(!(connFile >> secLoc2), "Could not read location of connection on parent branch.");

		std::map<std::string, SectionInfo>::const_iterator secIt1 = m_mSecInfo.find(secName1);
		UG_COND_THROW(secIt1 == m_mSecInfo.end(), "Section '" << secName1 << "' in connection declarations "
			"is not a known section name.");

		std::map<std::string, SectionInfo>::const_iterator secIt2 = m_mSecInfo.find(secName2);
		UG_COND_THROW(secIt2 == m_mSecInfo.end(), "Section '" << secName2 << "' in connection declarations "
			"is not a known section name.");

		UG_COND_THROW(secLoc1 != 0.0, "Connection local coordinate for child section is not zero."
			"This case is not implemented.");

		const size_t connectingIndex = secIt1->second.start;
		const size_t parentIndex = closest_pt_to_section_position(secIt2->second, secLoc2, m_vPtInfo);
		m_vPtInfo[connectingIndex].connectsTo = parentIndex;
		m_vPtInfo[connectingIndex].isDuplicate = true;
	}
	connFile.close();
}



void NeuronalTopologyImporter::extract_synapses_from_hoc(const std::string& text)
{
#if !defined(UG_FOR_VRL) && defined(NETI_ALL_FEATURES)
	// bi-exponential synapses
	// todo: using regex here is quite slow and demands a c++11-compatible compiler; rethink!
	// todo: we may also need to handle exponential notation in numbers matching
	const std::regex synapse_connectivity_pattern_exp2syn("(objectvar[^]*?Exp2Syn[^]*?NetCon.*\\))"); // sg
	const std::regex synapse_extract_information("\\((\\d+\\.\\d+)\\)[^]*tau1\\s=\\s(\\d+\\.\\d+)[^]*tau2\\s=\\s(\\d+\\.\\d+)[^]*e\\s=\\s(\\d+\\.\\d+)[^]*NetCon\\(.*?\\((\\d+\\.\\d+)\\).*?(\\d+\\.\\d+).*?(\\d+\\.\\d+).*?(\\d+\\.\\d+)\\)");
	const std::regex synapse_extract_from_compartment_and_synapse_name("(\\w+\\[?\\d*\\]?)\\s\\w+\\s=\\snew\\s*NetCon\\([^,]+, ([^,]+), .*?)\\)");
	const std::regex synapse_extract_to_compartment_and_synapse_name("(\\w+\\[?\\d*\\]?)\\s\\w+\\s=\\snew\\s*Exp2Syn\\([^,]+, ([^,]+), .*?)\\)");

	m_vSynInfo.clear();
	const std::sregex_token_iterator end;
	for (std::sregex_token_iterator i(text.cbegin(), text.cend(), synapse_connectivity_pattern_exp2syn, 1); i != end; ++i) {
		const std::string& synapseHunk = *i;

		// extract synapse physiological and location/connection properties
		std::vector<std::string> extractedSynInfos;
		for (size_t k = 1; k <= 8; k++)
			for (std::sregex_token_iterator j(synapseHunk.cbegin(), synapseHunk.cend(), synapse_extract_information, k); j != end; j++)
				extractedSynInfos.push_back(*j);

		for (size_t k = 1; k <= 2; k++)
			for (std::sregex_token_iterator j(synapseHunk.cbegin(), synapseHunk.cend(), synapse_extract_from_compartment_and_synapse_name, k); j != end; j++)
				extractedSynInfos.push_back(*j);

		for (size_t k = 1; k <= 2; k++)
			for (std::sregex_token_iterator j(synapseHunk.cbegin(), synapseHunk.cend(), synapse_extract_to_compartment_and_synapse_name, k); j != end; j++)
				extractedSynInfos.push_back(*j);

		//	populate synapse information object
		m_vSynInfo.resize(m_vSynInfo.size() + 1);
		SynapseInfo& s = m_vSynInfo.back();
		s.type = synapse_handler::EXP2_POST_SYNAPSE;

		s.secName_pre = extractedSynInfos[8];  // from section
		UG_COND_THROW(m_mSecInfo.find(s.secName_pre) == m_mSecInfo.end(),
			"Synapse presynaptic section '" << s.secName_pre << "' is not a known section name.");
		s.secName_post = extractedSynInfos[10]; // to section
		UG_COND_THROW(m_mSecInfo.find(s.secName_post) == m_mSecInfo.end(),
			"Synapse postsynaptic section '" << s.secName_post << "' is not a known section name.");

		s.secLocCoord_pre = boost::lexical_cast<number>(extractedSynInfos[0]);
		s.secLocCoord_post = boost::lexical_cast<number>(extractedSynInfos[4]);

		s.tau1 = m_scaleTime * boost::lexical_cast<number>(extractedSynInfos[1]);
		s.tau2 = m_scaleTime * boost::lexical_cast<number>(extractedSynInfos[2]);
		s.e = m_scalePotential * boost::lexical_cast<number>(extractedSynInfos[3]);
		s.threshold = m_scalePotential * boost::lexical_cast<number>(extractedSynInfos[5]);
		s.delay = m_scaleTime * boost::lexical_cast<number>(extractedSynInfos[6]);
		s.gMax = m_scaleConductance * boost::lexical_cast<number>(extractedSynInfos[7]);
	}

	// alpha synapses
	const std::regex synapse_connectivity_pattern_alphasyn("(objectvar[^]*?new AlphaSynapse\\(\\d*\\.\\d*\\)[^]*?onset[^]*?tau[^]*?gmax[^]*?e.*[^])");
	const std::regex alphasynapse_extract_pattern("objectvar (.*)[^](.*?) .*? = new AlphaSynapse\\((.*?)\\)[^]*?onset = (\\d*\\.?\\d*)[^]*?tau = (\\d*\\.?\\d*)[^]*?gmax = (\\d*\\.?\\d*)[^]*?e = (\\d*\\.?\\d*).*[^]");

	for (std::sregex_token_iterator a(text.cbegin(), text.cend(), synapse_connectivity_pattern_alphasyn, 1); a != end; ++a) {
		std::string synapseHunk = *a;

		std::vector<std::string> extractedSynInfos;
		for (size_t k = 1; k <= 7; k++)
			for (std::sregex_token_iterator b(synapseHunk.cbegin(), synapseHunk.cend(), alphasynapse_extract_pattern, k); b != end; ++b)
				extractedSynInfos.push_back(*b);

		m_vSynInfo.resize(m_vSynInfo.size() + 1);
		SynapseInfo& s = m_vSynInfo.back();
		s.type = synapse_handler::ALPHA_POST_SYNAPSE;

		s.secName_post = extractedSynInfos[1];
		UG_COND_THROW(m_mSecInfo.find(s.secName_post) == m_mSecInfo.end(),
			"Synapse postsynaptic section '" << s.secName_post << "' is not a known section name.");

		s.secLocCoord_post = boost::lexical_cast<number>(extractedSynInfos[2]);

		s.onset = m_scaleTime * boost::lexical_cast<number>(extractedSynInfos[3]);
		s.tau1 = m_scaleTime * boost::lexical_cast<number>(extractedSynInfos[4]);
		s.gMax = m_scaleConductance * boost::lexical_cast<number>(extractedSynInfos[5]);
		s.e = m_scalePotential * boost::lexical_cast<number>(extractedSynInfos[6]);
	}

	return;
#endif

	UG_THROW("Either your compiler does not support C++11 or this feature "
		"has been explicitly switched off at compile-time.");
}


void NeuronalTopologyImporter::extract_synapses_from_ngx(const boost::property_tree::ptree& tree)
{
	using boost::property_tree::ptree;

	const ptree& exp2SynTree = tree.get_child("ngx.exp2synapses");
	BOOST_FOREACH(const ptree::value_type& v, exp2SynTree)
	{
		if (v.first == "Exp2Syn")
		{
			//	populate synapse information object
			m_vSynInfo.resize(m_vSynInfo.size() + 1);
			SynapseInfo& s = m_vSynInfo.back();
			s.type = synapse_handler::EXP2_POST_SYNAPSE; // type of synapse

			s.secName_pre = v.second.get<std::string>("from");  // from section
			UG_COND_THROW(m_mSecInfo.find(s.secName_pre) == m_mSecInfo.end(),
				"Synapse presynaptic section '" << s.secName_pre << "' is not a known section name.");
			s.secName_post = v.second.get<std::string>("to"); // to section
			UG_COND_THROW(m_mSecInfo.find(s.secName_post) == m_mSecInfo.end(),
				"Synapse postsynaptic section '" << s.secName_post << "' is not a known section name.");

			s.secLocCoord_pre = v.second.get<number>("from__loc");  // from location
			s.secLocCoord_post = v.second.get<number>("to__loc");  // to location

			s.tau1 = m_scaleTime * v.second.get<number>("tau1"); // tau1
			s.tau2 = m_scaleTime * v.second.get<number>("tau2"); // tau2
			s.e = m_scalePotential * v.second.get<number>("e"); // reversal potential
			s.threshold = m_scalePotential * v.second.get<number>("threshold"); // threshold
			s.delay = m_scaleTime * v.second.get<number>("delay"); // delay
			s.gMax = m_scaleConductance * v.second.get<number>("gmax"); // weight
		}
	}

	const ptree& alphaSynTree = tree.get_child("ngx.alphasynapses");
	BOOST_FOREACH(const ptree::value_type& v, alphaSynTree)
	{
		if (v.first == "AlphaSyn")
		{
			m_vSynInfo.resize(m_vSynInfo.size() + 1);
			SynapseInfo& s = m_vSynInfo.back();
			s.type = synapse_handler::ALPHA_POST_SYNAPSE;  // type of synapse

			s.secName_post = v.second.get<std::string>("from"); // section of input
			UG_COND_THROW(m_mSecInfo.find(s.secName_post) == m_mSecInfo.end(),
				"Synapse postsynaptic section '" << s.secName_post << "' is not a known section name.");

			s.secLocCoord_post = v.second.get<number>("from__loc"); // location on section

			s.onset = m_scaleTime * v.second.get<number>("onset"); // onset
			s.tau1 = m_scaleTime * v.second.get<number>("tau"); // tau1=tau for alpha syn
			s.gMax = m_scaleConductance * v.second.get<number>("gmax"); // weight=gmax for alpha syn
			s.e = m_scalePotential * v.second.get<number>("e"); // reversal potential
		}
	}
}


void NeuronalTopologyImporter::extract_synapses_from_txt(const std::string& fileName)
{
	//UG_LOGN("Reading synapses from file '" << fileName << "'.");
	std::ifstream synFile(fileName.c_str());
	UG_COND_THROW(!synFile.is_open(), "Synapses file could not be opened!");

	std::map<Edge*, std::vector<std::pair<Edge*, size_t> > > pre_post_mapping;

	number from_loc;
	int type;
	int dummy;
	while (true)
	{
		if (!(synFile >> from_loc))
			break;  // regular end of file

		// create new SynapseInfo and fill
		m_vSynInfo.resize(m_vSynInfo.size() + 1);
		SynapseInfo& s = m_vSynInfo.back();

		s.secLocCoord_pre = from_loc;
		UG_COND_THROW(!(synFile >> s.secLocCoord_post), "Could not read location of post-synapse.");

		UG_COND_THROW(!(synFile >> type), "Could not read synapse type.");

		UG_COND_THROW(!(synFile >> s.gMax), "Could not read maximal conductance.");
		s.gMax *= m_scaleConductance;

		UG_COND_THROW(!(synFile >> s.secName_pre), "Could not read pre-synapse section name.");
		UG_COND_THROW(m_mSecInfo.find(s.secName_pre) == m_mSecInfo.end(),
			"Synapse presynaptic section '" << s.secName_pre << "' is not a known section name.");

		UG_COND_THROW(!(synFile >> s.secName_post), "Could not read post-synapse section name.");
		UG_COND_THROW(m_mSecInfo.find(s.secName_post) == m_mSecInfo.end(),
			"Synapse postsynaptic section '" << s.secName_post << "' is not a known section name.");

		// todo: these two remain unused; think about removing them from the file specification!
		UG_COND_THROW(!(synFile >> dummy), "Could not read pre-synapse index.");
		UG_COND_THROW(!(synFile >> dummy), "Could not read post-synapse index.");

		switch (type)
		{
			case 0:  // exp2
			{
				s.type = synapse_handler::EXP2_POST_SYNAPSE; // type of synapse
				s.tau1 = m_scaleTime * 0.2;
				s.tau2 = m_scaleTime * 1.7;
				s.e = m_scalePotential * 0.0;
				s.threshold = m_scalePotential * -10.0;
				break;
			}
			case 1:  // alpha
			{
				s.type = synapse_handler::ALPHA_POST_SYNAPSE;  // type of synapse
				s.tau1 = m_scaleTime * 1.7;
				s.e = m_scalePotential * 0.0;
				s.onset = m_scaleTime * 0.0;
				break;
			}
			default:
				UG_THROW("Unknown synapse type: " << type << ".");
		}
	}
	synFile.close();
}




bool NeuronalTopologyImporter::generate_grid()
{
	// grid
	m_grid = SmartPtr<MultiGrid>(new MultiGrid());

	// subset handler
	m_sh = SmartPtr<MGSubsetHandler>(new MGSubsetHandler(*m_grid));
	MGSubsetHandler& sh = *m_sh;

	// coordinates attachment (global attachment) and accessor
	m_grid->attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition3> aaPos(*m_grid, aPosition);

	m_aDiameterAttachment = GlobalAttachments::attachment<ANumber>("diameter");
	m_grid->attach_to_vertices(m_aDiameterAttachment);
	m_aaDiameter = Grid::VertexAttachmentAccessor<ANumber>(*m_grid, m_aDiameterAttachment);

	// synapse attachment
	const size_t nSyn = m_vSynInfo.size();
	if (nSyn)
	{
		m_aSynapse = GlobalAttachments::attachment<ASynapse>("synapses");
		m_grid->attach_to_edges(m_aSynapse);
		m_aaSynapse = Grid::EdgeAttachmentAccessor<ASynapse>(*m_grid, m_aSynapse);
	}

	// store pointer to vertices
	std::vector<Vertex*> vVrts;

	// create vertices and assign coordinates
	UG_DLOGN(NETIGrid, 6, "Generating vertices.");
	for (std::vector<PointInfo>::const_iterator it = m_vPtInfo.begin(); it != m_vPtInfo.end(); ++it) {
		// create and save grid vertex
		Vertex* vertex = *(m_grid->create<RegularVertex>());
		vVrts.push_back(vertex);

		// set coordinates and diameter
		aaPos[vertex] = it->coordinates;
		m_aaDiameter[vertex] = it->diameter;

		// assign to subset
		sh.assign_subset(vertex, it->index);
	}
	UG_DLOGN(NETIGrid, 6, "Vertex generation completed.");

	UG_DLOGN(NETIGrid, 6, "Generating edges.");
	const size_t nPts = m_vPtInfo.size();
	for (size_t p = 0; p < nPts; ++p)
	{
		const PointInfo& ptInfo = m_vPtInfo[p];
		const size_t& conn = ptInfo.connectsTo;
		if (conn == (size_t) -1)
			continue;

		if (ptInfo.isDuplicate)  // connection point with duplicate
		{
			MergeVertices(*m_grid, vVrts[conn], vVrts[p]);
			vVrts[p] = vVrts[conn];
		}
		else
		{
			EdgeDescriptor ed(vVrts[p], vVrts[conn]);
			Edge* edge = *(m_grid->create<RegularEdge>(ed));

			// assign edge to subset
			sh.assign_subset(edge, m_vPtInfo[p].index);  // end vertex of edge determines subset membership
		}
	}
	UG_DLOGN(NETIGrid, 6, "Edge generation completed.");

	// add synapses to edges
	UG_DLOGN(NETIGrid, 6, "Generating synapses.");
	synapse_id id = 0;
	size_t totalAlphaSyns = 0;
	size_t totalExp2Syns = 0;
	for (size_t s = 0; s < nSyn; ++s)
	{
		const SynapseInfo& synInfo = m_vSynInfo[s];

		switch (synInfo.type)
		{
			case cable_neuron::synapse_handler::ALPHA_POST_SYNAPSE:
			{
				Edge* preAndPostEdge = NULL;
				number locCoord = 0.0;
				try {get_synapse_edge_and_local_coords(&preAndPostEdge, locCoord, synInfo.secName_post,
					synInfo.secLocCoord_post, m_vPtInfo, m_mSecInfo, vVrts, *m_grid);}
				UG_CATCH_THROW("Could not determine edge for synapse.");

				// scaling synapse params from Neuron default units (mV, ms, uS) to cable_neuron units (V, s, S)
				IBaseSynapse* presyn = new OnsetPreSynapse(id, locCoord, 1e-3*synInfo.onset, 1e-3*6.0*synInfo.tau1);
				IBaseSynapse* postsyn = new AlphaPostSynapse(id, locCoord, 1e-6*synInfo.gMax, 1e-3*synInfo.tau1, 1e-3*synInfo.e);

				// both pre- and post-synapse go to the post-synaptic edge
				m_aaSynapse[preAndPostEdge].push_back(presyn);
				m_aaSynapse[preAndPostEdge].push_back(postsyn);

				++totalAlphaSyns;

				break;
			}
			case cable_neuron::synapse_handler::EXP2_POST_SYNAPSE:
			{
				Edge* pre = NULL;
				Edge* post = NULL;
				number locCoordPre = 0.0;
				number locCoordPost = 0.0;
				try
				{
					get_synapse_edge_and_local_coords(&pre, locCoordPre, synInfo.secName_pre,
						synInfo.secLocCoord_pre, m_vPtInfo, m_mSecInfo, vVrts, *m_grid);
					get_synapse_edge_and_local_coords(&post, locCoordPost, synInfo.secName_post,
						synInfo.secLocCoord_post, m_vPtInfo, m_mSecInfo, vVrts, *m_grid);
				}
				UG_CATCH_THROW("Could not determine edge for pre- or post-synapse.");

				// scaling synapse params from Neuron default units (mV, ms, uS) to cable_neuron units (V, s, S)
				const number& tau1 = synInfo.tau1;
				const number& tau2 = synInfo.tau2;
				const number tp = (tau1*tau2)/(tau2-tau1) * std::log(tau2/tau1);
				const number duration = tp + 3*tau2;
				IBaseSynapse* presyn = new ThresholdPreSynapse(id, locCoordPre, 1e-3*duration, 1e-3*synInfo.threshold);
				IBaseSynapse* postsyn = new Exp2PostSynapse(id, locCoordPost, 1e-6*synInfo.gMax, 1e-3*tau1, 1e-3*tau2, 1e-3*synInfo.e);

				// both pre- and post-synapse go to the post-synaptic edge
				m_aaSynapse[pre].push_back(presyn);
				m_aaSynapse[post].push_back(postsyn);

				++totalExp2Syns;

				break;
			}
			default:
				UG_LOGN("Encountered unknown synapse type '" << synInfo.type << "' (ignoring).");
		}

		++id;
	}
	UG_DLOGN(NETIGrid, 6, "Synapses generated (#exp2Syn: " << totalExp2Syns
		<< "  #alphaSyn: " << totalAlphaSyns << ").");


	m_grid->enable_hierarchical_insertion(false);  // otherwise, edges might end up in level 1 during RemoveDoubles
	RemoveDoubles<3>(*m_grid, m_grid->vertices_begin(), m_grid->vertices_end(), aPosition, REMOVE_DOUBLES_THRESHOLD);

	UG_DLOGN(NETIGrid, 6, "num_vertices of grid: " << m_grid->num_vertices());
	UG_DLOGN(NETIGrid, 6, "num_edges of grid: " << m_grid->num_edges());

	return true;
}




void NeuronalTopologyImporter::subset_processing()
{
	if (m_extension == "swc")
		subset_processing_SWC();
	else if (m_extension == "hoc")
		subset_processing_HOC();
	else if (m_extension == "ngx")
		subset_processing_NGX();
	else if (m_extension == "txt")
		subset_processing_TXT();
	else if (m_extension == "xml")
		subset_processing_NML();
	else
		UG_DLOGN(NETIGeometry, 1, "UNKNOWN filetype '" << m_extension << "'.");
}



void NeuronalTopologyImporter::subset_processing_SWC()
{
	MGSubsetHandler& sh = *m_sh;

	// rename subsets
	sh.subset_info(0).name = "undefined";
	sh.subset_info(1).name = "soma";
	sh.subset_info(2).name = "axon";
	sh.subset_info(3).name = "dendrite";
	sh.subset_info(4).name = "apical_dendrite";
	sh.subset_info(5).name = "fork";
	sh.subset_info(6).name = "end";
	sh.subset_info(7).name = "custom";
	for (int i = 8; i < sh.num_subsets(); ++i)
		sh.subset_info(i).name = "default";

	// erase empty subsets
	EraseEmptySubsets(sh);

	// assign colors for subsets
	AssignSubsetColors(sh);
}


void NeuronalTopologyImporter::subset_processing_HOC()
{
	MGSubsetHandler& sh = *m_sh;

	// name subsets
	for (std::map<std::string, SectionInfo>::iterator it = m_mSecInfo.begin(); it != m_mSecInfo.end(); ++it)
		sh.subset_info(it->second.index).name = it->first;

	// join subsets
	const size_t num_subsets = sh.num_subsets();
	for (std::vector<std::string>::const_iterator cit = m_joiningCriteria.begin(); cit != m_joiningCriteria.end(); ++cit) {
		int firstHitInd = -1;
		for (int i = 0; (size_t) i < num_subsets; ++i) {
			const std::string& name = sh.subset_info(i).name;
			size_t found = name.find(*cit);
			if (found != std::string::npos) {
				if (firstHitInd == -1)
					firstHitInd = i;
				else
				{
					sh.assign_subset(sh.begin<Vertex>(i, 0), sh.end<Vertex>(i, 0), firstHitInd);
					sh.assign_subset(sh.begin<Edge>(i, 0), sh.end<Edge>(i, 0), firstHitInd);
				}
			}
		}
		if (firstHitInd != -1)
			sh.subset_info(firstHitInd).name = *cit;
	}

	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);
}


void NeuronalTopologyImporter::subset_processing_NGX()
{
	subset_processing_HOC();
}


void NeuronalTopologyImporter::subset_processing_TXT()
{
	MGSubsetHandler& sh = *m_sh;

	// name and colorize subsets, erase empty ones
	if (!m_bWithCellType)
	{
		int si_dend = 0;
		int si_axon = 1;
		int si_soma = 2;
		sh.subset_info(si_dend).name = "Axon";
		sh.subset_info(si_axon).name = "Dendrite";
		sh.subset_info(si_soma).name = "Soma";
	}
	else
	{
		std::vector<std::string> vSectionNames(3);
		vSectionNames[0] = "AXON";
		vSectionNames[1] = "DEND";
		vSectionNames[2] = "SOMA";

		std::vector<std::string> vCellTypeNames(5);
		if (m_networkID == 0)
		{
			vCellTypeNames[0] = "L4_STELLATE";
			vCellTypeNames[1] = "L23_PYRAMIDAL";
			vCellTypeNames[2] = "L5A_PYRAMIDAL";
			vCellTypeNames[3] = "L5B_PYRAMIDAL";
			vCellTypeNames[4] = "L4_STAR_PYRAMIDAL";
		}
		else if (m_networkID == 1)
		{
			vCellTypeNames.resize(5);
			vCellTypeNames[0] = "CA1_PYRAMIDAL";
			vCellTypeNames[1] = "CB_CALBINDIN";
			vCellTypeNames[2] = "CCK_Cholecystokinin";
			vCellTypeNames[3] = "PV_Paralbumin";
			vCellTypeNames[4] = "SOM";
		}
		else
		{
			UG_THROW("Unknown network ID: '" << m_networkID << "'.");
		}

		std::stringstream ss3;
		const int numSectionTypes = 3;
		for (int i = 0; i < numSectionTypes; i++)
		{
			for (int j = 0; j < (int) m_numCellTypes; j++)
			{
				ss3 << vSectionNames[i] << "__" << vCellTypeNames[j];
				sh.subset_info(i*m_numCellTypes + j).name = ss3.str();
				ss3.clear(); ss3.str("");
			}
		}
	}

	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);
}



void NeuronalTopologyImporter::subset_processing_NML()
{
	MGSubsetHandler& sh = *m_sh;
	sh.set_subset_name("soma", 0);
	sh.set_subset_name("axon", 1);
	sh.set_subset_name("dend", 2);
	sh.set_subset_name("apic", 3);
	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);
}




}  // namespace neuronal_topology_importer
}  // namespace cable_neuron
}  // namespace ug

