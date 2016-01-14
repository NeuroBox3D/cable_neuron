/*
 * synapse_distributor.cpp
 *
 *  Created on: 03.11.2014
 *      Author: Lukas Reinhardt
 */

#include "synapse_distributor.h"

// system includes
#include <iostream>	// cout
#include <sstream>	// stringstream
#include <stdlib.h>	// rand()

// boost includes
#include "boost/lexical_cast.hpp"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

// ugcore includes
#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif
#include "lib_grid/global_attachments.h"

// plugin includes
#include "../synapse_handler/grid/synapse_info_io_traits.h"
#include "../synapse_handler/function/types.h"


#define SynapseDistributor_verbose

using namespace std;

namespace ug {
namespace cable_neuron {


////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::SynapseDistributor(string infile, string outfile)
/**
 * This constructor loads a grid from a given inputfile to work on.
 * outfile is the file, which the grid can be saved to.
 *
 * Note: Use this constructor, if you need preset activity timings.
 * Then ignore isActive attribute, only work with start_time and end_time attributes.
 */
SynapseDistributor::SynapseDistributor(string infile, string outfile, bool bRemoveExistingSynapses)
{
	m_bDomBased = false;

	m_InputFile = infile;
	m_OutputFile = outfile;
	m_LastMessage = "";

//	Grid setup
	pm_Grid = new MultiGrid;
	pm_SubsetHandler = new MGSubsetHandler(*pm_Grid);
	if (!LoadGridFromFile(*pm_Grid, *pm_SubsetHandler, infile.c_str()))
		UG_THROW("File '" << infile << "' could not be loaded.");

//	Global Attachment setup
//	Existence of a synapse is represented by a struct Synapse GlobalAttachment
//	Check existence
	if (!GlobalAttachments::is_declared("Synapses")) {
		UG_THROW("GlobalAttachment 'Synapses' not available.");
	}

//	Initialize attachment and attach to grid, if not already present
	m_aSynInfo = GlobalAttachments::attachment<AVSynapse>("Synapses");
	if (!pm_Grid->has_edge_attachment(m_aSynInfo)) {
		pm_Grid->attach_to_edges(m_aSynInfo);
	}

//	Initialize attachment accessors
	m_aaSynInfo = Grid::EdgeAttachmentAccessor<AVSynapse>(*pm_Grid, m_aSynInfo);
	m_aaPosition = Grid::VertexAttachmentAccessor<APosition>(*pm_Grid, aPosition);

//	Initialize state of synapses
	if(bRemoveExistingSynapses)
		clear();

	m_LastMessage = "SynapseDistributor: SynapseDistributor object successfully instantiated.";
	#ifdef SynapseDistributor_verbose
		UG_LOG("SynapseDistributor: SynapseDistributor object successfully instantiated.\n");
	#endif
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::SynapseDistributor(SmartPtr<Domain3d> dom, string outfile)
/**
 * This constructor gets a grid object reference to work on.
 * outfile is the file, which the grid can be saved to.
 *
 * Note: Use this constructor, if you need preset activity timings.
 * Then ignore isActive attribute, only work with start_time and end_time attributes.
 */
SynapseDistributor::SynapseDistributor(SmartPtr<Domain3d> dom, string outfile, bool bRemoveExistingSynapses)
{
	m_bDomBased = true;

	m_InputFile = "";
	m_OutputFile = outfile;
	m_LastMessage = "";

//	Grid setup
	pm_Grid = dom->grid().get();
	pm_SubsetHandler = dom->subset_handler().get(); //new SubsetHandler(*pm_Grid);

//	Global Attachment setup
//	Existence of a synapse is represented by a struct Synapse GlobalAttachment
//	Check existence
	if (!GlobalAttachments::is_declared("Synapses")) {
		UG_THROW("GlobalAttachment 'Synapses' not available.");
	}

//	Initialize attachment and attach to grid, if not already present
	m_aSynInfo = GlobalAttachments::attachment<AVSynapse>("Synapses");
	if (!pm_Grid->has_edge_attachment(m_aSynInfo)) {
		pm_Grid->attach_to_edges(m_aSynInfo);
	}

//	Initialize attachment accessors
	m_aaSynInfo = Grid::EdgeAttachmentAccessor<AVSynapse>(*pm_Grid, m_aSynInfo);
	m_aaPosition = Grid::VertexAttachmentAccessor<APosition>(*pm_Grid, aPosition);

//	Initialize state of synapses
	if(bRemoveExistingSynapses)
		clear();

	m_LastMessage = "SynapseDistributor: SynapseDistributor object successfully instantiated.";
	#ifdef SynapseDistributor_verbose
		UG_LOG("SynapseDistributor: SynapseDistributor object successfully instantiated.\n");
	#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Standard destr3ctor
SynapseDistributor::~SynapseDistributor()
{
	delete pm_Grid;
	delete pm_SubsetHandler;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// SynapseDistributor::CopySynapsesFromParentToChild(Edge* parent)
/**
 * Copies synapses from parent edge to corresponding child edge
 */
void SynapseDistributor::CopySynapsesFromParentToChild(Edge* parent)
{
	number localBaseCoord = 0.0;
	size_t childEdgeIndex = 0;

//	Iterate over all synapses of the parent edge
	for(size_t i = 0; i < m_aaSynInfo[parent].size(); ++i)
	{
	//	Determine index of child edge, which inherits the synapse, and get new local base coords for the child edge
		if(m_aaSynInfo[parent][i].m_locCoords < 0.5)
		{
			localBaseCoord = 2*m_aaSynInfo[parent][i].m_locCoords;
			childEdgeIndex = 0;
		}
		else{
			localBaseCoord = 2*m_aaSynInfo[parent][i].m_locCoords - 1;
			childEdgeIndex = 1;
		}

	//	Get correct child edge and the parent's synapse
		Edge* child = pm_Grid->get_child_edge(parent, childEdgeIndex);
		synapse_handler::SynapseInfo syn = m_aaSynInfo[parent][i];

	//	Copy synapse to child
		m_aaSynInfo[child].push_back(syn);

	//	Modify the local base coord and associated edge of synapse
		m_aaSynInfo[child].back().m_locCoords = localBaseCoord;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
// SynapseDistributor::CopySynapsesToAllLevels()
/**
 * Passes on synapse attachments to all refinement levels
 */
void SynapseDistributor::CopySynapsesToAllLevels()
{
//	Iterate over all grid levels
	for(size_t i = 0; i < pm_Grid->num_levels(); ++i)
	{

#ifdef UG_PARALLEL
	// 	copy from vertical masters to vertical slaves
		typedef GridLayoutMap::Types<Edge>::Layout layout_type;
		DistributedGridManager& dgm = *pm_Grid->distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
		pcl::InterfaceCommunicator<layout_type> icom;

		ComPol_CopyAttachment<layout_type, AVSynapse> compolCopySynapses(*pm_Grid, m_aSynInfo);
		icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopySynapses);
		icom.communicate();
#endif

	//	Iterate over all edges of level i
		for(EdgeIterator eIter = pm_Grid->begin<Edge>(i); eIter != pm_Grid->end<Edge>(i); ++eIter)
		{
			Edge* e = *eIter;

			if(m_aaSynInfo[e].size() > 0 && pm_Grid->has_children(e))
				CopySynapsesFromParentToChild(e);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::remove_synapse
/**
 * Removes a Synapse from Edge e.
 */
void SynapseDistributor::remove_synapse(Edge* e)
{
	if(m_aaSynInfo[e].size() > 0)
	{
		m_aaSynInfo[e].pop_back();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::remove_all_synapses
/**
 * Removes all Synapses from Edge e.
 */
void SynapseDistributor::remove_all_synapses(Edge* e)
{
	m_aaSynInfo[e].clear();
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::clear()
/**
 * Removes all synapses from all the edges of the grid.
 */
void SynapseDistributor::clear()
{
	for(size_t i = 0; i < pm_Grid->num_levels(); ++i)
	{
		for(EdgeIterator eIter = pm_Grid->begin<Edge>(0); eIter != pm_Grid->end<Edge>(0) ; eIter++)
		{
			Edge* e = *eIter;
			remove_all_synapses(e);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::clear()
/**
 * Removes all synapses from the subset.
 */
void SynapseDistributor::clear(int subsetIndex)
{
	for(size_t i = 0; i < pm_Grid->num_levels(); ++i)
	{
		for(EdgeIterator eIter = pm_SubsetHandler->begin<Edge>(subsetIndex, 0); eIter != pm_SubsetHandler->end<Edge>(subsetIndex, 0); eIter++)
		{
			Edge* e = *eIter;
			remove_all_synapses(e);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::place_synapse
/**
 * Places a Synapse on Edge e.
 */
void SynapseDistributor::place_synapse(Edge* e)
{
	/// TODO: rand() and RAND_MAX seem to be platform-dependent
	number localCoord = static_cast<number>(rand())  / RAND_MAX; //localCoord factor 0 means e[0], 1 means e[1]

//	Initialize an alpha synapse with dummy onset and end activity time
//	(end time is determined by alpha synapse specific parameter m_tau * 6)

	synapse_handler::SynapseInfo syn;
	syn.m_locCoords = localCoord;
	syn.m_type = synapse_handler::ALPHA_SYNAPSE;
	syn.m_onset = 0.0;				// to be set by set_activation_timing method
	syn.m_tau = 0.0;				// to be set by set_activation_timing method
	syn.m_gMax = 6e-4; 				// = 600pS
	syn.m_vRev = 0.0;

//	Add synapse to edge
	m_aaSynInfo[e].push_back(syn);

//	Log
	stringstream msg;
	msg<<"Synapse placed at edge_"<<e<<", Synapse is"
		<<((localCoord<0.5)?" near source node and ": " near destination node and ")<<endl;

	m_LastMessage =  msg.str();
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::set_activation_timing(number start_time, number duration, number start_time_dev, number duration_dev)
/**
 * Sets electrode activity
 *
 * start_time		: average startvalue
 * duration			: average duration of activity
 * start_time_dev	: deviation of startvalue
 * duration_dev		: deviation of duration
 */
void SynapseDistributor::set_activation_timing(number start_time, number duration, number start_time_dev, number duration_dev)
{
//	Activity timing setup
//	############## Random normal distribution
	boost::mt19937 rng;
	rng.seed(time(NULL));

	boost::normal_distribution<number> start_dist(start_time, start_time_dev);
	boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_start(rng, start_dist);

	boost::normal_distribution<number> duration_dist(duration, duration_dev);
	boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_duration(rng, duration_dist);
//	##################

//	Initialize with time dependency
// 	INFO: * end time is determined by start time AND alpha synapse specific parameter tau: t_end = t_start + m_tau * 6
//		  * correspondingly, m_tau is used in methods where end time is needed
	for(EdgeIterator eIter = pm_Grid->begin<Edge>(0); eIter != pm_Grid->end<Edge>(0) ; eIter++)
	{
		Edge* e = *eIter;

		for(size_t i = 0; i < m_aaSynInfo[e].size(); ++i)
		{
			number t_start = var_start();

			if(t_start < 0)
				t_start = 0;

			number dur = var_duration();
			dur = dur < 0 ? 0 : dur;

			number t_end = t_start + abs( var_duration());

			if(t_start > t_end)
				UG_THROW("ERROR in SynapseDistributor constructor: Synapse activity start time > end time.");

			m_aaSynInfo[e][i].m_onset = t_start;
			m_aaSynInfo[e][i].m_tau = dur / 6.0;
		}
	}

//	Handle multigrid transfer
	if(m_bDomBased)
		CopySynapsesToAllLevels();
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::place_synapses(vector<Edge*> vEdges)
/**
 * Wrapper for placing synapses in the grid on specified edges uniformly
 */
void SynapseDistributor::place_synapses_uniform(vector<Edge*> vEdges, size_t numSynapses)
{
//	Determine discrete random distribution of given edges by their lengths
	vector<number> probs;
	number totLength = 0.0;

	for(size_t i = 0; i < vEdges.size(); ++i)
	{
		Edge* e = vEdges[i];
		number edgeLength = EdgeLength(e, m_aaPosition);

	//	Push back current edge length as probability (still to be divided by total length of all edges)
		probs.push_back(edgeLength);

		totLength += edgeLength;
	}

//	Scale probabilities
	VecScale(probs, probs, 1./totLength);

//	Define random int generator for random access of edges vector
	boost::random::mt19937 gen;
#ifdef UG_PARALLEL
	gen.seed(pcl::ProcRank());
#endif
	boost::random::discrete_distribution<> dist(probs.begin(), probs.end());
	boost::variate_generator<boost::random::mt19937, boost::random::discrete_distribution<> > randomIndex(gen, dist);

//	Distribute specified number of synapses along the coarse grid edges randomly (uniformly s.t. individual edge lengths)
	size_t i = 0;
	while(i < numSynapses)
	{
		place_synapse(vEdges[randomIndex()]);
		i++;
	}

//	Handle multigrid transfer
	if(m_bDomBased)
		CopySynapsesToAllLevels();
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::place_synapses_uniform(int si)

/**
 * Places synapses in the grid on edges of given subset, uniform distributed according to set density.
 */
void SynapseDistributor::place_synapses_uniform(int si, size_t numSynapses)
{
//	Save subset coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_SubsetHandler->begin<Edge>(si, 0),
										 pm_SubsetHandler->end<Edge>(si, 0));

//	Call wrapper
	place_synapses_uniform(vEdges, numSynapses);
}

/**
 * Places synapses in the grid on edges, uniform distributed according to set density.
 */
void SynapseDistributor::place_synapses_uniform(size_t numSynapses)
{
//	Save coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_Grid->begin<Edge>(0), pm_Grid->end<Edge>(0));

//	Call wrapper
	place_synapses_uniform(vEdges, numSynapses);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::place_synapses_uniform(const char* subset, number density)

/**
 * Places synapses in the grid on edges of given subset, uniformly distributed according to set density.
 * Specification of density in [1/m]
 */
void SynapseDistributor::place_synapses_uniform(const char* subset, number density)
{
//	Get subset index from subset name
	int si = pm_SubsetHandler->get_subset_index(subset);

//	Get total length of given subset
	number length = this->get_subset_length(si);

//	Get number of synapses to place from density value
	size_t numSynapses = (size_t)(length * density);

//	Save subset coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_SubsetHandler->begin<Edge>(si, 0),
										 pm_SubsetHandler->end<Edge>(si, 0));

//	Call wrapper
	place_synapses_uniform(vEdges, numSynapses);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::place_synapses(vector<size_t> distr)
/**
 * Places synapses according to set density and given distribution in the single subsets.
 */
void SynapseDistributor::place_synapses(vector<number> distr, size_t numSynapses)
{
//	check validity of input
	number s = 0;

	for(size_t i=0; i<distr.size(); ++i)
	{
		s += distr[i];
	}

	if( abs(s - 1) > 1e-5 || (int)distr.size() != pm_SubsetHandler->num_subsets()) { //sum of probabilities has to be 1
		UG_THROW("SynapseDistributor::place_synapses(vector<number> distr): Specified distribution incorrect. Probabilities don't sum up to 1 or distr.size() != mg.num_subsets().");
	}

//	Placing synapses
	for(size_t i=0; i<distr.size(); ++i)
	{
		size_t numSynPerSubset = numSynapses*distr[i];

		vector<Edge*> vEdges = vector<Edge*>(pm_SubsetHandler->begin<Edge>(i, 0), pm_SubsetHandler->end<Edge>(i, 0));

		place_synapses_uniform(vEdges, numSynPerSubset);
	}

//	Handle multigrid transfer
	if(m_bDomBased)
		CopySynapsesToAllLevels();
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::degenerate_uniform(vector<Edge*> vEdges, size_t numSynapses)
/**
 * Wrapper for degenerate synapses in the grid on specified edges uniformly
 */
void SynapseDistributor::degenerate_uniform(vector<Edge*> vEdges, size_t numSynapses)
{
//	Determine number of synapses in the given vector of edges
	if(numSynapses > num_synapses(vEdges, false, 0.0))
		UG_THROW("SynapseDistributor::degenerate_uniform(vector<Edge*> vEdges, size_t numSynapses): Cannot degenerate more synapses than placed.");

//	Determine discrete random distribution of given edges by their lengths
	vector<number> probs;
	number totLength = 0.0;

	for(size_t i = 0; i < vEdges.size(); ++i)
	{
		Edge* e = vEdges[i];
		number edgeLength = EdgeLength(e, m_aaPosition);

	//	Push back current edge length as probability (still to be divided by total length of all edges)
		probs.push_back(edgeLength);

		totLength += edgeLength;
	}

//	Scale probabilities
	VecScale(probs, probs, 1./totLength);

//	Define random int generator for random access of edges vector
	boost::random::mt19937 gen;
	boost::random::discrete_distribution<> dist(probs.begin(), probs.end());
	boost::variate_generator<boost::random::mt19937, boost::random::discrete_distribution<> > randomIndex(gen, dist);

//	Distribute specified number of synapses along the coarse grid edges randomly (uniformly s.t. individual edge lengths)
	size_t i = 0;
	while(i < numSynapses)
	{
		size_t rID = randomIndex();

		if(m_aaSynInfo[vEdges[rID]].size() > 0)
			m_aaSynInfo[vEdges[rID]].pop_back();
		else
			continue;

		i++;
	}

//	Handle multigrid transfer
	if(m_bDomBased)
		CopySynapsesToAllLevels();
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::degenerate_uniform(number p, size_t numSynapses)
/**
 * Removes a percentage p of synapses from the whole grid
 */
void SynapseDistributor::degenerate_uniform(number p)
{
	if(p < 0 || p > 1)
		UG_THROW("SynapseDistributor::degenerate_uniform(number p, size_t numSynapses): Specified percentage incorrect.");

//	Save coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_Grid->begin<Edge>(0), pm_Grid->end<Edge>(0));

//	Determine number of synapses to be removed from subset
	size_t numSynToRemove = p * num_synapses(vEdges, false, 0.0);

	degenerate_uniform(vEdges, numSynToRemove);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::degenerate_uniform(number p, size_t numSynapses, int si)
/**
 * Removes a percentage p of synapses from subset subsetIndex
 */
void SynapseDistributor::degenerate_uniform(number p, int si)
{
	if(p < 0 || p > 1)
		UG_THROW("SynapseDistributor::degenerate_uniform(number p, int subsetIndex): Specified percentage incorrect.");

//	Save subset coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_SubsetHandler->begin<Edge>(si, 0), pm_SubsetHandler->end<Edge>(si, 0));

//	Determine number of synapses to be removed from subset
	size_t numSynToRemove = p * num_synapses(vEdges, false, 0.0);

	degenerate_uniform(vEdges, numSynToRemove);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::degenerate_uniform(number p, const char* subset)
/**
 * Removes a percentage p of synapses from subset subsetIndex
 */
void SynapseDistributor::degenerate_uniform(number p, const char* subset)
{
	if(p < 0 || p > 1)
		UG_THROW("SynapseDistributor::degenerate_uniform(number p, const char* subset): Specified percentage incorrect.");

//	Get subset index from subset name
	int si = pm_SubsetHandler->get_subset_index(subset);

	degenerate_uniform(p, si);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::get_subset_length(int si)
/**
 * Calculate and return length of specified subset in [m]
 */
number SynapseDistributor::get_subset_length(int si)
{
//	Save subset coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_SubsetHandler->begin<Edge>(si, 0),
										 pm_SubsetHandler->end<Edge>(si, 0));

	number totLength = 0.0;

	for(size_t i = 0; i < vEdges.size(); ++i)
	{
		Edge* e = vEdges[i];
		number edgeLength = EdgeLength(e, m_aaPosition);

		totLength += edgeLength;
	}

	return totLength;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::get_subset_length(const char* subset)
/**
 * Calculate and return length of specified subset in [m]
 */
number SynapseDistributor::get_subset_length(const char* subset)
{
//	Get subset index from subset name
	int si = pm_SubsetHandler->get_subset_index(subset);

	return get_subset_length(si);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::num_synapses(vector<Edge*> vEdges, bool bActive)
/**
 * Wrapper for num_synapses: Returns number of synapses in the given vector of edges
 */
size_t SynapseDistributor::num_synapses(vector<Edge*> vEdges, bool bActive, number time)
{
//	Determine number of synapses in the given vector of edges
	size_t numSynapses = 0;
	for(size_t i = 0; i < vEdges.size(); ++i)
	{
	//	Decide, if we want to count all or only active synapses
		if(!bActive)
			numSynapses += m_aaSynInfo[vEdges[i]].size();
		else
		{
			for(size_t j = 0; j < m_aaSynInfo[vEdges[i]].size(); ++j)
			{
				if(m_aaSynInfo[vEdges[i]][j].m_onset <= time && m_aaSynInfo[vEdges[i]][j].m_onset + m_aaSynInfo[vEdges[i]][j].m_tau*6.0 >= time)
					numSynapses += 1;
				else
					continue;
			}
		}
	}

	return numSynapses;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::num_synapses(int si)
/**
 * Returns the number of synapses in a given subset.
 */
size_t SynapseDistributor::num_synapses(int si)
{
//	Save subset coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_SubsetHandler->begin<Edge>(si, 0), pm_SubsetHandler->end<Edge>(si, 0));

	return num_synapses(vEdges, false, 0.0);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::num_synapses(const char* subset)
/**
 * Returns the number of synapses in a given subset.
 */
size_t SynapseDistributor::num_synapses(const char* subset)
{
//	Get subset index from subset name
	int si = pm_SubsetHandler->get_subset_index(subset);

	return num_synapses(si);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::num_synapses()
/**
 * Returns the number of synapses in the whole grid.
 */
size_t SynapseDistributor::num_synapses()
{
//	Save subset coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_Grid->begin<Edge>(0), pm_Grid->end<Edge>(0));

	return num_synapses(vEdges, false, 0.0);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::num_active_synapses(number time)
/**
 * Returns the number of active synapses in the whole grid.
 */
size_t SynapseDistributor::num_active_synapses(number time)
{
//	Save subset coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_Grid->begin<Edge>(0), pm_Grid->end<Edge>(0));

	return num_synapses(vEdges, true, time);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::num_active_synapses(number time)
/**
 * Returns the number of active synapses in the given subset.
 */
size_t SynapseDistributor::num_active_synapses(number time, int si)
{
//	Save subset coarse grid edges in a vector
	vector<Edge*> vEdges = vector<Edge*>(pm_Grid->begin<Edge>(0), pm_Grid->end<Edge>(0));

	return num_synapses(vEdges, true, time);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::activity_info()
/**
 * Prints start and end time for each synapse
 */
void SynapseDistributor::activity_info()
{
	size_t counter = 0;

	for(EdgeIterator eIter = pm_Grid->begin<Edge>(0); eIter != pm_Grid->end<Edge>(0); ++eIter)
	{
		Edge* e = *eIter;

		for(size_t i = 0; i < m_aaSynInfo[e].size(); ++i)
		{
			UG_LOG("Start: " << m_aaSynInfo[e][i].m_onset << " ; End: " << m_aaSynInfo[e][i].m_onset + m_aaSynInfo[e][i].m_tau*6.0 << endl);
			counter += 1;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
//	vector<Synapse*> SynapseDistributor::find_synapses(number x, number y, number z) const
/**
 * Returns vector with pointers to synapses at this position, if existent
 */
vector<synapse_handler::SynapseInfo*> SynapseDistributor::find_synapses(number x, number y, number z) const
{
///	TODO: Implement me!

	vector<synapse_handler::SynapseInfo*> vSyn;
	UG_THROW("vector<Synapse*> SynapseDistributor::find_synapses(number x, number y, number z) const. NOT IMPLEMENTED!");

	return vSyn;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::synapses_on_edge(const Edge* e, const size_t corner, const number time, number& current) const
/**
 * Returns whether any active synapse belongs to the given edge and corner
 */
bool SynapseDistributor::has_active_synapses(const Edge* e, const size_t corner, const number time, number& current) const
{
//	Check, if edge has synapses
	if(m_aaSynInfo[e].size() == 0)
		return false;

//	Determine number of active synapses on e
	size_t numActiveSynapses = 0;

	for(size_t i = 0; i < m_aaSynInfo[e].size(); ++i)
	{
		const synapse_handler::SynapseInfo* syn = &m_aaSynInfo[e][i];

	//	Check, if synapse is active at that time
		if((syn->m_onset > time) || (syn->m_onset + syn->m_tau*6.0 < time))
			continue;

	//	Check, if synapse belongs to given corner
		if ((corner==0 && syn->m_locCoords >= 0.5) || (corner==1 && syn->m_locCoords < 0.5))
			continue;

		numActiveSynapses += 1;
	}

//	Set electrode current times number of active synapses on e
	if(numActiveSynapses == 0)
		return false;
	else
	{
		current += -2e-14;
		current *= numActiveSynapses;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::print_status()
/**
 * Prints the status of the given grid, subsets and number of vertices in those subsets.
 * For internal testing purposes.
 */
void SynapseDistributor::print_status()
{
	//Iterates over all subsets and prints the number of edges in each subset
	cout<<"\n#####################################Status#####################################\n";
	for(int i=0; i<pm_SubsetHandler->num_subsets(); ++i) {
		cout<<"Subset '"<<pm_SubsetHandler->subset_info(i).name<<"("<<i<<")':\t\tEdges: "
			<<pm_SubsetHandler->num<Edge>(i)<<"\tSynapses: "<<num_synapses(i)<<endl;
	}
	cout<<endl<<"Total:\t\t\tEdges: "<<pm_Grid->num_edges()<<"\tSynapses: "<<num_synapses()<<endl
		<<"################################################################################"<<endl;

}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::get_last_message()
/**
 * Returns the last message generated.
 * For internal testing purposes.
 */
string SynapseDistributor::get_last_message()
{
	return m_LastMessage;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::set_output_filename()
/**
 * Sets the output filename in member m_OutputFile
 */
void SynapseDistributor::set_outfile(std::string outfile)
{
	m_OutputFile = outfile;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::save_changes()
/**
 * Saves all changes made to the input geometry file on disk.
 */
bool SynapseDistributor::export_grid()
{

	m_LastMessage = "export_grid(): Saving to "+m_OutputFile+"...";

	if(SaveGridToFile(*pm_Grid, *pm_SubsetHandler, m_OutputFile.c_str())) {
		m_LastMessage.append("done.");
		return true;
	}
	else {
		m_LastMessage.append("SaveGridToFile(*pm_Grid, *pm_SubsetHandler, m_outputFile.c_str()) failed.");
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::export_grid(string outfile)
/**
 * Saves all changes made to the input geometry file by exporting to specified outfile.
 */
bool SynapseDistributor::export_grid(std::string outfile)
{

	m_LastMessage = "export_grid(): Saving to "+outfile+"...";

	if(SaveGridToFile(*pm_Grid, *pm_SubsetHandler, outfile.c_str())) {
		m_LastMessage.append("done.");
		return true;
	}
	else {
		m_LastMessage.append("SaveGridToFile(*pm_Grid, *pm_SubsetHandler, outfile.c_str()) failed.");
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::get_grid()
/**
 * Returns pointer to current grid.
 */
MultiGrid* SynapseDistributor::get_grid()
{
	return pm_Grid;
}
////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::get_subset_handler()
/**
 * Returns pointer to current grid.
 */
MGSubsetHandler* SynapseDistributor::get_subset_handler()
{
	return pm_SubsetHandler;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::get_dummy_current()
/**
 * Returns dummy electrode current
 */
number SynapseDistributor::get_dummy_current() const { return -2e-14; /*2e-14C/ms=20pA*/}


} // namespace cable_neuron
} // namespace ug




