/*
 * synapse_distributor.h
 *
 *  Created on: 03.11.2014
 * 	Author: Lukas Reinhardt
 */


#ifndef __SynapseDistributor_h__
#define __SynapseDistributor_h__

/* system includes */
#include <stddef.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "lib_grid/lib_grid.h"
#include "lib_grid/global_attachments.h"

#include "../synapse_handler/grid/synapse_info_io_traits.h"

/// boost
#include "boost/lexical_cast.hpp"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>




using namespace std;


namespace ug {


////////////////////////////////////////////////////////////////////////////////////////////
///	SynapseDistributor
/// @todo: Make SynapseDistributor independent of dimension
///		   (using dim template and MathVector<dim> for coordinates)!
////////////////////////////////////////////////////////////////////////////////////////////
class SynapseDistributor
{
	typedef Attachment<vector<SynapseInfo> > AVSynapse;

public:
//	Constructor
	SynapseDistributor(string infile, string outfile, bool bRemoveExistingSynapses);
	SynapseDistributor(SmartPtr<Domain3d> dom, string outfile, bool bRemoveExistingSynapses);

//	Standard destructor
	~SynapseDistributor();

//	Multigrid Synapse transer methods
	void CopySynapsesToAllLevels();
	void CopySynapsesFromParentToChild(Edge* parent);

//	methods, general functionality
	void remove_synapse(Edge* e);
	void remove_all_synapses(Edge* e);
	void clear();
	void clear(int subsetIndex);

	void set_activation_timing(double start_time, double duration, double start_time_dev, double duration_dev);

	void place_synapse(Edge* e);
	void place_synapses_uniform(vector<Edge*> vEdges, size_t numSynapses);
	void place_synapses_uniform(size_t numSynapses);
	void place_synapses_uniform(int si, size_t numSynapses);
	void place_synapses(vector<double> distr, size_t numSynapses);

	void degenerate_uniform(vector<Edge*> vEdges, size_t numSynapses);
	void degenerate_uniform(double p);
	void degenerate_uniform(double p, int si);

	size_t num_synapses(vector<Edge*> vEdges, bool bActive, number time);
	size_t num_synapses(int si);
	size_t num_synapses();
	size_t num_active_synapses(number time);
	size_t num_active_synapses(number time, int si);

	void activity_info();

	vector<SynapseInfo*> find_synapses(double x, double y, double z) const;

	template <size_t dim>
	bool has_active_synapses(const MathVector<dim>& c, const number time, number& current) const;
	bool has_active_synapses(const Edge* e, const size_t corner, const number time, number& current) const;

//	testing and helper
	void print_status();
	string get_last_message();
	bool export_grid();
	MultiGrid* get_grid();
	MGSubsetHandler* get_subset_handler();

	/*!
	 * \brief dummy current
	 * \remarks to be exchanged by real implementation asap
	 */
	number get_dummy_current() const;


private:

	bool m_bDomBased;

	string m_InputFile;
	string m_OutputFile;
	string m_LastMessage;

//	Pointers to heap memory, have to be freed in destructor
	MultiGrid* pm_Grid;
	MGSubsetHandler* pm_SubsetHandler;
	AVSynapse m_aSynInfo;
	Grid::EdgeAttachmentAccessor<AVSynapse> m_aaSynInfo;
	Grid::VertexAttachmentAccessor<APosition> m_aaPosition;
};

}

#include "synapse_distributor_impl.h"

#endif

