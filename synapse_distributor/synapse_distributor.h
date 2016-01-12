/*
 * synapse_distributor.h
 *
 *  Created on: 03.11.2014
 * 	Author: Lukas Reinhardt
 */


#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_DISTRIBUTOR__SYNAPSE_DISTRIBUTOR_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_DISTRIBUTOR__SYNAPSE_DISTRIBUTOR_H__

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

#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "common/util/smart_pointer.h"

#include "lib_grid/lib_grid.h"
#include "lib_grid/global_attachments.h"

#include "../synapse_handler/grid/synapse_info_io_traits.h"

/// boost
#include "boost/lexical_cast.hpp"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


namespace ug {
namespace cable_neuron {


////////////////////////////////////////////////////////////////////////////////////////////
///	SynapseDistributor
/// @todo: Make SynapseDistributor independent of dimension
///		   (using dim template and MathVector<dim> for coordinates)!
////////////////////////////////////////////////////////////////////////////////////////////
class SynapseDistributor
{
	typedef Attachment<std::vector<SynapseInfo> > AVSynapse;

public:
//	Constructor
	SynapseDistributor(std::string infile, std::string outfile, bool bRemoveExistingSynapses);
	SynapseDistributor(SmartPtr<Domain3d> dom, std::string outfile, bool bRemoveExistingSynapses);

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

	void set_activation_timing(number start_time, number duration, number start_time_dev, number duration_dev);

	void place_synapse(Edge* e);
	void place_synapses_uniform(std::vector<Edge*> vEdges, size_t numSynapses);
	void place_synapses_uniform(size_t numSynapses);
	void place_synapses_uniform(int si, size_t numSynapses);
	void place_synapses_uniform(const char*, number density);
	void place_synapses(std::vector<number> distr, size_t numSynapses);

	void degenerate_uniform(std::vector<Edge*> vEdges, size_t numSynapses);
	void degenerate_uniform(number p);
	void degenerate_uniform(number p, int si);

	number get_subset_length(int si);
	size_t num_synapses(std::vector<Edge*> vEdges, bool bActive, number time);
	size_t num_synapses(int si);
	size_t num_synapses();
	size_t num_active_synapses(number time);
	size_t num_active_synapses(number time, int si);

	void activity_info();

	std::vector<SynapseInfo*> find_synapses(number x, number y, number z) const;

	template <size_t dim>
	bool has_active_synapses(const MathVector<dim>& c, const number time, number& current) const;
	bool has_active_synapses(const Edge* e, const size_t corner, const number time, number& current) const;

//	testing and helper
	void print_status();
	std::string get_last_message();
	void set_outfile(std::string outfile);
	bool export_grid();
	bool export_grid(std::string outfile);
	MultiGrid* get_grid();
	MGSubsetHandler* get_subset_handler();

	/*!
	 * \brief dummy current
	 * \remarks to be exchanged by real implementation asap
	 */
	number get_dummy_current() const;


private:

	bool m_bDomBased;

	std::string m_InputFile;
	std::string m_OutputFile;
	std::string m_LastMessage;

//	Pointers to heap memory, have to be freed in destructor
	MultiGrid* pm_Grid;
	MGSubsetHandler* pm_SubsetHandler;
	AVSynapse m_aSynInfo;
	Grid::EdgeAttachmentAccessor<AVSynapse> m_aaSynInfo;
	Grid::VertexAttachmentAccessor<APosition> m_aaPosition;
};


} // namespace cable_neuron
} // namespace ug

#include "synapse_distributor_impl.h"

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_DISTRIBUTOR__SYNAPSE_DISTRIBUTOR_H__

