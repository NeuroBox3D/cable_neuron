/*
 * split_synapse_distributor.h
 *
 *  Created on: Mar 24, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_DISTRIBUTOR_SPLIT_SYNAPSE_DISTRIBUTOR_H_
#define SPLIT_SYNAPSE_DISTRIBUTOR_SPLIT_SYNAPSE_DISTRIBUTOR_H_

// system includes
#include <vector>
#include <string>

// ugcore includes
#include "common/util/smart_pointer.h"
#include "lib_grid/lib_grid.h"
#include "common/types.h"
#include "lib_grid/global_attachments.h"

//SplitSynapses
#include "../split_synapse_handler/base_synapse.h"
#include "../split_synapse_handler/pre_synapse.h"
#include "../split_synapse_handler/post_synapse.h"
#include "../split_synapse_handler/synapse_dealer.h"
#include "../split_synapse_handler/split_synapse_info_io_traits.h"
#include "../split_synapse_handler/alpha_post_synapse.h"

// boost includes
#include "boost/lexical_cast.hpp"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class SplitSynapseDistributor {

private:
	typedef Attachment<std::vector<IBaseSynapse*> > AVSSynapse;
	bool m_bDomBased;

	std::string m_InputFile;
	std::string m_OutputFile;
	std::string m_LastMessage;

	/// balls / stimulation regions
	std::vector<std::pair<vector3, number> >m_balls;

//	Pointers to heap memory, have to be freed in destructor
	MultiGrid* pm_Grid;
	MGSubsetHandler* pm_SubsetHandler;
	AVSSynapse m_aSSyn;
	Grid::EdgeAttachmentAccessor<AVSSynapse> m_aaSSyn;
	Grid::VertexAttachmentAccessor<APosition> m_aaPosition;

public:
	//	Constructor
	SplitSynapseDistributor(std::string infile, std::string outfile, bool bRemoveExistingSynapses);
	SplitSynapseDistributor(SmartPtr<Domain3d> dom, std::string outfile, bool bRemoveExistingSynapses);
	virtual ~SplitSynapseDistributor();

	//	Multigrid Synapse transer methods
	void CopySynapsesToAllLevels();
	void CopySynapsesFromParentToChild(Edge* parent);

	//	methods, general functionality
	void remove_synapse(Edge* e);
	void remove_all_synapses(Edge* e);
	void clear();
	void clear(int subsetIndex);

	void place_synapse(Edge* e, IBaseSynapse* s);
	void place_synapses_uniform(std::vector<Edge*> vEdges, std::vector<IBaseSynapse*> s);

	void place_synapses_uniform(std::vector<IBaseSynapse*> s);
	void place_synapses_uniform(int si, std::vector<IBaseSynapse*> s);
	void place_synapses_uniform(std::vector<number> s);
	void place_synapses_uniform(int si, std::vector<number> s);

	void degenerate_uniform(std::vector<Edge*> vEdges, size_t numSynapses);
	void degenerate_uniform(number p);
	void degenerate_uniform(number p, int si);
	void degenerate_uniform(number p, const char* subset);

	number get_subset_length(int si);
	number get_subset_length(const char* subset);
	//todo:
	size_t num_synapses(std::vector<Edge*> vEdges, bool bActive, number time);

	size_t num_synapses(int si);
	size_t num_synapses(const char* subset);
	size_t num_synapses();
	size_t num_active_synapses(number time);
	size_t num_active_synapses(number time, int si);

//	void activity_info();

//	bool has_active_synapses(const Edge* e, const size_t corner, const number time, number& current) const;

	//	testing and helper
	void print_status();
	std::string get_last_message();
	void set_outfile(std::string outfile);
	bool export_grid();
	bool export_grid(std::string outfile);
	MultiGrid* get_grid();
	MGSubsetHandler* get_subset_handler();

	number get_dummy_current() const;


};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_DISTRIBUTOR_SPLIT_SYNAPSE_DISTRIBUTOR_H_ */
