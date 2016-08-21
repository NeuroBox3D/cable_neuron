/*
 * SplitSynapseHandler.h
 *
 *  Created on: Mar 23, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_
#define SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_

// includes
#include <map>
#include <string>
#include <algorithm>

// boost includes for random numbers
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

// ug
#include "common/log.h"
#include "lib_grid/lib_grid.h"
#include "lib_grid/global_attachments.h"
#include "lib_grid/tools/copy_attachment_handler.h"
#include "synapse_dealer.h"
#include "split_synapse_info_io_traits.h"
#include "../cable_disc/cable_equation.h"
#include "../split_synapse_distributor/split_synapse_distributor.h"
#include "split_synapse_attachment_handler.h"
#include "base_synapse.h"
#include "pre_synapse.h"
#include "post_synapse.h"
#include "pcl/pcl_process_communicator.h"

namespace ug {
namespace cable_neuron {

template <typename TDomain>
class CableEquation;
template <typename TDomain>
class ICableMembraneTransport;

namespace synapse_handler {

template <typename TDomain>
class SplitSynapseHandler
{
private:
	typedef Attachment<std::vector<IBaseSynapse*> > AVSSynapse;

	AVSSynapse m_aSSyn;
	Grid::EdgeAttachmentAccessor<AVSSynapse> m_aaSSyn;
	bool m_bInited;
	SmartPtr<MultiGrid> m_spGrid;
	SmartPtr<CableEquation<TDomain> > m_spCEDisc;
	SplitSynapseAttachmentHandler m_ssah;
	SmartPtr<ApproximationSpace<TDomain> > m_spApprox;


	/*	//probably unused
	std::vector<IBaseSynapse*> m_vAllSynapses;
	std::vector<IPreSynapse*> m_vPreSynapses;
	std::vector<IPostSynapse*> m_vPostSynapses;*/	// map with key ID instead of vector!

	//use these maps
	//std::map<SYNAPSE_ID, IBaseSynapse*> m_mAllSynapses; //contains all local synapses on grid by their respective id->IBaseSynapse pointer
	std::map<SYNAPSE_ID, IPreSynapse*> m_mPreSynapses; // " ... "
	std::map<SYNAPSE_ID, IPostSynapse*> m_mPostSynapses; // " ... "

	std::map<SYNAPSE_ID, IPreSynapse*> m_mActivePreSynapses; //internal memory of currently active presynapses
	std::map<SYNAPSE_ID, Edge*> m_mPreSynapseIdToEdge; //maps presynapses to their edge


public:
	SplitSynapseHandler();

	virtual ~SplitSynapseHandler() {}

	void set_ce_object(SmartPtr<CableEquation<TDomain> > disc) {
		UG_COND_THROW(m_bInited, "The CableEquation object associated to this synapse handler "
					  "must not be changed\nafter addition of the original CableEquation object "
					  "to the domain discretization.");

		// set cable equation disc object
		m_spCEDisc = disc;
	}

	number current_on_edge(const Edge* e, size_t scv, number t);

	void all_synapses();

	void update_presyn(number time);

	//interface for cable_equation
	bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current);
	void grid_first_available();
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_ */
