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
#include "alpha_pre_synapse.h"
#include "alpha_post_synapse.h"
#include "exp2_post_synapse.h"
#include "exp2_pre_synapse.h"
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
	std::vector<IBaseSynapse*> m_vAllSynapses;

	/*	//probably unused
	std::vector<IPreSynapse*> m_vPreSynapses;
	std::vector<IPostSynapse*> m_vPostSynapses;*/	// map with key ID instead of vector!

	//use these maps
	//std::map<SYNAPSE_ID, IBaseSynapse*> m_mAllSynapses; //contains all local synapses on grid by their respective id->IBaseSynapse pointer
	std::map<SYNAPSE_ID, IPreSynapse*> m_mPreSynapses; // " ... "
	std::map<SYNAPSE_ID, IPostSynapse*> m_mPostSynapses; // " ... "

	std::map<SYNAPSE_ID, IPreSynapse*> m_mActivePreSynapses; //internal memory of currently active presynapses
	std::map<SYNAPSE_ID, Edge*> m_mPreSynapseIdToEdge; //maps presynapses to their edge

	std::string m_iterator_type;
	std::vector<IBaseSynapse*>::iterator m_split_synapse_it;

	//TODO: don't use copy ctor at the moment
	SplitSynapseHandler(const SplitSynapseHandler& sh);


public:
	/**
	 * Declares global attachment AVSSynapse and attaches to split synapse attachment handler.
	 */
	SplitSynapseHandler();

	/**
	 * Default dtor
	 */
	virtual ~SplitSynapseHandler() {}

	/**
	 * Set cable equation disc object
	 */
	void set_ce_object(SmartPtr<CableEquation<TDomain> > disc) {
		UG_COND_THROW(m_bInited, "The CableEquation object associated to this synapse handler "
					  "must not be changed\nafter addition of the original CableEquation object "
					  "to the domain discretization.");
		m_spCEDisc = disc;
	}

	/**
	 * Returns total synaptic current on Edge e
	 */
	number current_on_edge(const Edge* e, size_t scv, number t);

	/**
	 * Fills the member m_vAllSynapses with every synapse currently on the grid.
	 * Called by grid_first_available.
	 */
	void all_synapses();


	/**
	 * update on synapses becoming active/inactive in timestep time
	 */
	void update_presyn(number time);


	/**
	 * Dummy interface method for compatibility reasons, check current on edge for functionality
	 */
	bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current);

	/**
	 * Retrieves approximation space and grid when available.
	 * Also gathers every synapse in m_vAllSynapses for quick access/iterators
	 */
	void grid_first_available();






	/**
	 * Iterator related methods
	 */


	/**
	 * Resets the iterator and sets the type of to be selected synapses to tsyn
	 */
	void reset_iterator(std::string tsyn) {
		m_split_synapse_it = m_vAllSynapses.begin();
		m_iterator_type = tsyn;
	}

	/**
	 * Points to the next synapse until processed (with set_onset etc...) of the specified type.
	 * Returns true in case of a valid next synapse.
	 * Returns false if the iterator reaches the end
	 */
	bool next() {
		while(m_split_synapse_it != m_vAllSynapses.end()) {
			IBaseSynapse* s = *m_split_synapse_it;
			if( s->name() == m_iterator_type)
				return true;
			++m_split_synapse_it;
		}
		return false;
	}

	/**
	 * Set onset of the following synapse types:
	 *
	 * -ALPHA_PRE_SYNAPSE
	 * -EXP2_PRE_SYNAPSE
	 *
	 * Advances the iterator by one BaseSynapse obejct.
	 */
	void set_onset(const number& t) {
		IBaseSynapse* s = *m_split_synapse_it;
		if(m_iterator_type == "ALPHA_PRE_SYNAPSE") {
			( (AlphaPreSynapse*) s)->set_onset(t);
		} else if(m_iterator_type == "EXP2_PRE_SYNAPSE") {
			( (Exp2PreSynapse*) s)->set_onset(t);
		}
		++m_split_synapse_it;
	}




};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_ */
