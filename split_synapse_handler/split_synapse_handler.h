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
#include "synapse_attachment_serializer.h"
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

/**
 * Provides iterator-like access to synapses of type TSyn.
 * Use SplitSynapseHandler::begin<TSyn> and SplitSynapseHandler::end<TSyn>
 * to get iterators
 */
template <typename TSyn>
class SynapseIter
{
	std::vector<IBaseSynapse*>::iterator m_it;
public:
	SynapseIter(std::vector<IBaseSynapse*>::iterator it) :m_it(it) {}

	TSyn* operator*() {
		return dynamic_cast<TSyn*>(*m_it);
	}

	TSyn* operator++() {
		return dynamic_cast<TSyn*>(*(++m_it) );
	}

	bool operator!=(SynapseIter* rhs) {
		return m_it != rhs->m_it;
	}
};

/**
 * SplitSynapseHandler is the central management class for everything related to SplitSynapses.
 * Provides iterators for each registered synapse type.
 */
template <typename TDomain>
class SplitSynapseHandler
{
private:
	typedef Attachment<std::vector<IBaseSynapse*> > AVSSynapse;
	Grid::VertexAttachmentAccessor<APosition> m_aaPosition;


	AVSSynapse m_aSSyn;
	Grid::EdgeAttachmentAccessor<AVSSynapse> m_aaSSyn;

	bool m_bInited;
	SmartPtr<MultiGrid> m_spGrid;
	SmartPtr<CableEquation<TDomain> > m_spCEDisc;
	SplitSynapseAttachmentHandler m_ssah;
	SmartPtr<ApproximationSpace<TDomain> > m_spApprox;

	std::vector<IBaseSynapse*> m_vAllSynapses;

	std::map<SYNAPSE_ID, IPreSynapse*> m_mPreSynapses; // " ... "
	std::map<SYNAPSE_ID, IPostSynapse*> m_mPostSynapses; // " ... "

	std::map<SYNAPSE_ID, IPreSynapse*> m_mActivePreSynapses; //internal memory of currently active presynapses
	std::map<SYNAPSE_ID, Edge*> m_mPreSynapseIdToEdge; //maps presynapses to their edge

	SmartPtr<GeomObjDataSerializer<Edge> > m_spSAS;
	MessageHub::SPCallbackId m_spGridDistributionCallbackID;

	//TODO: don't use copy ctor at the moment
	SplitSynapseHandler(const SplitSynapseHandler& sh);

	/**
	 * Fills the member m_vAllSynapses with every synapse currently on the grid.
	 * Called by grid_first_available.
	 */
	void collect_synapses_from_grid();
	/**
	 * used for sorting
	 */
	struct __comp{
		bool operator() (IBaseSynapse* a, IBaseSynapse* b) {return a->type() < b->type();}
	};

public:
	/**
	 * Retrieves approximation space and grid when available.
	 * Also gathers every synapse in m_vAllSynapses for quick access/iterators
	 */
	void grid_first_available();

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
	 * todo: don't use
	 */
	std::vector<IBaseSynapse*>& get_synapses() {return m_vAllSynapses;}


	std::vector<IBaseSynapse*>& get_synapses_on_edge(Edge* e) {return m_aaSSyn[e];}


	/**
	 * update on synapses becoming active/inactive in timestep time
	 */
	void update_presyn(number time);


	/**
	 * Dummy interface method for compatibility reasons, check current on edge for functionality
	 */
	bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current);

	/**
	 * Returns iterator to the desired synapse type.
	 * Templates have to be instantiated for use in LUA.
	 */
	template <typename TSyn>
	SmartPtr<SynapseIter<TSyn> > begin() {
		std::vector<IBaseSynapse*>::iterator it = m_vAllSynapses.begin();
		for(; it != m_vAllSynapses.end(); ++it ) {
			if(dynamic_cast<TSyn*>(*it) ) {//found the first TSyn*
				 break;
			}
		}
		return SmartPtr<SynapseIter<TSyn> >(new SynapseIter<TSyn>(it));
	}

	/**
	 * Returns an end iterator to the desired synapse type.
	 * Templates have to be instantiated for use in LUA.
	 */
	template <typename TSyn>
	SmartPtr<SynapseIter<TSyn> > end() {
		std::vector<IBaseSynapse*>::iterator it = m_vAllSynapses.begin();
		for(; it != m_vAllSynapses.end(); ++it ) {
			if(dynamic_cast<TSyn*>(*it) ) {//found the first TSyn*
				//std::cout << *it << std::endl;
				break;
			}
		}

		while(it != m_vAllSynapses.end()) {
			if(!dynamic_cast<TSyn*>(*it) ) {//found last TSyn*
				break;
			}
			++it;
		}
		return SmartPtr<SynapseIter<TSyn> >(new SynapseIter<TSyn>(it));
	}

	/**
	 * Prints list with all synapses and their parameters/id's etc.
	 */
	void show_status() {
		for(size_t i=0; i<m_vAllSynapses.size(); ++i) {
			std::cout << i << ": " << " " << m_vAllSynapses[i] << std::endl;
		}
		std::cout << std::endl;
		//std::cout << "Synapse iterator begin(): " << &*(m_vAllSynapses.begin()) << std::endl;
		//std::cout << "Synapse iterator end(): " << &*(m_vAllSynapses.end()) << std::endl;
		//std::cout << "AlphaPreSynapse iterator end(): " << &*end<AlphaPreSynapse>() << std::endl;
	}


protected:
	void grid_distribution_callback(const GridMessage_Distribution& gmd);

};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_ */
