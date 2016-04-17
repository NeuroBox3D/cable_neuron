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

namespace ug {
namespace cable_neuron {
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

	std::vector<IBaseSynapse*> m_vSynapses;

public:
	SplitSynapseHandler();
	virtual ~SplitSynapseHandler() {}

	void set_ce_object(SmartPtr<CableEquation<TDomain> > disc);

	number current_on_edge(const Edge* e, size_t scv, number t);

//	template <typename TSyn>
//	std::vector<TSyn*>::iterator begin()
//	{
//		std::vector<IBaseSynapse*> vSyn = all_synapses();
//		std::vector<TSyn*> vTSyn;
//		TSyn* s = new TSyn;
//
//		for(size_t i = 0; i<vSyn.size(); ++i) {
//			if(vSyn[i]->type() == s->type()) {
//				vTSyn.push_back(vSyn[i]);
//			}
//		}
//
//		delete s;
//		return vTSyn.begin();
//	}
//
//	template <typename TSyn>
//	std::vector<TSyn*>::iterator end()
//	{
//		std::vector<IBaseSynapse*> vSyn = all_synapses();
//		std::vector<TSyn*> vTSyn;
//		TSyn* s = new TSyn;
//
//		for(size_t i = 0; i<vSyn.size(); ++i) {
//			if(vSyn[i]->type() == s->type()) {
//				vTSyn.push_back(vSyn[i]);
//			}
//		}
//
//		delete s;
//		return vTSyn.end();
//	}

	std::vector<IBaseSynapse*> all_synapses();

	//interface for cable_equation
	bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current);
	void grid_first_available();
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_ */
