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

	std::vector<IBaseSynapse*> m_vAllSynapses;
	std::vector<IPreSynapse*> m_vPreSynapses;
	std::vector<IPostSynapse*> m_vPostSynapses;

	std::map<unsigned long long,IPostSynapse*> m_mPostSynapses;
	std::map<unsigned long long,IPreSynapse*> m_mActivePreSynapses;

public:
	SplitSynapseHandler();
	virtual ~SplitSynapseHandler() {}

	void set_ce_object(SmartPtr<CableEquation<TDomain> > disc);

	number current_on_edge(const Edge* e, size_t scv, number t);

	std::vector<IBaseSynapse*> all_synapses();
	std::vector<IPreSynapse*> all_pre_synapses();
	std::vector<IPostSynapse*> all_post_synapses();

	void update_presyn(number time);

	//interface for cable_equation
	bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current);
	void grid_first_available();
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_ */
