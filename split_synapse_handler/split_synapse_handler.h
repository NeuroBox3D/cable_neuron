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


#include "split_synapse_info_io_traits.h"
#include "../cable_disc/cable_equation.h"
#include "../split_synapse_distributor/split_synapse_distributor.h"
#include "split_synapse_attachment_handler.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {


class SplitSynapseHandler
{
private:
	typedef Attachment<std::vector<IBaseSynapse*> > AVSSynapse;


	AVSSynapse m_aSSyn;
	Grid::EdgeAttachmentAccessor<AVSSynapse> m_aaSSyn;
	bool m_bInited;
	SmartPtr<MultiGrid> m_spGrid;
	SplitSynapseAttachmentHandler m_ssah;

public:
	SplitSynapseHandler();
	virtual ~SplitSynapseHandler() {}

	number current_on_edge(Edge* e, number t);

	//interface for cable_equation
	bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current);
	void grid_first_available();
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_ */
