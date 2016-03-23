/*
 * SplitSynapseHandler.h
 *
 *  Created on: Mar 23, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_
#define SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_

#include "split_synapse_info_io_traits.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {


class SplitSynapseHandler {
public:
	SplitSynapseHandler();
	virtual ~SplitSynapseHandler();
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_HANDLER_H_ */
