/*
 * synapse_container.h
 *
 *  Created on: Apr 17, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SYNAPSE_CONTAINER_H_
#define SPLIT_SYNAPSE_HANDLER_SYNAPSE_CONTAINER_H_

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class ISynapseContainer {
public:
	virtual ~ISynapseContainer() {}
	virtual std::vector<IBaseSynapse*> get_synapses() = 0;
	virtual size_t size() = 0;
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SYNAPSE_CONTAINER_H_ */
