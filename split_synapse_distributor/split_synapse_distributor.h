/*
 * split_synapse_distributor.h
 *
 *  Created on: Mar 24, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_DISTRIBUTOR_SPLIT_SYNAPSE_DISTRIBUTOR_H_
#define SPLIT_SYNAPSE_DISTRIBUTOR_SPLIT_SYNAPSE_DISTRIBUTOR_H_

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class SplitSynapseDistributor {
public:
	SplitSynapseDistributor();
	virtual ~SplitSynapseDistributor();
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_DISTRIBUTOR_SPLIT_SYNAPSE_DISTRIBUTOR_H_ */
