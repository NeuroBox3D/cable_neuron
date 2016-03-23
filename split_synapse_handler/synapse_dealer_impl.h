/*
 * SynapseDealer_impl.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

namespace ug {
namespace cable_neuron {
namespace synapse_handler {


template <typename TSyn>
void SynapseDealer::register_synapse_type(std::string t)
{
	// TODO: check if name already exists and throw if that is the case
	ISynapseGenerator* sg = new SynapseGenerator<TSyn>();
	m_register[t] = sg;
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
