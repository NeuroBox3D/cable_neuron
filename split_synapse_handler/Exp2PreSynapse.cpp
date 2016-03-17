/*
 * PreExp2Synapse.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#include "Exp2PreSynapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

/**
 * Synapse with id=0
 */
PreExp2Synapse::PreExp2Synapse(
		const number& location,
		const number& onset)

:IPreSynapse(0, 0, location),
 m_onset(onset)
{
}

PreExp2Synapse::PreExp2Synapse(
		const unsigned long id,
		const unsigned long postsynapse_id,
		const number& location,
		const number& onset)

:IPreSynapse(id, postsynapse_id, location),
 m_onset(onset)
{
}

PreExp2Synapse::~PreExp2Synapse()
{
}

void PreExp2Synapse::update(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy update atm
}

bool PreExp2Synapse::is_active(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy
	return (t >= m_onset);
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
