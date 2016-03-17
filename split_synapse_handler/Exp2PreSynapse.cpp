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
Exp2PreSynapse::Exp2PreSynapse(
		const number& location,
		const number& onset)

:IPreSynapse(0, 0, location),
 m_onset(onset)
{
}

Exp2PreSynapse::Exp2PreSynapse(
		const unsigned long id,
		const unsigned long postsynapse_id,
		const number& location,
		const number& onset)

:IPreSynapse(id, postsynapse_id, location),
 m_onset(onset)
{
}

Exp2PreSynapse::~Exp2PreSynapse()
{
}

void Exp2PreSynapse::update(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy update atm
}

bool Exp2PreSynapse::is_active(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy
	return (t >= m_onset);
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
