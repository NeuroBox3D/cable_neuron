/*
 * PreAlphaSynapse.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#include "AlphaPreSynapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

AlphaPreSynapse::AlphaPreSynapse(
		const number& location,
		const number& onset)

:IPreSynapse(0, 0, location),
 m_onset(onset)
{
}

AlphaPreSynapse::AlphaPreSynapse(
		const unsigned long id,
		const unsigned long postsynapse_id,
		const number& location,
		const number& onset)

:IPreSynapse(id, postsynapse_id, location),
 m_onset(onset)
{
}

AlphaPreSynapse::~AlphaPreSynapse()
{
}

void AlphaPreSynapse::update(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy update atm
}


bool AlphaPreSynapse::is_active(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy
	return (t >= m_onset);
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
