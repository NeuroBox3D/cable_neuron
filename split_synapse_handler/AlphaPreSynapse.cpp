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

PreAlphaSynapse::PreAlphaSynapse(
		const number& location,
		const number& onset)

:IPreSynapse(0, 0, location),
 m_onset(onset)
{
}

PreAlphaSynapse::PreAlphaSynapse(
		const unsigned long id,
		const unsigned long postsynapse_id,
		const number& location,
		const number& onset)

:IPreSynapse(id, postsynapse_id, location),
 m_onset(onset)
{
}

PreAlphaSynapse::~PreAlphaSynapse()
{
}

void PreAlphaSynapse::update(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy update atm
}


bool PreAlphaSynapse::is_active(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy
	return (t >= m_onset);
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
