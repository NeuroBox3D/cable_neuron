/*
 * PreAlphaSynapse.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#include "PreAlphaSynapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

PreAlphaSynapse::PreAlphaSynapse(
		const number& location,
		const number& onset)

:IPreSynapse(0, location),
 m_onset(onset)
{
	//std::cout<<"PreAlphaSynapse()"<<std::endl;
}

PreAlphaSynapse::PreAlphaSynapse(
		const unsigned long id,
		const number& location,
		const number& onset)

:IPreSynapse(id, location),
 m_onset(onset)
{
	//std::cout<<"PreAlphaSynapse()"<<std::endl;
}

PreAlphaSynapse::~PreAlphaSynapse() {
	//std::cout<<"~PreAlphaSynapse()"<<std::endl;
}

void PreAlphaSynapse::update(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy update atm
}


bool PreAlphaSynapse::active(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy
	return (t >= m_onset);
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
