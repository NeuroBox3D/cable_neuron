/*
 * IPreSynapse.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: lreinhardt
 */

#include "pre_synapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

IPreSynapse::IPreSynapse()
:IBaseSynapse(0),
 m_id(0),
 m_postsynapse_id(0)
{
}

IPreSynapse::IPreSynapse(
		const unsigned long long id,
		const unsigned long long postsynapse_id,
		const number& location)

:IBaseSynapse(location),
 m_id(id),
 m_postsynapse_id(postsynapse_id)
{
}


IPreSynapse::~IPreSynapse()
{
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
