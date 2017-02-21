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
:IBaseSynapse(0,0)
{
}

IPreSynapse::IPreSynapse(
		const synapse_id id,
		const number location
)
:IBaseSynapse(id, location)
{
}


IPreSynapse::~IPreSynapse()
{
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
