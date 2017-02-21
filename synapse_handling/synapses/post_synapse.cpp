/*
 * IPostSynapse.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#include "post_synapse.h"


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

IPostSynapse::IPostSynapse()
:IBaseSynapse(0, 0)
{
}

IPostSynapse::IPostSynapse(
		const synapse_id id,
		const number location
)
:IBaseSynapse(id, location)
{
}

IPostSynapse::~IPostSynapse()
{
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
