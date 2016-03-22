/*
 * IPostSynapse.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#include "../split_synapse_handler/IPostSynapse.h"

#include <common/types.h> 											//number
#include "../synapse_handler/function/types.h" 						//SynapseType
#include <string>													//std::string
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 	//VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

IPostSynapse::IPostSynapse(
		const unsigned long long id,
		const unsigned long long presynapse_id,
		const number& location)

:m_id(id),
 m_presynapse_id(presynapse_id),
 m_location(location)
{
}

IPostSynapse::~IPostSynapse()
{
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
