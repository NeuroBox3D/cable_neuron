/*
 * IPreSynapse.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: lreinhardt
 */

#include "../split_synapse_handler/IPreSynapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

IPreSynapse::IPreSynapse()
:m_id(0),
 m_postsynapse_id(0),
 m_location(0.0)
{

}


/**
 * Synapse with id's=0
 */
IPreSynapse::IPreSynapse(const number& location)

:m_id(0),
 m_postsynapse_id(0),
 m_location(location)
{
}

IPreSynapse::IPreSynapse(
		const unsigned long id,
		const unsigned long postsynapse_id,
		const number& location)

:m_id(id),
 m_postsynapse_id(postsynapse_id),
 m_location(location)
{
}


IPreSynapse::~IPreSynapse()
{
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
