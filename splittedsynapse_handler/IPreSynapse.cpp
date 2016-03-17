/*
 * IPreSynapse.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: lreinhardt
 */

#include "IPreSynapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

/**
 * Synapse with id=0
 */
IPreSynapse::IPreSynapse(const number& location)

:m_id(0),
 m_location(location)
{
	//std::cout<<"IPreSynapse()"<<std::endl;
}

IPreSynapse::IPreSynapse(
		const unsigned long id,
		const number& location)

:m_id(id),
 m_location(location)
{
	//std::cout<<"IPreSynapse()"<<std::endl;
}


IPreSynapse::~IPreSynapse()
{
	//std::cout<<"~IPreSynapse()"<<std::endl;
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
