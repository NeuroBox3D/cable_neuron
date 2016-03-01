/*
 * IPreSynapse.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: lreinhardt
 */

#include "IPreSynapse.h"
#include <iostream>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {


IPreSynapse::
IPreSynapse(number location)
:m_location(location)
{
	std::cout<<"IPreSynapse()"<<std::endl;
}


IPreSynapse::~IPreSynapse()
{
	std::cout<<"~IPreSynapse()"<<std::endl;
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
