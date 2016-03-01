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

PreAlphaSynapse::PreAlphaSynapse(number location)
:IPreSynapse(location)
{
	std::cout<<"PreAlphaSynapse()"<<std::endl;

}

PreAlphaSynapse::~PreAlphaSynapse() {
	std::cout<<"~PreAlphaSynapse()"<<std::endl;
}

void
PreAlphaSynapse::
update(number t, VectorProxyBase* up)
{
	//dummy update atm
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
