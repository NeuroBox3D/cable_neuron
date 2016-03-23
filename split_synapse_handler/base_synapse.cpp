/*
 * IBaseSynapse.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#include "base_synapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

std::ostream& operator<<(std::ostream& os, const IBaseSynapse* s)
{
	s->put_to(os);
	return os;
}

std::istream& operator>>(std::istream& is, IBaseSynapse* s)
{
	s->get_from(is);
	return is;
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
