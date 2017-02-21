/*
 * IBaseSynapse.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#include "base_synapse.h"
#include "common/error.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

std::ostream& operator<<(std::ostream& os, const IBaseSynapse* s)
{
	try	{s->put_to(os);}
	UG_CATCH_THROW("Could not write synapse to stream.");

	return os;
}

std::istream& operator>>(std::istream& is, IBaseSynapse* s)
{
	try {s->get_from(is);}
	UG_CATCH_THROW("Could not read synapse from stream.");
	return is;
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
