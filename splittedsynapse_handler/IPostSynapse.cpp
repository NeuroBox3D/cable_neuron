/*
 * IPostSynapse.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#include "IPostSynapse.h"
#include <common/types.h> //number

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

IPostSynapse::
IPostSynapse()
:m_presynapse_id(0)
{
}

IPostSynapse::~IPostSynapse() {
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
