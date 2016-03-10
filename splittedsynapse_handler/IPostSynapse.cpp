/*
 * IPostSynapse.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#include "IPostSynapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

IPostSynapse::
IPostSynapse(const number& gMax, const number& tau, const number& vm, const number& e)
:m_gMax(gMax),m_tau(tau),m_vm(vm),m_e(e),m_pre_synapse_id(0)
{
}

IPostSynapse::~IPostSynapse() {
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
