/*
 * PostAlphaSynapse.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#include "PostAlphaSynapse.h"
#include <cmath>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

PostAlphaSynapse::
PostAlphaSynapse(const number& gMax, const number& onset, const number& tau, const number& tau, const number& e)
:m_gMax(gMax), m_onset(onset), m_tau(tau), m_vm(vm), m_e(e)
{
	// TODO Auto-generated constructor stub

}

PostAlphaSynapse::
~PostAlphaSynapse() {
	// TODO Auto-generated destructor stub
}

number
PostAlphaSynapse::
current(const number& t)
{
	return m_gMax * (t - m_onset)/m_tau * std::exp(-(t - m_onset - m_tau)/m_tau) * (m_vm - m_e);	// current (in units of A)
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
