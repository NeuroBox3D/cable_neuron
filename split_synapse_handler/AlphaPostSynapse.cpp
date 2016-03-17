/*
 * PostAlphaSynapse.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#include "AlphaPostSynapse.h"

#include <cmath>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

AlphaPostSynapse::AlphaPostSynapse(
		const number& location,
		const number& gMax,
		const number& onset,
		const number& tau,
		const number& vm,
		const number& e)

:IPostSynapse(0, 0, location),
 m_gMax(gMax),
 m_onset(onset),
 m_tau(tau),
 m_vm(vm),
 m_e(e)
{
}

AlphaPostSynapse::AlphaPostSynapse(
		const unsigned long id,
		const unsigned long presynapse_id,
		const number& location,
		const number& gMax,
		const number& onset,
		const number& tau,
		const number& vm,
		const number& e)

:IPostSynapse(id, presynapse_id, location),
 m_gMax(gMax),
 m_onset(onset),
 m_tau(tau),
 m_vm(vm),
 m_e(e)
{
}

AlphaPostSynapse::~AlphaPostSynapse()
{
}

number AlphaPostSynapse::current(const number& t)
{
	return m_gMax * (t - m_onset)/m_tau * std::exp(-(t - m_onset - m_tau)/m_tau) * (m_vm - m_e);	// current (in units of A)
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
