/*
 * PostExp2Synapse.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#include "Exp2PostSynapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

PostExp2Synapse::PostExp2Synapse(
		const number& location,
		const number& tau1,
		const number& tau2,
		const number& e,
		const number& w,
		const number& vm)

:IPostSynapse(0, 0, location),
 m_tau1(tau1),
 m_tau2(tau2),
 m_e(e),
 m_w(w),
 m_vm(vm)
{
}

PostExp2Synapse::PostExp2Synapse(
		const unsigned long id,
		const unsigned long presynapse_id,
		const number& location,
		const number& tau1,
		const number& tau2,
		const number& e,
		const number& w,
		const number& vm)

:IPostSynapse(id, presynapse_id, location),
 m_tau1(tau1),
 m_tau2(tau2),
 m_e(e),
 m_w(w),
 m_vm(vm)
{
}

PostExp2Synapse::~PostExp2Synapse()
{
}

number PostExp2Synapse::current(const number& t)
{
	number tp = (m_tau1*m_tau2)/(m_tau2 - m_tau1) * std::log(m_tau2/m_tau1);	// time of maximal current
	number factor = 1.0 / (std::exp(-tp/m_tau2) - std::exp(-tp/m_tau1));		// normalization factor
	number i = m_w * factor * (m_vm - m_e) * (std::exp(-t/m_tau2) - std::exp(-t/m_tau1));
	return i; //!< i: current (in units of A)
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
