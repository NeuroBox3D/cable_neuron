/*
 * PostExp2Synapse.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#include "exp2_post_synapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

Exp2PostSynapse::Exp2PostSynapse()
:IPostSynapse(0, 0, 0),
 m_tau1(0),
 m_tau2(0),
 m_e(0),
 m_w(0),
// m_vm(0)
 m_onset(0)
{}

Exp2PostSynapse::Exp2PostSynapse(
		const number& location,
		const number& onset,
		const number& tau1,
		const number& tau2,
		const number& e,
		const number& w
//		const number& vm
		)

:IPostSynapse(0, 0, location),
 m_tau1(tau1),
 m_tau2(tau2),
 m_e(e),
 m_w(w),
// m_vm(vm)
 m_onset(onset)
{
}

Exp2PostSynapse::Exp2PostSynapse(
		const unsigned long long id,
		const unsigned long long presynapse_id,
		const number& onset,
		const number& location,
		const number& tau1,
		const number& tau2,
		const number& e,
		const number& w
//		const number& vm
		)

:IPostSynapse(id, presynapse_id, location),
 m_tau1(tau1),
 m_tau2(tau2),
 m_e(e),
 m_w(w),
// m_vm(vm)
 m_onset(onset)
{
}

Exp2PostSynapse::~Exp2PostSynapse()
{
}

number Exp2PostSynapse::current(const number& t, const number &vm)
{
	number tp = (m_tau1*m_tau2)/(m_tau2 - m_tau1) * std::log(m_tau2/m_tau1);	// time of maximal current
	number factor = 1.0 / (std::exp(-tp/m_tau2) - std::exp(-tp/m_tau1));		// normalization factor
	number i = m_w * factor * (vm - m_e) * (std::exp(-t/m_tau2) - std::exp(-t/m_tau1));
	return i; //!< i: current (in units of A)
}

void Exp2PostSynapse::put_to(std::ostream& os) const
{
	using std::ostringstream;
	ostringstream strs;
	strs << name() << " ";					//identifier for reconstruction
	strs << id() << " ";
	strs << presynapse_id() << " ";
	strs << location() << " ";
	strs << m_onset << " ";
	strs << m_tau1 << " ";
	strs << m_tau2 << " ";
	strs <<  m_e << " ";
	strs << m_w << " ";
//	strs << m_vm;
	os << strs.str();
}

void Exp2PostSynapse::get_from(std::istream& is)
{
	//std::string t; is >> t;
	unsigned long long id; is >> id; set_id(id);
	unsigned long long presyn_id; is >> presyn_id; set_presynapse_id(presyn_id);
	number loc; is >> loc; set_location(loc);
	is >> m_onset;
	is >> m_tau1;
	is >> m_tau2;
	is >> m_e;
	is >> m_w;
//	is >> m_vm;
}

bool Exp2PostSynapse::is_active(number time)
{
	return m_onset == m_onset; //if false: m_onset is nan -> inactive
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
