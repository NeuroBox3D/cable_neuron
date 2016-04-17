/*
 * PostAlphaSynapse.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#include "alpha_post_synapse.h"

#include <cmath>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

AlphaPostSynapse::AlphaPostSynapse()
:IPostSynapse(0, 0, 0),
 m_gMax(1.2e-3),
 m_onset(0),
 m_tau(0),
// m_vm(0),
 m_e(0)
{}

AlphaPostSynapse::AlphaPostSynapse(
		const number& location,
		const number& gMax,
		const number& onset,
		const number& tau,
//		const number& vm,
		const number& e)

:IPostSynapse(0, 0, location),
 m_gMax(gMax),
 m_onset(onset),
 m_tau(tau),
// m_vm(vm),
 m_e(e)
{
}

AlphaPostSynapse::AlphaPostSynapse(
		const unsigned long long id,
		const unsigned long long presynapse_id,
		const number& location,
		const number& gMax,
		const number& onset,
		const number& tau,
//		const number& vm,
		const number& e)

:IPostSynapse(id, presynapse_id, location),
 m_gMax(gMax),
 m_onset(onset),
 m_tau(tau),
// m_vm(vm),
 m_e(e)
{
}

AlphaPostSynapse::~AlphaPostSynapse()
{
}

number AlphaPostSynapse::current(const number& t, const number& vm)
{
	return m_gMax * (t - m_onset)/m_tau * std::exp(-(t - m_onset - m_tau)/m_tau) * (vm - m_e);	// current (in units of A)
}

/**
 * ploymorphic part of << serialization
 */
void AlphaPostSynapse::put_to(std::ostream& os) const
{
	using std::ostringstream;
	ostringstream strs;
	strs << name() << " ";					//identifier for reconstruction
	strs << id() << " ";
	strs << presynapse_id() << " ";
	strs << location() << " ";
	strs << m_gMax << " ";
	strs << m_onset << " ";
	strs <<  m_tau << " ";
//	strs << m_vm << " ";
	strs << m_e;
	os << strs.str();
}

/**
 * ploymorphic part of >> deserialization
 */
void AlphaPostSynapse::get_from(std::istream& is)
{
	//std::string t; is >> t;
	unsigned long long id; is >> id; set_id(id);
	unsigned long long presyn_id; is >> presyn_id; set_presynapse_id(presyn_id);
	number loc; is >> loc; set_location(loc);
	is >> m_gMax;
	is >> m_onset;
	is >> m_tau;
//	is >> m_vm;
	is >> m_e;
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
