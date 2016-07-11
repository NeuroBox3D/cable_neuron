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
		const SYNAPSE_ID id,
		const SYNAPSE_ID presynapse_id,
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
	if (t >= m_onset)	// this excludes onset == NaN
		return m_gMax * (t - m_onset)/m_tau * std::exp(-(t - m_onset - m_tau)/m_tau) * (vm - m_e);	// current (in units of A)

	return 0.0;
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
	using boost::lexical_cast;
	std::string tmp;
	//std::string t; is >> t;
	SYNAPSE_ID id; is >> id; set_id(id);
	SYNAPSE_ID presyn_id; is >> presyn_id; set_presynapse_id(presyn_id);

	number loc;
	is >> tmp;
	loc = lexical_cast<number>(tmp);
	set_location(loc);
	tmp.clear();

	number gMax;
	is >> tmp;
	gMax = lexical_cast<number>(tmp);
	set_gMax(gMax);
	tmp.clear();

	number onset;
	is >> tmp;
	onset = lexical_cast<number>(tmp);
	set_onset(onset);
	tmp.clear();

	number tau;
	is >> tmp;
	tau = lexical_cast<number>(tmp);
	set_tau(tau);
	tmp.clear();

	number e;
	is >> tmp;
	e = lexical_cast<number>(tmp);
	set_e(e);
	tmp.clear();

//	is >> m_vm;
}

bool AlphaPostSynapse::is_active(number time)
{
	return m_onset == m_onset;
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
