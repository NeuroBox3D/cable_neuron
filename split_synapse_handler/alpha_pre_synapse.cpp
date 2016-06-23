/*
 * PreAlphaSynapse.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#include "alpha_pre_synapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

AlphaPreSynapse::AlphaPreSynapse()
:IPreSynapse(0, 0, 0),
 m_onset(0),m_duration(0)
{}

AlphaPreSynapse::AlphaPreSynapse(
		const number& location,
		const number& onset,
		const number& duration)

:IPreSynapse(0, 0, location),
 m_onset(onset),m_duration(duration)
{
}

AlphaPreSynapse::AlphaPreSynapse(
		const SYNAPSE_ID id,
		const SYNAPSE_ID postsynapse_id,
		const number& location,
		const number& onset,
		const number& duration)

:IPreSynapse(id, postsynapse_id, location),
 m_onset(onset),m_duration(duration)
{
}

AlphaPreSynapse::~AlphaPreSynapse()
{
}

void AlphaPreSynapse::update(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy update atm
}


bool AlphaPreSynapse::is_active(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy
	return (t >= m_onset);
}

void AlphaPreSynapse::put_to(std::ostream& os) const
{
	using std::ostringstream;
	ostringstream strs;
	strs << name() << " ";
	strs << id() << " ";
	strs << postsynapse_id() << " ";
	strs << location() << " ";
	strs << m_onset << " ";
	strs << m_duration;
	os << strs.str();
}

void AlphaPreSynapse::get_from(std::istream& is)
{
	//std::string t; is >> t;
	SYNAPSE_ID id; is >> id; set_id(id);
	SYNAPSE_ID postsyn_id; is >> postsyn_id; set_postsynapse_id(postsyn_id);
	number loc; is >> loc; set_location(loc);
	is >> m_onset;
	is >> m_duration;
}

bool AlphaPreSynapse::fire(number time, SYNAPSE_ID& postsyn_id)
{
	if(time >= m_onset) {
		postsyn_id = postsynapse_id();
		return true;
	}
	return false;
}

bool AlphaPreSynapse::cooldown(number time, SYNAPSE_ID& postsyn_id)
{
	number dt = time - m_onset;
	if( (time >= m_onset) && (dt >= m_duration)) {
		postsyn_id = postsynapse_id();
		return true;
	}
	return false;
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
