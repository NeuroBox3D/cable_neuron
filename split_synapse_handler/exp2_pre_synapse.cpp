/*
 * PreExp2Synapse.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#include "exp2_pre_synapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

/**
 * Synapse with id=0
 */
Exp2PreSynapse::Exp2PreSynapse(
		const number& location,
		const number& onset)

:IPreSynapse(0, 0, location),
 m_onset(onset)
{
}

Exp2PreSynapse::Exp2PreSynapse(
		const unsigned long long id,
		const unsigned long long postsynapse_id,
		const number& location,
		const number& onset)

:IPreSynapse(id, postsynapse_id, location),
 m_onset(onset)
{
}

Exp2PreSynapse::~Exp2PreSynapse()
{
}

void Exp2PreSynapse::update(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy update atm
}

bool Exp2PreSynapse::is_active(const number& t, VectorProxyBase* up)
{
	//todo:
	//dummy
	return (t >= m_onset);
}

void Exp2PreSynapse::put_to(std::ostream& os) const
{
	using std::ostringstream;
	ostringstream strs;
	strs << name() << " ";
	strs << id() << " ";
	strs << postsynapse_id() << " ";
	strs << location() << " ";
	strs << m_onset;
	os << strs.str();
}


void Exp2PreSynapse::get_from(std::istream& is)
{
	std::string t; is >> t;
	unsigned long long id; is >> id; set_id(id);
	unsigned long long postsyn_id; is >> postsyn_id; set_postsynapse_id(postsyn_id);
	number loc; is >> loc; set_location(loc);
	is >> m_onset;
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
