/*
 * PreAlphaSynapse.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#include "AlphaPreSynapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

AlphaPreSynapse::AlphaPreSynapse(
		const number& location,
		const number& onset)

:IPreSynapse(0, 0, location),
 m_onset(onset)
{
}

AlphaPreSynapse::AlphaPreSynapse(
		const unsigned long long id,
		const unsigned long long postsynapse_id,
		const number& location,
		const number& onset)

:IPreSynapse(id, postsynapse_id, location),
 m_onset(onset)
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
	strs << id() << " ";
	strs << postsynapse_id() << " ";
	strs << location() << " ";
	strs << m_onset << " ";
	os << strs.str();
}

void AlphaPreSynapse::get_from(std::istream& is)
{
	unsigned long long id; is >> id; set_id(id);
	unsigned long long postsyn_id; is >> postsyn_id; set_postsynapse_id(postsyn_id);
	number loc; is >> loc; set_location(loc);
	is >> m_onset;
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
