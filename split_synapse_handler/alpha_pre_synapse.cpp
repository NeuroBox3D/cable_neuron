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

void AlphaPreSynapse::update(const number& t, const std::vector<number>& u)
{
	//todo:
	//dummy update atm
}


bool AlphaPreSynapse::is_active(const number& t)
{
	//todo:
	//t in [onset, onset+duration] ?
	return (t >= m_onset) && (t <= (m_onset + m_duration) ) ;
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
	using boost::lexical_cast;
	std::string tmp;
	//std::string t; is >> t;
	SYNAPSE_ID id; is >> id; set_id(id);
	SYNAPSE_ID postsyn_id; is >> postsyn_id; set_postsynapse_id(postsyn_id);

	number loc;
	is >> tmp;
	loc = lexical_cast<number>(tmp);
	set_location(loc);
	tmp.clear();

	number onset;
	is >> tmp;
	onset = lexical_cast<number>(tmp);
	set_onset(onset);
	tmp.clear();

	number duration;
	is >> tmp;
	duration = lexical_cast<number>(tmp);
	set_duration(duration);
	tmp.clear();
}

//bool AlphaPreSynapse::fire(number time, SYNAPSE_ID& postsyn_id)
//{
//	if(time >= m_onset) {
//		postsyn_id = postsynapse_id();
//		return true;
//	}
//	return false;
//}
//
//bool AlphaPreSynapse::cooldown(number time, SYNAPSE_ID& postsyn_id)
//{
//	number dt = time - m_onset;
//	if( (time >= m_onset) && (dt >= m_duration)) {
//		postsyn_id = postsynapse_id();
//		return true;
//	}
//	return false;
//}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
