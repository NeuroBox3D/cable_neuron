/*
 * PreExp2Synapse.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#include "threshold_pre_synapse.h"

#include <boost/lexical_cast.hpp>
#include <limits>


namespace ug {
namespace cable_neuron {
namespace synapse_handler {


const std::string ThresholdPreSynapse::name_string = "ThresholdPreSynapse";


ThresholdPreSynapse::ThresholdPreSynapse()
:IPreSynapse(0, 0),
 m_onset(std::numeric_limits<number>::quiet_NaN()), m_duration(0), m_threshold(0)
{
}

/**
 * Synapse with id=0
 */
ThresholdPreSynapse::ThresholdPreSynapse(
		const number location,
		const number duration,
		const number threshold
		)

:IPreSynapse(0, location),
 m_onset(std::numeric_limits<number>::quiet_NaN()), m_duration(duration), m_threshold(threshold)
{
}

ThresholdPreSynapse::ThresholdPreSynapse(
		const synapse_id id,
		const number location,
		const number duration,
		const number threshold
		)
:IPreSynapse(id, location),
 m_onset(std::numeric_limits<number>::quiet_NaN()), m_duration(duration), m_threshold(threshold)
{
}

ThresholdPreSynapse::~ThresholdPreSynapse()
{
}

void ThresholdPreSynapse::update(const number& t, const std::vector<number>& u)
{
	// TODO: at the moment, u contains exactly one value, which is the potential at the synapse location
	// this might be changed in the future!

	// activate if potential above threshold
	if (u[0] > m_threshold)
		m_onset = t;

	// deactivate if duration expired
	if (m_onset + m_duration < t)
		m_onset = std::numeric_limits<number>::quiet_NaN();
}

bool ThresholdPreSynapse::is_active(const number& t)
{
	return m_onset == m_onset;
}

/**
 * serialisation
 */
void ThresholdPreSynapse::put_to(std::ostream& os) const
{
	using std::ostringstream;
	ostringstream strs;
	strs << name() << " ";
	strs << id() << " ";
	strs << location() << " ";
	strs << m_duration << " ";
	strs << m_threshold;
	os << strs.str();
}

/**
 * deserialisation
 */
void ThresholdPreSynapse::get_from(std::istream& is)
{
	using boost::lexical_cast;
	std::string tmp;

	//std::string t; is >> t;
	synapse_id id; is >> id; set_id(id);
	number loc;
	is >> tmp;
	loc = lexical_cast<number>(tmp);
	set_location(loc);
	tmp.clear();

	number duration;
	is >> tmp;
	duration = lexical_cast<number>(tmp);
	set_duration(duration);
	tmp.clear();

	number threshold;
	is >> tmp;
	threshold = lexical_cast<number>(tmp);
	set_threshold(threshold);
}

//bool Exp2PreSynapse::fire(number time, SYNAPSE_ID& postsyn_id)
//{
//	if(time >= m_onset) {
//		postsyn_id = postsynapse_id();
//		return true;
//	}
//	return false;
//}
//
//bool Exp2PreSynapse::cooldown(number time, SYNAPSE_ID& postsyn_id)
//{
//	number dt = time - m_onset;
//	if( (time >= m_onset) && (dt >= m_duration)) {
//		postsyn_id = postsynapse_id();
//		return true;
//	}
//	return false;
//}

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
