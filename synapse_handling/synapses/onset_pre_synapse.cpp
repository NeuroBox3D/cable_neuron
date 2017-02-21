/*
 * AlphaPreSynapse.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#include "onset_pre_synapse.h"

#include <boost/lexical_cast.hpp>


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

const std::string OnsetPreSynapse::name_string = "OnsetPreSynapse";


OnsetPreSynapse::OnsetPreSynapse()
:IPreSynapse(0, 0),
 m_onset(0),m_duration(0)
{}

OnsetPreSynapse::OnsetPreSynapse(
		const number location,
		const number onset,
		const number duration)

:IPreSynapse(0, location),
 m_onset(onset),m_duration(duration)
{
}

OnsetPreSynapse::OnsetPreSynapse(
		const synapse_id id,
		const number location,
		const number onset,
		const number duration)

:IPreSynapse(id, location),
 m_onset(onset),m_duration(duration)
{
}

OnsetPreSynapse::~OnsetPreSynapse()
{
}

void OnsetPreSynapse::update(const number& t, const std::vector<number>& u)
{
	//interface dummy
}


bool OnsetPreSynapse::is_active(const number& t)
{
	//t in [onset, onset+duration] ?
	return (t >= m_onset) && (t <= (m_onset + m_duration) ) ;
}

void OnsetPreSynapse::put_to(std::ostream& os) const
{
	using std::ostringstream;
	ostringstream strs;
	strs << name() << " ";
	strs << id() << " ";
	strs << location() << " ";
	strs << m_onset << " ";
	strs << m_duration;
	os << strs.str();
}

void OnsetPreSynapse::get_from(std::istream& is)
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


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
