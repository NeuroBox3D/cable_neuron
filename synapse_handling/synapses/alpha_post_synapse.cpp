/*
 * AlphaPostSynapse.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#include "alpha_post_synapse.h"

#include <cmath>
#include <limits> // numeric_limits
#include <boost/lexical_cast.hpp>


namespace ug {
namespace cable_neuron {
namespace synapse_handler {


const std::string AlphaPostSynapse::name_string = "AlphaPostSynapse";


AlphaPostSynapse::AlphaPostSynapse()
:IPostSynapse(0, 0),
 m_onset(std::numeric_limits<number>::quiet_NaN()),
 m_gMax(0),
 m_tau(0),
 m_rev(0)
{}

AlphaPostSynapse::AlphaPostSynapse(
		const number location,
		const number gMax,
		const number tau,
		const number rev)

:IPostSynapse(0, location),
 m_onset(std::numeric_limits<number>::quiet_NaN()),
 m_gMax(gMax),
 m_tau(tau),
 m_rev(rev)
{}

AlphaPostSynapse::AlphaPostSynapse(
		const synapse_id id,
		const number location,
		const number gMax,
		const number tau,
		const number rev)

:IPostSynapse(id, location),
 m_onset(std::numeric_limits<number>::quiet_NaN()),
 m_gMax(gMax),
 m_tau(tau),
 m_rev(rev)
{}

AlphaPostSynapse::~AlphaPostSynapse()
{}


number AlphaPostSynapse::current(const number& t, const number& vm)
{

	if (t >= m_onset) {	// this excludes onset == NaN
		double curr = m_gMax * (t - m_onset) / m_tau
		              * std::exp(-(t - m_onset - m_tau)/m_tau) * (vm - m_rev);	// current (in units of A);
//		std::cout << "AlphaPostSynapse" << id() << ":" << std::endl
//				  << "location: " << location() << std::endl
//				  << "onset: " << onset() << std::endl
//				  << "gMax: " << gMax() << std::endl
//				  << "tau: " << tau() << std::endl
//				  << "rev:" << rev() << std::endl
//				  << "is_active(" << t << "): " << is_active(t) << std::endl
//				  << "current: " << curr << std::endl << std::endl;
		return curr;
	}

	return 0.0;
}


void AlphaPostSynapse::put_to(std::ostream& os) const
{
	using std::ostringstream;
	ostringstream strs;
	strs << name() << " ";					//identifier for reconstruction
	strs << id() << " ";
	strs << location() << " ";
	strs << m_gMax << " ";
	strs << m_tau << " ";
	strs << m_rev;
	os << strs.str();
}


void AlphaPostSynapse::get_from(std::istream& is)
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

	number gMax;
	is >> tmp;
	gMax = lexical_cast<number>(tmp);
	set_gMax(gMax);
	tmp.clear();

	number tau;
	is >> tmp;
	tau = lexical_cast<number>(tmp);
	set_tau(tau);
	tmp.clear();

	number rev;
	is >> tmp;
	rev = lexical_cast<number>(tmp);
	set_rev(rev);
	tmp.clear();

}

void AlphaPostSynapse::activate(number time)
{
    m_onset = time;
}

void AlphaPostSynapse::deactivate()
{
    m_onset = std::numeric_limits<number>::quiet_NaN();
}


bool AlphaPostSynapse::is_active(number time)
{
	return m_onset == m_onset;
}

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
