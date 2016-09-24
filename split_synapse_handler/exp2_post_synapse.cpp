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
:IPostSynapse(0, 0),
 m_onset(nan("")),
 m_gMax(0),
 m_tau1(0),
 m_tau2(0),
 m_rev(0)
 {}

Exp2PostSynapse::Exp2PostSynapse(
		const number location,
		const number onset,
		const number gMax,
		const number tau1,
		const number tau2,
		const number rev
		)
:IPostSynapse(0, location),
 m_onset(nan("")),
 m_gMax(gMax),
 m_tau1(tau1),
 m_tau2(tau2),
 m_rev(rev)
{
}

Exp2PostSynapse::Exp2PostSynapse(
		const SYNAPSE_ID id,
		const number location,
		const number onset,
		const number gMax,
		const number tau1,
		const number tau2,
		const number rev
		)
:IPostSynapse(id, location),
 m_onset(nan("")),
 m_gMax(gMax),
 m_tau1(tau1),
 m_tau2(tau2),
 m_rev(rev)
{
}

Exp2PostSynapse::~Exp2PostSynapse()
{
}

number Exp2PostSynapse::current(const number& t, const number &vm)
{
	if (t >= m_onset)	// this excludes onset == NaN
	{
		number tp = (m_tau1*m_tau2)/(m_tau2 - m_tau1) * std::log(m_tau2/m_tau1);	// time of maximal current
		number factor = 1.0 / (std::exp(-tp/m_tau2) - std::exp(-tp/m_tau1));		// normalization factor
		number i = m_gMax * factor * (vm - m_rev) * (std::exp(-(t-m_onset)/m_tau2) - std::exp(-(t-m_onset)/m_tau1));
/*
		std::cout << "Exp2PostSynapse" << id() << ":" << std::endl
				  << "location: " << location() << std::endl
				  << "onset: " << onset() << std::endl
				  << "gMax: " << gMax() << std::endl
				  << "tau1: " << tau1() << std::endl
				  << "tau2: " << tau2() << std::endl
				  << "rev:" << rev() << std::endl
				  << "is_active(" << t << "): " << is_active(t) << std::endl
				  << "current: " << i << std::endl << std::endl;
*/
		return i; //!< i: current (in units of A)
	}

	return 0.0;
}

void Exp2PostSynapse::put_to(std::ostream& os) const
{
	using std::ostringstream;
	ostringstream strs;
	strs << name() << " ";					//identifier for reconstruction
	strs << id() << " ";
	strs << location() << " ";
	strs << m_onset << " ";
	strs << m_gMax << " ";
	strs << m_tau1 << " ";
	strs << m_tau2 << " ";
	strs << m_rev;
	os << strs.str();
}

void Exp2PostSynapse::get_from(std::istream& is)
{
	using boost::lexical_cast;
	std::string tmp;

	SYNAPSE_ID id; is >> id; set_id(id);

	number loc;
	is >> tmp;
	loc = lexical_cast<number>(tmp);
	set_location(loc);
	tmp.clear();

	number onset;
	is >> tmp;
	onset = lexical_cast<number>(tmp);
	set_onset(nan(""));
	tmp.clear();

	number gMax;
	is >> tmp;
	gMax = lexical_cast<number>(tmp);
	set_gMax(gMax);
	tmp.clear();

	number tau1;
	is >> tmp;
	tau1 = lexical_cast<number>(tmp);
	set_tau1(tau1);
	tmp.clear();

	number tau2;
	is >> tmp;
	tau2 = lexical_cast<number>(tmp);
	set_tau2(tau2);
	tmp.clear();

	number rev;
	is >> tmp;
	rev = lexical_cast<number>(tmp);
	set_rev(rev);
	tmp.clear();

}

bool Exp2PostSynapse::is_active(number time)
{
	return m_onset == m_onset; //if false: m_onset is nan -> inactive
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
