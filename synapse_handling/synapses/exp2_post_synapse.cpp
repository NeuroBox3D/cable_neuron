/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Lukas Reinhardt
 * Creation date: 2016-03-16
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "exp2_post_synapse.h"

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <limits> // numeric_limits


namespace ug {
namespace cable_neuron {
namespace synapse_handler {


const std::string Exp2PostSynapse::name_string = "Exp2PostSynapse";


Exp2PostSynapse::Exp2PostSynapse()
: IPostSynapse(0, 0),
  m_onset(std::numeric_limits<number>::quiet_NaN()),
  m_gMax(0),
  m_tau1(0),
  m_tau2(0),
  m_rev(0)
{}

Exp2PostSynapse::Exp2PostSynapse
(
    const number location,
    const number gMax,
    const number tau1,
    const number tau2,
	const number rev
)
: IPostSynapse(0, location),
  m_onset(std::numeric_limits<number>::quiet_NaN()),
  m_gMax(gMax),
  m_tau1(tau1),
  m_tau2(tau2),
  m_rev(rev)
{}

Exp2PostSynapse::Exp2PostSynapse
(
    const synapse_id id,
    const number location,
    const number gMax,
    const number tau1,
    const number tau2,
    const number rev
)
: IPostSynapse(id, location),
  m_onset(std::numeric_limits<number>::quiet_NaN()),
  m_gMax(gMax),
  m_tau1(tau1),
  m_tau2(tau2),
  m_rev(rev)
{}

Exp2PostSynapse::~Exp2PostSynapse()
{}

number Exp2PostSynapse::current(const number& t, const number &vm) const
{
	if (t >= m_onset)	// this excludes onset == NaN
	{
		// in case tau1 = tau2, the current degenerates to alpha current
		if (fabs(1.0 - m_tau1/m_tau2) < 1e-8)
			return m_gMax * (vm - m_rev) * (t-m_onset)/m_tau2 * std::exp(-(t-m_onset-m_tau2)/m_tau2);

		number tp = (m_tau1*m_tau2)/(m_tau2 - m_tau1) * std::log(m_tau2/m_tau1);	// time of maximal current
		number factor = 1.0 / (std::exp(-tp/m_tau2) - std::exp(-tp/m_tau1));		// normalization factor
		return m_gMax * factor * (vm - m_rev) * (std::exp(-(t-m_onset)/m_tau2) - std::exp(-(t-m_onset)/m_tau1));

#if 0
		std::cout << "Exp2PostSynapse " << id() << ":" << std::endl
				  << "location: " << location() << std::endl
				  << "onset: " << onset() << std::endl
				  << "gMax: " << gMax() << std::endl
				  << "tau1: " << tau1() << std::endl
				  << "tau2: " << tau2() << std::endl
				  << "rev:" << rev() << std::endl
				  << "is_active(" << t << "): " << is_active(t) << std::endl
				  << "current: " << i << std::endl << std::endl;
#endif
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


void Exp2PostSynapse::activate(number time) {m_onset = time;}
void Exp2PostSynapse::deactivate() {m_onset = std::numeric_limits<number>::quiet_NaN();}


bool Exp2PostSynapse::is_active(number time) const
{
	return m_onset == m_onset; //if false: m_onset is nan -> inactive
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
