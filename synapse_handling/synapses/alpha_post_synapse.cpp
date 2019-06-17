/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Lukas Reinhardt
 * Creation date: 2016-03-10
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


number AlphaPostSynapse::current(const number& t, const number& vm) const
{

	if (t >= m_onset) {	// this excludes onset == NaN
		double curr = m_gMax * (t - m_onset) / m_tau
		              * std::exp(-(t - m_onset - m_tau)/m_tau) * (vm - m_rev);	// current (in units of A);
#if 0
		// debug
		std::cout << "AlphaPostSynapse" << id() << ":" << std::endl
				  << "location: " << location() << std::endl
				  << "onset: " << onset() << std::endl
				  << "gMax: " << gMax() << std::endl
				  << "tau: " << tau() << std::endl
				  << "rev:" << rev() << std::endl
				  << "is_active(" << t << "): " << is_active(t) << std::endl
				  << "current: " << curr << std::endl << std::endl;
#endif
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


bool AlphaPostSynapse::is_active(number time) const
{
	return m_onset == m_onset;
}

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
