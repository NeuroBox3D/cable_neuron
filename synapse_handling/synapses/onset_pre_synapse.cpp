/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Lukas Reinhardt
 * Creation date: 2016-03-01
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

#include "onset_pre_synapse.h"

#include <boost/lexical_cast.hpp>
#include <limits>


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

const std::string OnsetPreSynapse::name_string = "OnsetPreSynapse";


OnsetPreSynapse::OnsetPreSynapse()
:IPreSynapse(0, 0),
 m_onset(std::numeric_limits<number>::max()), m_duration(0)
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
