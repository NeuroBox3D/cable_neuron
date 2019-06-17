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

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__THRESHOLD_PRE_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__THRESHOLD_PRE_SYNAPSE_H

#include <string>           // std::string
#include "pre_synapse.h"    // IPreSynapse

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class ThresholdPreSynapse : public IPreSynapse
{
public:
    static const std::string name_string;

private:
	number m_onset;
	number m_duration;
	number m_threshold;

public:
	//ctor & dtor
	ThresholdPreSynapse(); // needed for template generator
	ThresholdPreSynapse(
			const number location,
			const number duration,
			const number threshold
			);

	ThresholdPreSynapse(
			const synapse_id id,
			const number location,
			const number duration,
			const number threshold
			);

	virtual ~ThresholdPreSynapse();

	// setter & getter
	void set_onset(number onset) {m_onset = onset;}
	void set_duration(number duration) {m_duration = duration;}
	void set_threshold(number val) {m_threshold = val;}

	number onset() const {return m_onset;}
	number duration() const {return m_duration;}
	number threshold() const {return m_threshold;}

	SynapseType type() const {return THRESHOLD_PRE_SYNAPSE;}
	const std::string& name() const {return name_string;}

	//interface methods
	void update(const number& t, const std::vector<number>& u);
	bool is_active(const number& t);

	//serialization interface methods
	void put_to(std::ostream& os) const;        //'put_to' == operator<<
	void get_from(std::istream& is);            //'get_from' == operator>>
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__THRESHOLD_PRE_SYNAPSE_H
