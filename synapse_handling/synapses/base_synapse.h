/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Lukas Reinhardt
 * Creation date: 2016-03-22
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

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__BASE_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__BASE_SYNAPSE_H

#include "../types.h"
#include "common/types.h"
#include <iostream>
#include <string>

typedef unsigned long long synapse_id;

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IBaseSynapse {
private:
	synapse_id m_id;							//own presynapse id / alternative: shared id if there is a 1:1 relation?
	number m_location;							//location on edge

public:

	IBaseSynapse()
	:m_id(0),
	 m_location(0)
	{}

	IBaseSynapse(synapse_id id, number location)
	:m_id(id),
	 m_location(location)
	{}

	virtual ~IBaseSynapse(){}

	//todo: maybe protected for polymorphic ops and private for global friend ops?
	friend std::ostream& operator<<(std::ostream& os, const IBaseSynapse* s);
	friend std::istream& operator>>(std::istream& is, IBaseSynapse* s);
	virtual void put_to(std::ostream& os) const = 0;							//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;								//'get_from' == operator>>

	virtual const std::string& name() const = 0;
	virtual SynapseType type() const = 0;
	virtual bool is_presynapse() const = 0;
	virtual bool is_postsynapse() const = 0;

	void set_location(const number& loc) {m_location = loc;}
	number location() const {return m_location;}

	void set_id(const synapse_id id) {m_id = id;}
	synapse_id id() const {return m_id;}
};


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__BASE_SYNAPSE_H
