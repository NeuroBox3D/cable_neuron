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

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__EXP2_POST_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__EXP2_POST_SYNAPSE_H

#include "post_synapse.h"   // IPostSynapse


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class Exp2PostSynapse : public IPostSynapse
{
    public:
        static const std::string name_string;

private:
	number m_onset; ///<onset (s)
	number m_gMax;  ///< conductance (S)
	number m_tau1;  ///< time constant (s)
	number m_tau2;  ///< time constant (s)
	number m_rev;   ///< reversal potential (V)


public:
	//ctor & dtor
	Exp2PostSynapse();  // needed for template generator
	Exp2PostSynapse(
			const number location,
			const number gMax,
			const number tau1,
			const number tau2,
			const number rev
			);

	Exp2PostSynapse(
			const synapse_id id,
			const number location,
			const number gMax,
			const number tau1,
			const number tau2,
			const number rev
			);

	virtual ~Exp2PostSynapse();

	//setter & getter
	void set_onset(number onset) {m_onset=onset;}
	void set_gMax(number gmax) {m_gMax=gmax;}
	void set_tau1(number tau1) {m_tau1=tau1;}
	void set_tau2(number tau2) {m_tau2=tau2;}
	void set_rev(number rev) {m_rev=rev;}


	number onset() const {return m_onset;}
	number gMax() const {return m_gMax;}
	number tau1() const {return m_tau1;}
	number tau2() const {return m_tau2;}
	number rev() const {return m_rev;}


	SynapseType type() const {return EXP2_POST_SYNAPSE;}
	const std::string& name() const {return name_string;}

	//functionality
	number current(const number& t, const number& vm) const;

	//serialization interface methods
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>

	void activate(number time);
	void deactivate();

	bool is_active(number time) const;
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSES__EXP2_POST_SYNAPSE_H
