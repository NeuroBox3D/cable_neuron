/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Pascal Gottmann
 * Creation date: 2014-10-29
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

#ifndef UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__CABLE_MEMBRANE_TRANSPORT_INTERFACE_H
#define UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__CABLE_MEMBRANE_TRANSPORT_INTERFACE_H


#include "../cable_disc/cable_equation.h"


namespace ug {
namespace cable_neuron {


// forward declaration
template <typename TDomain>
class CableEquation;


template <typename TDomain>
class ICableMembraneTransport
{
	public:
		// TODO: We do not need functions here!
		///	constructor with comma-separated c-string
		ICableMembraneTransport(const char* functions, const char* subsets);

		/// constructor with vector of string
		ICableMembraneTransport(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		///	destructor
		virtual ~ICableMembraneTransport() {};

		/// name
		virtual std::string name() {return std::string("unknown");};

		/**
		 * @brief Initializes the defined channel type.
		 *
		 * During the initialization, the necessary attachments are attached to the vertices
		 * and their values calculated by the equilibrium state for the initial membrane potential
		 * given in the grid function passed
		 *
		 * @param time			initial time
		 * @param spGridFct		initial solution (containing membrane potential and ion concentrations)
		 */
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values) = 0;

		/// updates the gating parameters
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) = 0;

		/// provides the current densities for flowing quantities (C/(m^2*s) or mol/(m^2*s)) at a given vertex
		virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) = 0;

		/// called when approximation space is available
		void approx_space_available();

		/// called when access to the underlying CableEquation object is possible
		virtual void ce_obj_available() {};

		/// getting values of internal channel states
		virtual std::vector<number> state_values(number x, number y, number z) const
		{return std::vector<number>(0);}

		// TODO: think about generalizing this to a real Jacobian;
		// and about making implementation of this method mandatory (for time step security!)
		/**
		 * 	@brief calculate (estimate) linear dependency on potential
		 *
		 *	This method is useful for automatic time step size calculation. It is supposed to calculate
		 *	(or at least estimate) the linear dependency (first term in the taylor expansion)
		 *	of the electric current through the mechanism on the given potential. This is used in
		 *	the CableEquation to get an estimate for the CFL condition that has to be fulfilled by the
		 *	time step size.
		 *	The default implementation in this base class is to call the current() function with the given
		 *	value for the potential (if at all dependent on it) and with a value slightly bigger. Both values
		 *	are then used to estimate the dependency through a finite difference.
		 *
		 * @param vrt			vertex to compute dependency for
		 * @param vrt_values	current solution at this vertex
		 * @return				linear dependency
		 */
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);


		const std::vector<size_t>& fct_indices() const {return m_vWFctInd;}
		const std::vector<std::string>& write_subsets() {return m_vSubset;}

		/**
		 * @brief check definition on a given subset
		 *
		 * @param	si subset index to compare against
		 * @return	true iff channel is defined on subset index si
		 */
		bool is_def_on_subset(int si) const;

		void set_ce_object(CableEquation<TDomain>* pCE) {m_pCE = pCE;}

	private:
		virtual void specify_write_function_indices() {UG_THROW("specify_write_function_indices() not implemented!");}

	protected:
		void subsetCString2Vector(std::vector<std::string>& outVec, const char* cstr);
		void subsetNames2Indices(std::vector<int>& ind, const std::vector<std::string>& names);

	protected:
		/// indices in underlying CableEquation for functions whose defect will be written to by this channel
		std::vector<size_t> m_vWFctInd;

		/// vector of subsets this channel is declared on
		std::vector<std::string> m_vSubset;

		/// vector of subset indices this channel is declared on (sorted)
		std::vector<int> m_vSI;

		/// joint CableEquation
		CableEquation<TDomain>* m_pCE;

};


} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__CABLE_MEMBRANE_TRANSPORT_INTERFACE_H
