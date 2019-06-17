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

#ifndef UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__LEAKAGE_H
#define UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__LEAKAGE_H


#include "cable_membrane_transport_interface.h"


namespace ug {
namespace cable_neuron {


template <typename TDomain>
class ChannelLeak
	: public ICableMembraneTransport<TDomain>
{
	public:
		using ICableMembraneTransport<TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*)
		ChannelLeak(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&)
		ChannelLeak(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~ChannelLeak();

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// set leakage conductance
		void set_cond(number g);

		/// set leakage conductance on specific subsets
		/// \{
		void set_cond(number g, const char* subsets);
		void set_cond(number g, const std::vector<std::string>& subsets);
		/// \}

#ifdef UG_FOR_LUA
		/// set leakage conductance by Lua function
		/// \{
		void set_cond(SmartPtr<LuaUserData<number, TDomain::dim> > fct);
		void set_cond(const char* fct);
		/// \}
#endif

		/// set leakage equilibrium potential
		void set_rev_pot(number e);
		/// set leakage equilibrium potential on specific subsets
		void set_rev_pot(number e, const char* subsets);
		void set_rev_pot(number e, const std::vector<std::string>& subsets);


		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		struct Params
		{
			Params() : g(1.0), E(-0.065) {};
			number g; ///< membrane conductance (S/m^2)
			number E; ///< reversal potential (V)
		};

#ifdef UG_FOR_LUA
		SmartPtr<LuaUserData<number, TDomain::dim> > m_spCondFct;
		ANumber m_aG;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaG;
#endif

		std::map<std::string, Params> m_mSubsetParams2Save;
		std::map<int, Params> m_mSubsetParams;
};


} // namespace cable_neuron
} // namespace ug

#endif  // UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__LEAKAGE_H
