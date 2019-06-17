/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-06-12
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


#ifndef UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__KDR_GOLDING01_H
#define UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__KDR_GOLDING01_H

#include "cable_membrane_transport_interface.h"


namespace ug {
namespace cable_neuron {


// forward declaration 
template <typename TDomain> 
class CableEquation; 


/**
 * @brief K delayed rectifier channel from Klee, Ficker and Heinemann
 * modified to account for Dax et al.
 * taken from Golding 2001
 */
template <typename TDomain>
class KDR_Golding01
: public ICableMembraneTransport<TDomain>
{ 
	public:
		using ICableMembraneTransport <TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*, const char*)
		KDR_Golding01(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&, const std::vector<std::string>&)
		KDR_Golding01
		(
			const std::vector<std::string>& functions,
			const std::vector<std::string>& subsets
		);

		/// destructor
		virtual ~KDR_Golding01();

		/// name
		virtual std::string name();


		// setters
		void set_gkbar(number val);
		void set_ikmax(number val);
		void set_vhalfn(number val);
		void set_a0n(number val);
		void set_zetan(number val);
		void set_gmn(number val);
		void set_nmin(number val);
		void set_q10(number val);

#ifdef UG_FOR_LUA
		// function setter for conductance
		void set_gkbar(SmartPtr<LuaUserData<number, TDomain::dim> > fct);
		void set_gkbar(const char* fct);
#endif

		void set_logNGate(bool bLognGate);


		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		virtual std::vector<number> state_values(number x, number y, number z) const;

	private:
		virtual void specify_write_function_indices();

	protected:
		/// create attachments and accessors
		void init_attachments();


	private:
		ANumber m_aNGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaNGate;

#ifdef UG_FOR_LUA
		bool m_bConductanceDependsOnCoordinates;
		SmartPtr<LuaUserData<number, TDomain::dim> > m_spCondFct;
		ANumber m_aGKBar;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaGKBar;
#endif
		number m_gkbar;  // S/m^2
		number m_ikmax;  // A/m^2

		number m_vhalfn;  // V
		number m_a0n;  // s^-1
		number m_zetan;  // 1
		number m_gmn;  // 1
		number m_nmin;  // s
		number m_q10;  // 1

		bool m_bLogNGate;
}; 

} // namespace cable_neuron
} // namespace ug


#endif  // UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__KDR_GOLDING01_H
