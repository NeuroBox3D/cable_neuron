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

#ifndef UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__NAX_GOLDING01_H
#define UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__NAX_GOLDING01_H

#include "cable_membrane_transport_interface.h"


namespace ug {
namespace cable_neuron {


// forward declaration 
template <typename TDomain> 
class CableEquation; 


/**
 * @brief Na channel for axons
 * without slow inactivation
 * taken from Golding 2001
 */
template <typename TDomain>
class Nax_Golding01
: public ICableMembraneTransport<TDomain>
{ 
	public:
		using ICableMembraneTransport <TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*, const char*)
		Nax_Golding01(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&, const std::vector<std::string>&)
		Nax_Golding01
		(
			const std::vector<std::string>& functions,
			const std::vector<std::string>& subsets
		);

		/// destructor
		virtual ~Nax_Golding01();

		/// name
		virtual std::string name();


		// setters
		void set_gbar(number val);
		void set_tha(number val);
		void set_qa(number val);
		void set_Ra(number val);
		void set_Rb(number val);
		void set_thi1(number val);
		void set_thi2(number val);
		void set_qd(number val);
		void set_qg(number val);
		void set_mmin(number val);
		void set_hmin(number val);
		void set_q10(number val);
		void set_Rg(number val);
		void set_Rd(number val);
		void set_thinf(number val);
		void set_qinf(number val4);

#ifdef UG_FOR_LUA
		// function setter for conductance
		void set_gbar(SmartPtr<LuaUserData<number, TDomain::dim> > fct);
		void set_gbar(const char* fct);
#endif

		void set_logMGate(bool b);
		void set_logHGate(bool b);


		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		virtual std::vector<number> state_values(number x, number y, number z) const;

	private:
		virtual void specify_write_function_indices();

	protected:
		number alpn(number v);
		number betn(number v);

		/// create attachments and accessors
		void init_attachments();


	private:
		ANumber m_aMGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaMGate;
		ANumber m_aHGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaHGate;

#ifdef UG_FOR_LUA
		bool m_bConductanceDependsOnCoordinates;
		SmartPtr<LuaUserData<number, TDomain::dim> > m_spCondFct;
		ANumber m_aGBar;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaGBar;
#endif
		number m_gbar;  // S/m^2

		number m_tha;  // V
		number m_qa;  // V
		number m_Ra;  // s^-1
		number m_Rb;  // s^-1

		number m_thi1;  // V
		number m_thi2;  // V
		number m_qd;  // V
		number m_qg;  // V

		number m_mmin;  // s
		number m_hmin;  // s
		number m_q10;  // 1
		number m_Rg;  // s^-1
		number m_Rd;  // s^-1

		number m_thinf;  // V
		number m_qinf;  // V

		bool m_bLogMGate;
		bool m_bLogHGate;
}; 

} // namespace cable_neuron
} // namespace ug


#endif  // UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__NAX_GOLDING01_H
