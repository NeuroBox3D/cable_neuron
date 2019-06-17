/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Pascal Gottmann
 * Creation date: 2015-08-24
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

#ifndef UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__VDCC_BG_H
#define UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__VDCC_BG_H


#include "cable_membrane_transport_interface.h"

namespace ug {
namespace cable_neuron {


template <typename TDomain>
class VDCC_BG_cable
	: public ICableMembraneTransport<TDomain>
{
	public:
		using ICableMembraneTransport<TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*)
		VDCC_BG_cable(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&)
		VDCC_BG_cable(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);;

		/// destructor
		virtual ~VDCC_BG_cable() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// setter for output behavior of gatings
		void set_log_nGate(bool bLogNGate);
		void set_log_hGate(bool bLogHGate);
		void set_log_mGate(bool bLogMGate);

	public:
		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		virtual std::vector<number> state_values(number x, number y, number z) const;
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		// membrane conductivities
		int m_mp;
		number m_z_m;
		number m_v12_m;		// in units of mV
		number m_tau0_m;	// in units of ms

		int m_hp;
		number m_z_h;
		number m_v12_h;		// in units of mV
		number m_tau0_h;	// in units of ms

		// permeability
		number m_perm;		// in units of m/s

		// one attachment per state variable
		ANumber m_MGate;
		ANumber m_HGate;

		Grid::AttachmentAccessor<Vertex, ANumber> m_aaMGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaHGate;

		// attachment log_files
		bool m_log_hGate, m_log_mGate;
};



} // namespace cable_neuron
} // namespace ug

#endif  // UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__VDCC_BG_H
