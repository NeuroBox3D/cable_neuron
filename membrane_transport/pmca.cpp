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

#include "pmca.h"


namespace ug {
namespace cable_neuron {


template<typename TDomain>
PMCA_cable<TDomain>::PMCA_cable(const char* functions, const char* subsets)
try : ICableMembraneTransport<TDomain>(functions, subsets),
m_kd(6.0e-5), m_maxFlux(8.5e-9) {}
UG_CATCH_THROW("Error in PMCA_cable initializer list.");

template<typename TDomain>
PMCA_cable<TDomain>::PMCA_cable
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : ICableMembraneTransport<TDomain>(functions, subsets),
m_kd(6.0e-5), m_maxFlux(8.5e-9) {}
UG_CATCH_THROW("Error in PMCA_cable initializer list.");


template<typename TDomain>
std::string PMCA_cable<TDomain>::
name()
{
	return std::string("PMCA_cable");
}

template<typename TDomain>
void PMCA_cable<TDomain>::
set_kd(number kd)
{
	m_kd = kd;
}

template<typename TDomain>
void PMCA_cable<TDomain>::
set_max_flux(number maxFlux)
{
	m_maxFlux = maxFlux;
}



template<typename TDomain>
void PMCA_cable<TDomain>::ce_obj_available()
{

}


template<typename TDomain>
void PMCA_cable<TDomain>::init_attachments()
{

}


// Methods for using gatings
template<typename TDomain>
void PMCA_cable<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{

}

template<typename TDomain>
void PMCA_cable<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{

}


template<typename TDomain>
void PMCA_cable<TDomain>::current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	number ca = vrt_values[CableEquation<TDomain>::_ca_];

	number ionic_current = ca*ca / (m_kd*m_kd + ca*ca) * m_maxFlux;
	outCurrentValues.push_back(ionic_current * 2*this->m_pCE->F);
	outCurrentValues.push_back(ionic_current); // mol/(m^2*ms)
}



template<typename TDomain>
number PMCA_cable<TDomain>::
lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values)
{
	return 0;
}


template<typename TDomain>
void PMCA_cable<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing CableEquation fct indices which this channel writes to
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_);
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_ca_);
}



////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////
#ifdef UG_DIM_3
	template class PMCA_cable<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
