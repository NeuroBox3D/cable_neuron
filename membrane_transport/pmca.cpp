/*
 * leakage.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann, mbreit
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

#ifdef UG_DIM_1
	template class PMCA_cable<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class PMCA_cable<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class PMCA_cable<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
