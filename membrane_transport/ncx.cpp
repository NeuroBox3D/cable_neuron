/*
 * leakage.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann, mbreit
 */

#include "ncx.h"


namespace ug {
namespace cable_neuron {


template<typename TDomain>
NCX_cable<TDomain>::NCX_cable(const char* functions, const char* subsets)
try : ICableMembraneTransport<TDomain>(functions, subsets),
m_kd(1.8e-3), m_max_flux(3.75e-8) {}
UG_CATCH_THROW("Error in NCX_cable initializer list.");

template<typename TDomain>
NCX_cable<TDomain>::NCX_cable
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : ICableMembraneTransport<TDomain>(functions, subsets),
m_kd(1.8e-3), m_max_flux(3.75e-8) {}
UG_CATCH_THROW("Error in NCX_cable initializer list.");


template<typename TDomain>
std::string NCX_cable<TDomain>::
name()
{
	return std::string("NCX_cable");
}

template<typename TDomain>
void NCX_cable<TDomain>::
set_kd(number kd)
{
	m_kd = kd;
}

template<typename TDomain>
void NCX_cable<TDomain>::
set_max_flux(number maxFlux)
{
	m_max_flux = maxFlux;
}


template<typename TDomain>
void NCX_cable<TDomain>::ce_obj_available()
{

}


template<typename TDomain>
void NCX_cable<TDomain>::init_attachments()
{

}


// Methods for using gatings
template<typename TDomain>
void NCX_cable<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{

}

template<typename TDomain>
void NCX_cable<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{

}


template<typename TDomain>
void NCX_cable<TDomain>::current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	number ca = vrt_values[CableEquation<TDomain>::_ca_];

	// Ca flux in mol/(m^2*s)
	outCurrentValues.push_back(ca / (m_kd + ca) * m_max_flux);
}



template<typename TDomain>
number NCX_cable<TDomain>::
lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values)
{
	return 0;
}


template<typename TDomain>
void NCX_cable<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing CableEquation fct indices which this channel writes to
	//this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_);
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_ca_);
}



////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class NCX_cable<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class NCX_cable<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class NCX_cable<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
