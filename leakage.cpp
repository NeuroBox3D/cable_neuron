/*
 * leakage.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann
 */

#include "leakage.h"


namespace ug {
namespace cable {

template<typename TDomain>
ChannelLeak<TDomain>::ChannelLeak(const char* functions, const char* subsets)
try : IChannel<TDomain>(functions, subsets),
m_g(1.0e-6), m_E(-65.0) {}
UG_CATCH_THROW("Error in ChannelHH initializer list.");

template<typename TDomain>
ChannelLeak<TDomain>::ChannelLeak
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : IChannel<TDomain>(functions, subsets),
m_g(1.0e-6), m_E(-65.0) {}
UG_CATCH_THROW("Error in ChannelHH initializer list.");


template<typename TDomain>
std::string ChannelLeak<TDomain>::
name()
{
	return std::string("Leakage");
}


template<typename TDomain>
void ChannelLeak<TDomain>::
set_cond(number g)
{
	m_g = g;
}

template<typename TDomain>
void ChannelLeak<TDomain>::set_rev_pot(number e)
{
	m_E = e;
}


template<typename TDomain>
void ChannelLeak<TDomain>::vm_disc_available()
{

}


template<typename TDomain>
void ChannelLeak<TDomain>::init_attachments()
{

}

template<typename TDomain>
std::vector<number> ChannelLeak<TDomain>::state_values(number x, number y, number z)
{
	std::vector<number> GatingAccesors;


	return GatingAccesors;
}

// Methods for using gatings
template<typename TDomain>
void ChannelLeak<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	// nothing to do
}

template<typename TDomain>
void ChannelLeak<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{
	// nothing to do
}


template<typename TDomain>
void ChannelLeak<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	number VM 	 = vrt_values[m_pVMDisc->_v_];

	const number leakage_part_of_flux = m_g * (VM - m_E);

	number flux_value = (leakage_part_of_flux);
	outCurrentValues.push_back(flux_value);
}


#if 0
template<typename TDomain>
void ChannelLeak<TDomain>::Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux)
{
	number NGate = m_aaNGate[vrt];
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];


	number Jac = (m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g);

	outJFlux.push_back(Jac);

}

#endif




////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class ChannelLeak<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class ChannelLeak<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class ChannelLeak<Domain3d>;
#endif


} // namespace cable
} // namespace ug
