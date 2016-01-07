/*
 * leakage.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann, mbreit
 */

#include "Ca_PMCA.h"


namespace ug {
namespace cable {


template<typename TDomain>
Ca_PMCA<TDomain>::Ca_PMCA(const char* functions, const char* subsets)
try : IChannel<TDomain>(functions, subsets),
KD_P(6.0e-5), IMAX_P(8.5e-12) {}
UG_CATCH_THROW("Error in Ca_PMCA initializer list.");

template<typename TDomain>
Ca_PMCA<TDomain>::Ca_PMCA
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : IChannel<TDomain>(functions, subsets),
KD_P(6.0e-5), IMAX_P(8.5e-12) {}
UG_CATCH_THROW("Error in Ca_PMCA initializer list.");


template<typename TDomain>
std::string Ca_PMCA<TDomain>::
name()
{
	return std::string("Ca_PMCA");
}

template<typename TDomain>
void Ca_PMCA<TDomain>::
set_KD_P(number KD)
{
	KD_P = KD;
}

template<typename TDomain>
void Ca_PMCA<TDomain>::
set_IMAX_P(number IMAX)
{
	IMAX_P = IMAX;
}



template<typename TDomain>
void Ca_PMCA<TDomain>::vm_disc_available()
{

}


template<typename TDomain>
void Ca_PMCA<TDomain>::init_attachments()
{

}


// Methods for using gatings
template<typename TDomain>
void Ca_PMCA<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{

}

template<typename TDomain>
void Ca_PMCA<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{

}


template<typename TDomain>
void Ca_PMCA<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	number ca = vrt_values[CableEquation<TDomain>::_ca_];

	//outCurrentValues.push_back(0);	// implement!?
	outCurrentValues.push_back(ca*ca / (KD_P*KD_P + ca*ca) * IMAX_P); // mol/(m^2*ms)
}



template<typename TDomain>
number Ca_PMCA<TDomain>::
lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values)
{
	return 0;
}


template<typename TDomain>
void Ca_PMCA<TDomain>::
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
	template class Ca_PMCA<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class Ca_PMCA<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class Ca_PMCA<Domain3d>;
#endif


} // namespace cable
} // namespace ug
