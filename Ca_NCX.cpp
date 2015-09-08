/*
 * leakage.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann, mbreit
 */

#include "Ca_NCX.h"


namespace ug {
namespace cable {


template<typename TDomain>
Ca_NCX<TDomain>::Ca_NCX(const char* functions, const char* subsets)
try : IChannel<TDomain>(functions, subsets),
KD_N(1.8e-3), IMAX_N(3.75e-11) {}
UG_CATCH_THROW("Error in Ca_NCX initializer list.");

template<typename TDomain>
Ca_NCX<TDomain>::Ca_NCX
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : IChannel<TDomain>(functions, subsets),
KD_N(1.8e-3), IMAX_N(3.75e-11) {}
UG_CATCH_THROW("Error in Ca_NCX initializer list.");


template<typename TDomain>
std::string Ca_NCX<TDomain>::
name()
{
	return std::string("Ca_NCX");
}

template<typename TDomain>
void Ca_NCX<TDomain>::
set_KD_N(number KD)
{
	KD_N = KD;
}

template<typename TDomain>
void Ca_NCX<TDomain>::
set_IMAX_N(number IMAX)
{
	IMAX_N = IMAX;
}


template<typename TDomain>
void Ca_NCX<TDomain>::vm_disc_available()
{

}


template<typename TDomain>
void Ca_NCX<TDomain>::init_attachments()
{

}


// Methods for using gatings
template<typename TDomain>
void Ca_NCX<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{

}

template<typename TDomain>
void Ca_NCX<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{

}


template<typename TDomain>
void Ca_NCX<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	number ca = vrt_values[VMDisc<TDomain>::_ca_];

	//outCurrentValues.push_back(0);
	// Ca flux in mol/(m^2*ms)
	outCurrentValues.push_back(ca / (KD_N + ca) * IMAX_N);
}



template<typename TDomain>
number Ca_NCX<TDomain>::
lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values)
{
	return 0;
}


template<typename TDomain>
void Ca_NCX<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing VMDisc fct indices which this channel writes to
	//this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_);
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_ca_);
}



////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class Ca_NCX<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class Ca_NCX<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class Ca_NCX<Domain3d>;
#endif


} // namespace cable
} // namespace ug