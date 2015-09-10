/*
 * Na_K_Pump.cpp
 *
 *  Created on: 09.09.2015
 *      Author: pgottmann, mbreit
 *      Out of Paper: Tissue-specific Versus Isoform-specific Differences in Cation
						Activation Kinetics of the Na,K-ATPase
						Therien et all (2006)
 */

#include "Na_K_Pump.h"


namespace ug {
namespace cable {


template<typename TDomain>
Na_K_Pump<TDomain>::Na_K_Pump(const char* functions, const char* subsets)
try : IChannel<TDomain>(functions, subsets),
K_K(1.37), K_Na(5.74), IMAX_P(3.6e-5) {}
UG_CATCH_THROW("Error in Na_K_Pump initializer list.");

template<typename TDomain>
Na_K_Pump<TDomain>::Na_K_Pump
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : IChannel<TDomain>(functions, subsets),
K_K(1.37), K_Na(5.74), IMAX_P(3.6e-5) {}
UG_CATCH_THROW("Error in Na_K_Pump initializer list.");


template<typename TDomain>
std::string Na_K_Pump<TDomain>::
name()
{
	return std::string("Na_K_Pump");
}

template<typename TDomain>
void Na_K_Pump<TDomain>::
set_K_K(number K)
{
	K_K = K;
}

template<typename TDomain>
void Na_K_Pump<TDomain>::
set_K_Na(number Na)
{
	K_Na = Na;
}

template<typename TDomain>
void Na_K_Pump<TDomain>::
set_IMAX_P(number IMAX)
{
	IMAX_P = IMAX;
}



template<typename TDomain>
void Na_K_Pump<TDomain>::vm_disc_available()
{

}


template<typename TDomain>
void Na_K_Pump<TDomain>::init_attachments()
{

}


// Methods for using gatings
template<typename TDomain>
void Na_K_Pump<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{

}

template<typename TDomain>
void Na_K_Pump<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{

}


template<typename TDomain>
void Na_K_Pump<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	//number ca = vrt_values[VMDisc<TDomain>::_ca_];
	number na = vrt_values[VMDisc<TDomain>::_na_];
	number k = vrt_values[VMDisc<TDomain>::_k_];

	number napump = 1.0 /(1.0 + K_Na/na * (1.0 + k/K_K));
	napump *= napump*napump;
	napump *= IMAX_P;

	//std::cout << "pumping: " << napump << std::endl;
	//std::cout << "pumping: " << ((-0.66666666)*napump) << std::endl;
	//outCurrentValues.push_back(0);	// implement!?
	outCurrentValues.push_back(napump); // mol/(m^2*ms)
	// 3na vs 2k
	outCurrentValues.push_back((-2.0/3.0) * napump); // mol/(m^2*ms)
}



template<typename TDomain>
number Na_K_Pump<TDomain>::
lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values)
{
	return 0;
}


template<typename TDomain>
void Na_K_Pump<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing VMDisc fct indices which this channel writes to
	//this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_);
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_na_);
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_k_);
}



////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class Na_K_Pump<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class Na_K_Pump<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class Na_K_Pump<Domain3d>;
#endif


} // namespace cable
} // namespace ug
