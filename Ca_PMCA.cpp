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
KD_P(6.0e-5), IMAX_P(1.7e-23) {}
UG_CATCH_THROW("Error in Ca_PMCA initializer list.");

template<typename TDomain>
Ca_PMCA<TDomain>::Ca_PMCA
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : IChannel<TDomain>(functions, subsets),
KD_P(6.0e-5), IMAX_P(1.7e-23) {}
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
void Ca_PMCA<TDomain>::
set_scaling(number scale)
{
	m_scaling = scale;
}




template<typename TDomain>
void Ca_PMCA<TDomain>::vm_disc_available()
{
	init_attachments();
}


template<typename TDomain>
void Ca_PMCA<TDomain>::init_attachments()
{
	SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid();


	if (spGrid->has_vertex_attachment(this->gatingFactorGate))
	UG_THROW("Attachment necessary (mGate) for ca_converted_allNernst_UG channel dynamics "
	"could not be made, since it already exists.");
	spGrid->attach_to_vertices(this->gatingFactorGate);
	this->aagatingFactorGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->gatingFactorGate);
}

template<typename TDomain>
std::vector<number> Ca_PMCA<TDomain>::state_values(number x, number y, number z)
{
	std::vector<number> GatingAccesors;


	return GatingAccesors;
}

// Methods for using gatings
template<typename TDomain>
void Ca_PMCA<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	number ca = vrt_values[VMDisc<TDomain>::_ca_];



	this->aagatingFactorGate[vrt] = ca*ca / (KD_P*KD_P + ca*ca);

}

template<typename TDomain>
void Ca_PMCA<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{


	number dt = newTime - m_pVMDisc->time();
	number v = vrt_values[VMDisc<TDomain>::_v_];
	number ca = vrt_values[VMDisc<TDomain>::_ca_];

	aagatingFactorGate[vrt] += ((2*KD_P*KD_P*ca / std::pow(KD_P*KD_P + ca*ca, 2))*dt);


}


template<typename TDomain>
void Ca_PMCA<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	number gatingFactor = aagatingFactorGate[vrt];

	outCurrentValues.push_back(0);
	outCurrentValues.push_back(gatingFactor * IMAX_P * m_scaling);
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
	// prepare vector containing VMDisc fct indices which this channel writes to
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_);
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_ca_);
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
