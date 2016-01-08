/*
 * leakage.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann, mbreit
 */

#include "leakage.h"


namespace ug {
namespace cable_neuron {

template<typename TDomain>
ChannelLeak<TDomain>::ChannelLeak(const char* functions, const char* subsets)
try : ICableMembraneTransport<TDomain>(functions, subsets) {}
UG_CATCH_THROW("Error in ChannelHH initializer list.");

template<typename TDomain>
ChannelLeak<TDomain>::ChannelLeak
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : ICableMembraneTransport<TDomain>(functions, subsets) {}
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
	set_cond(g, this->m_vSubset);
}
template<typename TDomain>
void ChannelLeak<TDomain>::
set_cond(number g, const char* subsets)
{
	std::vector<std::string> ssVec;
	try	{this->subsetCString2Vector(ssVec, subsets);}
	UG_CATCH_THROW("Error while setting conductance value for leakage term.");

	set_cond(g, ssVec);
}
template<typename TDomain>
void ChannelLeak<TDomain>::
set_cond(number g, const std::vector<std::string>& subsets)
{
	size_t sz = subsets.size();
	for (size_t i = 0; i < sz; ++i)
		m_mSubsetParams2Save[subsets[i]].g = g;
}

template<typename TDomain>
void ChannelLeak<TDomain>::set_rev_pot(number e)
{
	set_rev_pot(e, this->m_vSubset);
}
template<typename TDomain>
void ChannelLeak<TDomain>::set_rev_pot(number e, const char* subsets)
{
	std::vector<std::string> ssVec;
	try	{this->subsetCString2Vector(ssVec, subsets);}
	UG_CATCH_THROW("Error while setting reversal potential value for leakage term.");

	set_rev_pot(e, ssVec);
}
template<typename TDomain>
void ChannelLeak<TDomain>::set_rev_pot(number e, const std::vector<std::string>& subsets)
{
	size_t sz = subsets.size();
	for (size_t i = 0; i < sz; ++i)
		m_mSubsetParams2Save[subsets[i]].E = e;
}


template<typename TDomain>
void ChannelLeak<TDomain>::ce_obj_available()
{
// save parameters for subset indices
	ConstSmartPtr<MGSubsetHandler> ssh = m_pCE->approx_space()->domain()->subset_handler();

	// if special params saved for individual subsets, take these
	typedef typename std::map<std::string, Params>::const_iterator MapIter;
	MapIter it = m_mSubsetParams2Save.begin();
	MapIter itEnd = m_mSubsetParams2Save.end();
	for (; it != itEnd; ++it)
	{
		int si = ssh->get_subset_index(it->first.c_str());
		if (si == -1)
		{
			UG_THROW("Unknown subset '" << it->first << "' in '"
					 << name() << "' channel mechanism.");
		}
		m_mSubsetParams[si].g = it->second.g;
		m_mSubsetParams[si].E = it->second.E;
	}

	// for the subsets this channel is defined on, but where no individual
	// parameterization is given, take default params, but warn
	size_t sz = this->m_vSI.size();
	for (size_t i = 0; i < sz; ++i)
	{
		if (m_mSubsetParams.find(this->m_vSI[i]) == m_mSubsetParams.end())

		// warn
		UG_LOG_ALL_PROCS("WARNING: Leakage defined on subset '" << ssh->get_subset_name(this->m_vSI[i])
				<< "', but no parameterization provided there. Taking default parameters." << std::endl);

		// set default params
		m_mSubsetParams[this->m_vSI[i]];
	}

	m_mSubsetParams2Save.clear();
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
void ChannelLeak<TDomain>::current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	number VM 	 = vrt_values[m_pCE->_v_];

	// params for this subset
	int si = m_pCE->current_subset_index();
	const number g = m_mSubsetParams[si].g;
	const number E = m_mSubsetParams[si].E;

	const number leakage_part_of_flux = g * (VM - E);

	outCurrentValues.push_back(leakage_part_of_flux);
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


template<typename TDomain>
number ChannelLeak<TDomain>::
lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values)
{
	return m_mSubsetParams[m_pCE->current_subset_index()].g;
}


template<typename TDomain>
void ChannelLeak<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing CableEquation fct indices which this channel writes to
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_);
}


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


} // namespace cable_neuron
} // namespace ug
