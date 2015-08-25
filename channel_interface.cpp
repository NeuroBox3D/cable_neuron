/*
 * channel_interface.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann, mbreit
 */

#include "channel_interface.h"


namespace ug {
namespace cable {

template <typename TDomain>
IChannel<TDomain>::IChannel(const char* functions, const char* subsets)
: m_pVMDisc(NULL)
{
	m_vSubset = TokenizeString(subsets);

	// remove white space
	for (size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);

	// if no subsets passed, clear
	if (m_vSubset.size() == 1 && m_vSubset[0].empty()) m_vSubset.clear();

	// if subsets passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < m_vSubset.size(); ++i)
	{
		if (m_vSubset.empty())
		{
			UG_THROW("Error while setting subsets in IChannel: Passed subset string lacks\n"
					 "a subset specification at position " << i << " (of " << m_vSubset.size()-1 << ")");
		}
	}
};


template <typename TDomain>
IChannel<TDomain>::IChannel(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
: m_pVMDisc(NULL)
{
	m_vSubset = subsets;
};


template <typename TDomain>
void IChannel<TDomain>::approx_space_available()
{
	ConstSmartPtr<ApproximationSpace<TDomain> > approx = m_pVMDisc->approx_space();

	// get indices for subsets
	ConstSmartPtr<MGSubsetHandler> ssh = approx->domain()->subset_handler();

	size_t sz = m_vSubset.size();
	m_vSI.resize(sz);
	for (size_t i = 0; i < sz; ++i)
		m_vSI[i] = ssh->get_subset_index(m_vSubset[i].c_str());

	// sort for faster access
	std::sort(m_vSI.begin(), m_vSI.end());

	vm_disc_available();

	// we want to do this _after_ vm_disc_available()
	// as some implementations (e.g. IonLeakage) need this
	specify_write_function_indices();
}


template <typename TDomain>
bool IChannel<TDomain>::
is_def_on_subset(int si) const
{
	size_t sz = m_vSI.size();
	size_t i = 0;
	for (; i < sz && m_vSI[i] < si; ++i)
		;
	if (i == sz) return false;
	if (m_vSI[i] == si) return true;
	return false;
}


//	explicit template instantiations //

#ifdef UG_DIM_1
	template class IChannel<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class IChannel<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class IChannel<Domain3d>;
#endif

} // namespace cable
} // namespace ug
