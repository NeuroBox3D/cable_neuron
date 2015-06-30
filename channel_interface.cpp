/*
 * channel_interface.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann
 */

#include "channel_interface.h"


namespace ug {
namespace cable {

template <typename TDomain>
IChannel<TDomain>::IChannel(const char* functions, const char* subsets)
: m_pVMDisc(NULL)
{
	m_vSubset = TokenizeString(subsets);
	m_vWFct = TokenizeString(functions);

	// remove white space
	for (size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);
	for (size_t i = 0; i < m_vWFct.size(); ++i)
		RemoveWhitespaceFromString(m_vWFct[i]);

	// if no function/subsets passed, clear functions
	if (m_vWFct.size() == 1 && m_vWFct[0].empty()) m_vWFct.clear();
	if (m_vSubset.size() == 1 && m_vSubset[0].empty()) m_vSubset.clear();

	// if functions/subsets passed with separator, but not all tokens filled, throw error
	for (size_t i = 0; i < m_vWFct.size(); ++i)
	{
		if (m_vWFct.empty())
		{
			UG_THROW("Error while setting functions in IChannel: Passed function string lacks\n"
					 "a function specification at position " << i << " (of " << m_vWFct.size()-1 << ")");
		}
	}
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
	m_vWFct = functions;
};


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
