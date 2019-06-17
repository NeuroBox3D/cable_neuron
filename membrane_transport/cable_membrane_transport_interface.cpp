/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Pascal Gottmann
 * Creation date: 2014-10-29
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "cable_membrane_transport_interface.h"

#include <algorithm>


namespace ug {
namespace cable_neuron {

template <typename TDomain>
ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char* functions, const char* subsets)
: m_pCE(NULL)
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
			UG_THROW("Error while setting subsets in ICableMembraneTransport: Passed subset string lacks\n"
					 "a subset specification at position " << i << " (of " << m_vSubset.size()-1 << ")");
		}
	}
};


template <typename TDomain>
ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
: m_pCE(NULL)
{
	m_vSubset = subsets;
};


template <typename TDomain>
void ICableMembraneTransport<TDomain>::approx_space_available()
{
	ConstSmartPtr<ApproximationSpace<TDomain> > approx = m_pCE->approx_space();

	// get indices for subsets
	ConstSmartPtr<MGSubsetHandler> ssh = approx->domain()->subset_handler();

	size_t sz = m_vSubset.size();
	m_vSI.resize(sz);
	for (size_t i = 0; i < sz; ++i)
	{
		m_vSI[i] = ssh->get_subset_index(m_vSubset[i].c_str());
		if (m_vSI[i] == -1)
		{
			UG_THROW("Unknown subset '" << m_vSubset[i] << "' in '"
					 << name() << "' channel mechanism.");
		}
	}
	// sort for faster access
	std::sort(m_vSI.begin(), m_vSI.end());

	ce_obj_available();

	// we want to do this _after_ ce_obj_available()
	// as some implementations (e.g. IonLeakage) need this
	specify_write_function_indices();
}


template <typename TDomain>
number ICableMembraneTransport<TDomain>::lin_dep_on_pot
(
	Vertex* vrt,
	const std::vector<number>& vrt_values
)
{
	// for easy access to membrane potential index
	size_t vi = CableEquation<TDomain>::_v_;

	// find index in current() output vector corresponding to membrane potential
	size_t i = 0;
	size_t sz = m_vWFctInd.size();
	while (i < sz)
	{
		if (m_vWFctInd[i++] == vi)
		{
			--i;
			break;
		}
	}

	// no influence on membrane potential; return 0
	if (i == sz) return 0.0;

	// call current twice
	std::vector<number> currentDensity;
	std::vector<number> vrt_vals2(vrt_values);
	vrt_vals2[vi] *= 1.0+1e-8;
	number dv = vrt_vals2[vi] - vrt_values[vi];

	current(vrt, vrt_vals2, currentDensity);
	number dcd = currentDensity[i];
	currentDensity.clear();
	current(vrt, vrt_values, currentDensity);
	dcd -= currentDensity[i];

	return dcd / dv;
}


template <typename TDomain>
bool ICableMembraneTransport<TDomain>::
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


template<typename TDomain>
void ICableMembraneTransport<TDomain>::
subsetCString2Vector(std::vector<std::string>& outVec, const char* cstr)
{
	// tokenize
	outVec = TokenizeString(cstr);

	// remove white space
	for (size_t i = 0; i < outVec.size(); ++i)
		RemoveWhitespaceFromString(outVec[i]);

	// if no subsets passed, clear
	if (outVec.size() == 1 && outVec[0].empty()) outVec.clear();

	// if subsets passed with separator, but not all tokens filled, throw error
	for (size_t i = 0; i < outVec.size(); ++i)
	{
		if (outVec.empty())
		{
			UG_THROW("Passed subset string lacks a subset specification at position "
					 << i << " (of " << outVec.size()-1 << ")");
		}
	}
}


template <typename TDomain>
void ICableMembraneTransport<TDomain>::
subsetNames2Indices(std::vector<int>& ind, const std::vector<std::string>& names)
{
	//
}


//	explicit template instantiations
#ifdef UG_DIM_3
	template class ICableMembraneTransport<Domain3d>;
#endif

} // namespace cable_neuron
} // namespace ug
