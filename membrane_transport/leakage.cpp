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

#include "leakage.h"


namespace ug {
namespace cable_neuron {

template<typename TDomain>
ChannelLeak<TDomain>::ChannelLeak(const char* functions, const char* subsets)
try : ICableMembraneTransport<TDomain>(functions, subsets)
#ifdef UG_FOR_LUA
	, m_spCondFct(SPNULL)
#endif
{}
UG_CATCH_THROW("Error in ChannelHH initializer list.");

template<typename TDomain>
ChannelLeak<TDomain>::ChannelLeak
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : ICableMembraneTransport<TDomain>(functions, subsets)
#ifdef UG_FOR_LUA
	, m_spCondFct(SPNULL)
#endif
{}
UG_CATCH_THROW("Error in ChannelHH initializer list.");

template <typename TDomain>
ChannelLeak<TDomain>::~ChannelLeak()
{
#ifdef UG_FOR_LUA
	if (m_spCondFct.valid())
	{
		SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();
		if (spGrid->has_vertex_attachment(m_aG))
			spGrid->detach_from_vertices(m_aG);
	}
#endif
}


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

#ifdef UG_FOR_LUA
template <typename TDomain>
void ChannelLeak<TDomain>::set_cond(SmartPtr<LuaUserData<number, TDomain::dim> > fct)
{
	m_spCondFct = fct;
}

template <typename TDomain>
void ChannelLeak<TDomain>::set_cond(const char* fct)
{
	m_spCondFct = LuaUserDataFactory<number, TDomain::dim>::create(fct);
}
#endif


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
	init_attachments();

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
#ifdef UG_FOR_LUA
	// attach conductance attachment if necessary
	if (m_spCondFct.valid())
	{
		SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();
		if (spGrid->has_vertex_attachment(m_aG))
			UG_THROW("Attachment necessary (g) for " << name() << " channel dynamics "
				"could not be made, since it already exists.");
		spGrid->attach_to_vertices(m_aG);
		m_aaG = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aG);
	}
#endif
}

// Methods for using gatings
template<typename TDomain>
void ChannelLeak<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{
#ifdef UG_FOR_LUA
	// write conductance value to attachment
	if (m_spCondFct.valid())
	{
		const int si = this->m_pCE->approx_space()->domain()->subset_handler()->get_subset_index(vrt);
		const MathVector<TDomain::dim>& coords =
			this->m_pCE->approx_space()->domain()->position_accessor()[vrt];
		m_spCondFct->evaluate(m_aaG[vrt], coords, 0.0, si);
	}
#endif
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
	number VM = vrt_values[m_pCE->_v_];

	// params for this subset
	int si = m_pCE->current_subset_index();
	number g = m_mSubsetParams[si].g;
	const number E = m_mSubsetParams[si].E;

#ifdef UG_FOR_LUA
	if (m_spCondFct.valid())
		g = m_aaG[vrt];
#endif

	const number leakage_part_of_flux = g * (VM - E);

	outCurrentValues.push_back(leakage_part_of_flux);
}


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
#ifdef UG_DIM_3
	template class ChannelLeak<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
