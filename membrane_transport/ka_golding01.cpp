/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-06-12
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

#include <cmath>
#include <limits>

#include "ka_golding01.h"

namespace ug { 
namespace cable_neuron { 


template <typename TDomain>
KA_Golding01<TDomain>::
KA_Golding01(const char* functions, const char* subsets)
try :
	ICableMembraneTransport<TDomain>(functions, subsets),
#ifdef UG_FOR_LUA
	m_bConductanceDependsOnCoordinates(false),
	m_spCondFct(SPNULL),
#endif
	m_gkbar(80.0),
	m_vhalfn(0.011),
	m_vhalfnDist(-0.001),
	m_vhalfl(-0.056),
	m_a0n(50.0),
	m_a0nDist(100.0),
	m_zetan(-1.5),
	m_zetanDist(-1.8),
	m_zetal(3.0),
	m_gmn(0.55),
	m_gmnDist(0.39),
	m_gml(1.0),
	m_lmin(2e-3),
	m_nmin(1e-4),
	m_pw(-1.0),
	m_tq(-0.04),
	m_qq(5e-3),
	m_q10(5.0),
	m_qtl(1.0),
#ifdef UG_FOR_LUA
	m_bDistinguishDistalRegions(false),
	m_spProximalityFct(SPNULL),
#endif
	m_bLogNGate(false),
	m_bLogLGate(false)
{}
UG_CATCH_THROW("Error in KA_Golding01 initializer list.");


template <typename TDomain>
KA_Golding01<TDomain>::
KA_Golding01
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try :
	ICableMembraneTransport<TDomain>(functions, subsets),
#ifdef UG_FOR_LUA
	m_bConductanceDependsOnCoordinates(false),
	m_spCondFct(SPNULL),
#endif
	m_gkbar(80.0),
	m_vhalfn(0.011),
	m_vhalfnDist(-0.001),
	m_vhalfl(-0.056),
	m_a0n(50.0),
	m_a0nDist(100.0),
	m_zetan(-1.5),
	m_zetanDist(-1.8),
	m_zetal(3.0),
	m_gmn(0.55),
	m_gmnDist(0.39),
	m_gml(1.0),
	m_lmin(2e-3),
	m_nmin(1e-4),
	m_pw(-1.0),
	m_tq(-0.04),
	m_qq(5e-3),
	m_q10(5.0),
	m_qtl(1.0),
#ifdef UG_FOR_LUA
	m_bDistinguishDistalRegions(false),
	m_spProximalityFct(SPNULL),
#endif
	m_bLogNGate(false),
	m_bLogLGate(false)
{}
UG_CATCH_THROW("Error in KA_Golding01 initializer list.");


template <typename TDomain>
KA_Golding01<TDomain>::
~KA_Golding01()
{
#ifdef UG_FOR_LUA
	if (m_bConductanceDependsOnCoordinates)
	{
		SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();
		if (spGrid->has_vertex_attachment(m_aGKBar))
			spGrid->detach_from_vertices(m_aGKBar);
	}
	if (m_bDistinguishDistalRegions)
	{
		SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();
		if (spGrid->has_vertex_attachment(m_aProx))
			spGrid->detach_from_vertices(m_aProx);
	}
#endif
}


template <typename TDomain>
std::string KA_Golding01<TDomain>::
name()
{
	return std::string("KA_Golding01");
}





template <typename TDomain>
void KA_Golding01<TDomain>::set_gkbar(number val)
{ 
	m_gkbar = val;
}

#ifdef UG_FOR_LUA
template <typename TDomain>
void KA_Golding01<TDomain>::set_gkbar(SmartPtr<LuaUserData<number, TDomain::dim> > fct)
{
	m_spCondFct = fct;
	m_bConductanceDependsOnCoordinates = true;
}

template <typename TDomain>
void KA_Golding01<TDomain>::set_gkbar(const char* fct)
{
	m_spCondFct = LuaUserDataFactory<number, TDomain::dim>::create(fct);
	m_bConductanceDependsOnCoordinates = true;
}
#endif


template <typename TDomain>
void KA_Golding01<TDomain>::set_vhalfn(number val)
{ 
	m_vhalfn = val;
}

template <typename TDomain>
void KA_Golding01<TDomain>::set_a0n(number val)
{ 
	m_a0n = val;
}

template <typename TDomain>
void KA_Golding01<TDomain>::set_zetan(number val)
{ 
	m_zetan = val;
}

template <typename TDomain>
void KA_Golding01<TDomain>::set_gmn(number val)
{ 
	m_gmn = val;
}


template <typename TDomain>
void KA_Golding01<TDomain>::set_vhalfn_dist(number val)
{
	m_vhalfnDist = val;
}

template <typename TDomain>
void KA_Golding01<TDomain>::set_a0n_dist(number val)
{
	m_a0nDist = val;
}

template <typename TDomain>
void KA_Golding01<TDomain>::set_zetan_dist(number val)
{
	m_zetanDist = val;
}

template <typename TDomain>
void KA_Golding01<TDomain>::set_gmn_dist(number val)
{
	m_gmnDist = val;
}


#ifdef UG_FOR_LUA
template <typename TDomain>
void KA_Golding01<TDomain>::set_proximality_fct(SmartPtr<LuaUserData<number, TDomain::dim> > fct)
{
	m_spProximalityFct = fct;
	m_bDistinguishDistalRegions = true;
}

template <typename TDomain>
void KA_Golding01<TDomain>::set_proximality_fct(const char* fct)
{
	m_spProximalityFct = LuaUserDataFactory<number, TDomain::dim>::create(fct);
	m_bDistinguishDistalRegions = true;
}
#endif


template <typename TDomain>
void KA_Golding01<TDomain>::
set_logLGate(bool bLoglGate)
{
	m_bLogLGate = bLoglGate;
}

template <typename TDomain>
void KA_Golding01<TDomain>::
set_logNGate(bool bLognGate)
{
	m_bLogNGate = bLognGate;
}




template <typename TDomain>
void KA_Golding01<TDomain>::
init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	// init gating values
	const number F = m_pCE->F;
	const number R = m_pCE->R;
	const number T = m_pCE->temperature();
	const number v = vrt_values[CableEquation<TDomain>::_v_];

	number vhalfn = m_vhalfn;
	number zetan = m_zetan;
#ifdef UG_FOR_LUA
	if (m_bDistinguishDistalRegions)
	{
		if (!m_aaProx[vrt])
		{
			vhalfn = m_vhalfnDist;
			zetan = m_zetanDist;
		}
	}
#endif

	const number zeta = zetan + m_pw / (1.0 + exp((v-m_tq) / m_qq));
	number a = exp(zeta * (v - vhalfn) * F/(R*T));
	m_aaNGate[vrt] = 1.0 / (1.0 + a);

	a = exp(m_zetal * (v - m_vhalfl) * F/(R*T));
	m_aaLGate[vrt] = 1.0 / (1.0 + a);

#ifdef UG_FOR_LUA
	// write conductance value to attachment
	if (m_bConductanceDependsOnCoordinates)
	{
		const int si = m_pCE->subset_handler().get_subset_index(vrt);
		const MathVector<TDomain::dim>& coords =
			this->m_pCE->approx_space()->domain()->position_accessor()[vrt];
		m_spCondFct->evaluate(m_aaGKBar[vrt], coords, 0.0, si);
	}

	// write proximality information to attachment
	if (m_bDistinguishDistalRegions)
	{
		const int si = m_pCE->subset_handler().get_subset_index(vrt);
		const MathVector<TDomain::dim>& coords =
			this->m_pCE->approx_space()->domain()->position_accessor()[vrt];
		number tmp = 1.0;
		m_spProximalityFct->evaluate(tmp, coords, 0.0, si);
		m_aaProx[vrt] = tmp != 0.0;
	}
#endif
}


template <typename TDomain>
void KA_Golding01<TDomain>::
update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{ 
	const number F = m_pCE->F;
	const number R = m_pCE->R;
	const number T = m_pCE->temperature();
	const number celsius = m_pCE->temperature_celsius();

	const number dt = newTime - m_pCE->time();
	const number v = vrt_values[CableEquation<TDomain>::_v_];
	number& n = m_aaNGate[vrt];
	number& l = m_aaLGate[vrt];

	// decide whether this is a proximal or a distal vertex
	number vhalfn = m_vhalfn;
	number zetan = m_zetan;
	number gmn = m_gmn;
	number a0n = m_a0n;
#ifdef UG_FOR_LUA
	if (m_bDistinguishDistalRegions)
	{
		if (!m_aaProx[vrt])
		{
			vhalfn = m_vhalfnDist;
			zetan = m_zetanDist;
			gmn = m_gmnDist;
			a0n = m_a0nDist;
		}
	}
#endif

	const number zeta = zetan + m_pw / (1.0 + exp((v-m_tq) / m_qq));
	number a = exp(zeta * (v - vhalfn) * F/(R*T));
	const number ninf = 1.0 / (1.0 + a);

	const number b = exp(zeta * gmn * (v - vhalfn) * F/(R*T));
	const number qt = pow(m_q10, (celsius - 24.0) / 10.0);
	const number taun = std::max(m_nmin, b / (qt * a0n * (1.0 + a)));

	// we give the exact solution (instead of implicit Euler)
	n = ninf + (n - ninf) * exp(-dt/taun);


	a = exp(m_zetal * (v - m_vhalfl) * F/(R*T));
	const number linf = 1.0 / (1.0 + a);
	const number taul = std::max(m_lmin/m_qtl, 0.26 * (v + 0.05) / m_qtl);

	// we give the exact solution (instead of implicit Euler)
	l = linf + (l - linf) * exp(-dt/taul);
}


template <typename TDomain>
void KA_Golding01<TDomain>::current
(
	Vertex* vrt,
	const std::vector<number>& vrt_values,
	std::vector<number>& outCurrentValues
)
{
	const number v = vrt_values[CableEquation<TDomain>::_v_];
	const number ek = this->m_pCE->rev_pot_k();

	number gkbar = m_gkbar;
#ifdef UG_FOR_LUA
	if (m_bConductanceDependsOnCoordinates)
		gkbar = m_aaGKBar[vrt];
#endif
	outCurrentValues.push_back(gkbar * m_aaNGate[vrt] * m_aaLGate[vrt] * (v - ek));
}


template <typename TDomain>
void KA_Golding01<TDomain>::ce_obj_available()
{
	// init the attachments
	init_attachments();
}


template <typename TDomain>
std::vector<number> KA_Golding01<TDomain>::state_values(number x, number y, number z) const
{ 
#ifdef UG_PARALLLEL
	if (pcl::NumProcs() > 1)
		UG_THROW("The state_values method is not parallelized.");
#endif

	std::vector<number> states;

	if (!(m_bLogNGate || m_bLogLGate))
		return states;

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	typedef ug::MathVector<TDomain::dim> position_type;

	position_type coord;
	if (coord.size() >= 1)
		coord[0] = x;
	if (coord.size() >= 2)
		coord[1] = y;
	if (coord.size() >= 3)
		coord[2] = z;

	// get position accessor
	const typename TDomain::position_accessor_type& aaPos =
		m_pCE->approx_space()->domain()->position_accessor();

	// get surface dof distribution
	ConstSmartPtr<DoFDistribution> dd = m_pCE->approx_space()->dof_distribution(GridLevel());

	SubsetGroup ssGrp;
	try {ssGrp = SubsetGroup(m_pCE->approx_space()->domain()->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");
	const size_t nSs = ssGrp.size();

	// iterating over all vertices to find the one closest to desired coords
	number bestDistSq = std::numeric_limits<number>::max();
	Vertex* bestVrt = NULL;
	for (size_t si = 0; si < nSs; ++si)
	{
		itType it = dd->template begin<Vertex>(ssGrp[si]);
		itType itEnd = dd->template end<Vertex>(ssGrp[si]);

		for (; it != itEnd; ++it)
		{
			const number distSq = VecDistanceSq(coord, aaPos[*it]);
			if (distSq < bestDistSq)
			{
				bestDistSq = distSq;
				bestVrt = *it;
			}
		}
	}

	UG_COND_THROW(!bestVrt, "No vertex found near the desired coordinates, in fact, none found at all.");

	if (m_bLogNGate)
		states.push_back(m_aaNGate[bestVrt]);
	if (m_bLogLGate)
		states.push_back(m_aaLGate[bestVrt]);

	return states;
}


template <typename TDomain>
void KA_Golding01<TDomain>::specify_write_function_indices()
{
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_);
}



template <typename TDomain>
void KA_Golding01<TDomain>::init_attachments()
{
	SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();

	// attach gatings
	if (spGrid->has_vertex_attachment(m_aNGate))
		UG_THROW("Attachment necessary (nGate) for KA_Golding01 channel dynamics "
			"could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_aNGate);
	m_aaNGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aNGate);

	if (spGrid->has_vertex_attachment(this->m_aLGate))
		UG_THROW("Attachment necessary (lGate) for KA_Golding01 channel dynamics "
			"could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_aLGate);
	m_aaLGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aLGate);

#ifdef UG_FOR_LUA
	// attach conductance attachment if necessary
	if (m_bConductanceDependsOnCoordinates)
	{
		if (spGrid->has_vertex_attachment(m_aGKBar))
			UG_THROW("Attachment necessary (gKBar) for " << name() << " channel dynamics "
				"could not be made, since it already exists.");
		spGrid->attach_to_vertices(m_aGKBar);
		m_aaGKBar = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aGKBar);
	}

	// attach proximality attachment if necessary
	if (m_bDistinguishDistalRegions)
	{
		if (spGrid->has_vertex_attachment(m_aProx))
			UG_THROW("Attachment necessary (prox) for " << name() << " channel dynamics "
				"could not be made, since it already exists.");
		spGrid->attach_to_vertices(m_aProx);
		m_aaProx = Grid::AttachmentAccessor<Vertex, ABool>(*spGrid, m_aProx);
	}
#endif
}



//	explicit template instantiations
#ifdef UG_DIM_3
template class KA_Golding01<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
