/*
 * kadist_g01.cpp
 *
 *  Created on: 2019-05-17
 *      Author: mbreit
 */

#include <cmath>
#include <limits>
#include "ka_dist_golding01.h"

namespace ug { 
namespace cable_neuron { 


template<typename TDomain>
KaDist_Golding01<TDomain>::
KaDist_Golding01(const char* functions, const char* subsets)
try :
	ICableMembraneTransport<TDomain>(functions, subsets),
	m_ek(-0.09),
	m_gkabar(80.0),
	m_vhalfn(-1e-3),
	m_vhalfl(-0.056),
	m_a0n(100.0),
	m_zetan(-1.8),
	m_zetal(3.0),
	m_gmn(0.39),
	m_gml(1.0),
	m_lmin(2e-3),
	m_nmin(1e-4),
	m_pw(-1.0),
	m_tq(-0.04),
	m_qq(5e-3),
	m_q10(5.0),
	m_bLogNGate(false),
	m_bLogLGate(false)
{}
UG_CATCH_THROW("Error in KaDist_Golding01 initializer list.");


template<typename TDomain>
KaDist_Golding01<TDomain>::
KaDist_Golding01
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try :
	ICableMembraneTransport<TDomain>(functions, subsets),
	m_ek(-0.09),
	m_gkabar(80.0),
	m_vhalfn(-1e-3),
	m_vhalfl(-0.056),
	m_a0n(100.0),
	m_zetan(-1.8),
	m_zetal(3.0),
	m_gmn(0.39),
	m_gml(1.0),
	m_lmin(2e-3),
	m_nmin(1e-4),
	m_pw(-1.0),
	m_tq(-0.04),
	m_qq(5e-3),
	m_q10(5.0),
	m_bLogNGate(false),
	m_bLogLGate(false)
{}
UG_CATCH_THROW("Error in KaDist_Golding01 initializer list.");


template<typename TDomain>
KaDist_Golding01<TDomain>::
~KaDist_Golding01()
{}



template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_ek(number val)
{ 
	m_ek = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_gkabar(number val)
{ 
	m_gkabar = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_vhalfn(number val)
{ 
	m_vhalfn = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_vhalfl(number val)
{ 
	m_vhalfl = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_a0n(number val)
{ 
	m_a0n = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_zetan(number val)
{ 
	m_zetan = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_zetal(number val)
{ 
	m_zetal = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_gmn(number val)
{ 
	m_gmn = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_gml(number val)
{ 
	m_gml = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_lmin(number val)
{ 
	m_lmin = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_nmin(number val)
{ 
	m_nmin = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_pw(number val)
{ 
	m_pw = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_tq(number val)
{ 
	m_tq = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_qq(number val)
{ 
	m_qq = val;
}

template<typename TDomain> 
void KaDist_Golding01<TDomain>::set_q10(number val)
{ 
	m_q10 = val;
}

template<typename TDomain>
void KaDist_Golding01<TDomain>::
set_logLGate(bool bLoglGate)
{
	m_bLogLGate = bLoglGate;
}

template<typename TDomain>
void KaDist_Golding01<TDomain>::
set_logNGate(bool bLognGate)
{
	m_bLogNGate = bLognGate;
}



template<typename TDomain> 
number KaDist_Golding01<TDomain>::ek() const
{ 
	return m_ek;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::gkabar() const
{ 
	return m_gkabar;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::vhalfn() const
{ 
	return m_vhalfn;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::vhalfl() const
{ 
	return m_vhalfl;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::a0n() const
{ 
	return m_a0n;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::zetan() const
{ 
	return m_zetan;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::zetal() const
{ 
	return m_zetal;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::gmn() const
{ 
	return m_gmn;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::gml() const
{ 
	return m_gml;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::lmin() const
{ 
	return m_lmin;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::nmin() const
{ 
	return m_nmin;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::pw() const
{ 
	return m_pw;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::tq() const
{ 
	return m_tq;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::qq() const
{ 
	return m_qq;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::q10() const
{ 
	return m_q10;
}





template<typename TDomain>
void KaDist_Golding01<TDomain>::
init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	const number v = vrt_values[CableEquation<TDomain>::_v_];

	number a = alpn(v);
	m_aaNGate[vrt] = 1.0 / (1.0 + a);

	a = alpl(v);
	m_aaLGate[vrt] = 1.0 / (1.0 + a);
}


template<typename TDomain> 
void KaDist_Golding01<TDomain>::
update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{ 
	const number celsius = m_pCE->temperature_celsius();

	const number dt = newTime - m_pCE->time();
	const number v = vrt_values[CableEquation<TDomain>::_v_];
	number& n = m_aaNGate[vrt];
	number& l = m_aaLGate[vrt];

	const number qt = pow(m_q10, (celsius - 24.0) / 10.0);
	number a = alpn(v);
	const number ninf = 1.0 / (1.0 + a);
	const number taun = std::max(m_nmin, betn(v) / (qt * m_a0n * (1.0 + a)));

	// we give the exact solution (instead of implicit Euler)
	n = ninf + (n - ninf) * exp(-dt/taun);


	a = alpl(v);
	const number linf = 1.0 / (1.0 + a);
	const number taul = std::max(m_lmin, 0.26 * (v + 0.05));

	// we give the exact solution (instead of implicit Euler)
	l = linf + (l - linf) * exp(-dt/taul);
}


template<typename TDomain>
void KaDist_Golding01<TDomain>::current
(
	Vertex* vrt,
	const std::vector<number>& vrt_values,
	std::vector<number>& outCurrentValues
)
{
	const number v = vrt_values[CableEquation<TDomain>::_v_];
	outCurrentValues.push_back(m_gkabar * m_aaNGate[vrt] * m_aaLGate[vrt] * (v - m_ek));
}


template<typename TDomain>
void KaDist_Golding01<TDomain>::ce_obj_available()
{
	init_attachments();
}


template<typename TDomain> 
std::vector<number> KaDist_Golding01<TDomain>::state_values(number x, number y, number z) const
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
	try { ssGrp = SubsetGroup(m_pCE->approx_space()->domain()->subset_handler(), this->m_vSubset);}
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


template<typename TDomain>
void KaDist_Golding01<TDomain>::specify_write_function_indices()
{
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_);
}




template<typename TDomain> 
number KaDist_Golding01<TDomain>::alpn(number v)
{
	const number F = m_pCE->F;
	const number R = m_pCE->R;
	const number T = m_pCE->temperature();
	const number zeta = m_zetan + m_pw / (1.0 + exp((v-m_tq) / m_qq));
	return exp(zeta * (v - m_vhalfn) * F/(R*T));
}

template<typename TDomain>
number KaDist_Golding01<TDomain>::betn(number v)
{
	const number F = m_pCE->F;
	const number R = m_pCE->R;
	const number T = m_pCE->temperature();
	const number zeta = m_zetan + m_pw / (1.0 + exp((v-m_tq) / m_qq));
	return exp(zeta * m_gmn * (v - m_vhalfn) * F/(R*T)) ;
}

template<typename TDomain> 
number KaDist_Golding01<TDomain>::alpl(number v)
{
	const number F = m_pCE->F;
	const number R = m_pCE->R;
	const number T = m_pCE->temperature();
	return exp(m_zetal * (v - m_vhalfl) * F/(R*T)) ;
}

template<typename TDomain>
number KaDist_Golding01<TDomain>::betl(number v)
{
	const number F = m_pCE->F;
	const number R = m_pCE->R;
	const number T = m_pCE->temperature();
	return exp(m_zetal * m_gml * (v - m_vhalfl) * F/(R*T)) ;
}


template<typename TDomain>
void KaDist_Golding01<TDomain>::init_attachments()
{
	SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();

	if (spGrid->has_vertex_attachment(m_aNGate))
		UG_THROW("Attachment necessary (nGate) for KaDist_Golding01 channel dynamics "
			"could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_aNGate);
	m_aaNGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aNGate);

	if (spGrid->has_vertex_attachment(this->m_aLGate))
		UG_THROW("Attachment necessary (lGate) for KaDist_Golding01 channel dynamics "
			"could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_aLGate);
	m_aaLGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aLGate);
}



//	explicit template instantiations
#ifdef UG_DIM_1
template class KaDist_Golding01<Domain1d>;
#endif
#ifdef UG_DIM_2
template class KaDist_Golding01<Domain2d>;
#endif
#ifdef UG_DIM_3
template class KaDist_Golding01<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
