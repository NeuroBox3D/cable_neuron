/*
 * nax_g01.cpp
 *
 *  Created on: 2019-06-12
 *      Author: mbreit
 */

#include <cmath>
#include <limits>

#include "nax_golding01.h"

namespace ug { 
namespace cable_neuron { 


template<typename TDomain>
Nax_Golding01<TDomain>::
Nax_Golding01(const char* functions, const char* subsets)
try :
	ICableMembraneTransport<TDomain>(functions, subsets),
#ifdef UG_FOR_LUA
	m_bConductanceDependsOnCoordinates(false),
	m_spCondFct(SPNULL),
#endif
	m_gbar(100.0),
	m_tha(-0.03),
	m_qa(0.0072),
	m_Ra(400.0),
	m_Rb(124.0),
	m_thi1(-0.045),
	m_thi2(-0.045),
	m_qd(0.0015),
	m_qg(0.0015),
	m_mmin(2e-5),
	m_hmin(5e-4),
	m_q10(2.0),
	m_Rg(10.0),
	m_Rd(30.0),
	m_thinf(-0.05),
	m_qinf(0.004),
	m_bLogMGate(false),
	m_bLogHGate(false)
{}
UG_CATCH_THROW("Error in Nax_Golding01 initializer list.");


template<typename TDomain>
Nax_Golding01<TDomain>::
Nax_Golding01
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
	m_gbar(100.0),
	m_tha(-0.03),
	m_qa(0.0072),
	m_Ra(400.0),
	m_Rb(124.0),
	m_thi1(-0.045),
	m_thi2(-0.045),
	m_qd(0.0015),
	m_qg(0.0015),
	m_mmin(2e-5),
	m_hmin(5e-4),
	m_q10(2.0),
	m_Rg(10.0),
	m_Rd(30.0),
	m_thinf(-0.05),
	m_qinf(0.004),
	m_bLogMGate(false),
	m_bLogHGate(false)
{}
UG_CATCH_THROW("Error in Nax_Golding01 initializer list.");


template<typename TDomain>
Nax_Golding01<TDomain>::
~Nax_Golding01()
{
#ifdef UG_FOR_LUA
	if (m_bConductanceDependsOnCoordinates)
	{
		SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();
		if (spGrid->has_vertex_attachment(m_aGBar))
			spGrid->detach_from_vertices(m_aGBar);
	}
#endif
}


template <typename TDomain>
std::string Nax_Golding01<TDomain>::
name()
{
	return std::string("Nax_Golding01");
}



template<typename TDomain> 
void Nax_Golding01<TDomain>::set_gbar(number val)
{ 
	m_gbar = val;
}

template<typename TDomain> 
void Nax_Golding01<TDomain>::set_tha(number val)
{ 
	m_tha = val;
}

template<typename TDomain> 
void Nax_Golding01<TDomain>::set_qa(number val)
{ 
	m_qa = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_Ra(number val)
{ 
	m_Ra = val;
}

template<typename TDomain> 
void Nax_Golding01<TDomain>::set_Rb(number val)
{ 
	m_Rb = val;
}

template<typename TDomain> 
void Nax_Golding01<TDomain>::set_thi1(number val)
{ 
	m_thi1 = val;
}

template<typename TDomain> 
void Nax_Golding01<TDomain>::set_thi2(number val)
{ 
	m_thi2 = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_qd(number val)
{
	m_qd = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_qg(number val)
{
	m_qg = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_mmin(number val)
{
	m_mmin = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_hmin(number val)
{
	m_hmin = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_q10(number val)
{
	m_q10 = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_Rg(number val)
{
	m_Rg = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_Rd(number val)
{
	m_Rd = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_thinf(number val)
{
	m_thinf = val;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_qinf(number val)
{
	m_qinf = val;
}


#ifdef UG_FOR_LUA
template <typename TDomain>
void Nax_Golding01<TDomain>::set_gbar(SmartPtr<LuaUserData<number, TDomain::dim> > fct)
{
	m_spCondFct = fct;
	m_bConductanceDependsOnCoordinates = true;
}

template <typename TDomain>
void Nax_Golding01<TDomain>::set_gbar(const char* fct)
{
	m_spCondFct = LuaUserDataFactory<number, TDomain::dim>::create(fct);
	m_bConductanceDependsOnCoordinates = true;
}
#endif


template<typename TDomain>
void Nax_Golding01<TDomain>::set_logMGate(bool bLogmGate)
{
	m_bLogMGate = bLogmGate;
}

template<typename TDomain>
void Nax_Golding01<TDomain>::set_logHGate(bool bLoghGate)
{
	m_bLogHGate = bLoghGate;
}





static number trap0(number v, number th, number a, number q)
{
	if (fabs(v-th) > 1e-6)
		return a * (v - th) / (1.0 - exp(-(v - th) / q));
	return a*q;
}


template<typename TDomain>
void Nax_Golding01<TDomain>::
init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	// init gating values
	const number v = vrt_values[CableEquation<TDomain>::_v_];

	const number a = trap0(v, m_tha, m_Ra, m_qa);
	const number b = trap0(-v, -m_tha, m_Rb, m_qa);
	const number minf = a / (a+b);
	m_aaMGate[vrt] = minf;

	const number hinf = 1.0 / (1.0 + exp((v-m_thinf) / m_qinf));
	m_aaHGate[vrt] = hinf;


#ifdef UG_FOR_LUA
	// write conductance value to attachment
	if (m_bConductanceDependsOnCoordinates)
	{
		const int si = m_pCE->subset_handler().get_subset_index(vrt);
		const MathVector<TDomain::dim>& coords =
			this->m_pCE->approx_space()->domain()->position_accessor()[vrt];
		m_spCondFct->evaluate(m_aaGBar[vrt], coords, 0.0, si);
	}
#endif
}


template<typename TDomain> 
void Nax_Golding01<TDomain>::
update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{ 
	const number celsius = m_pCE->temperature_celsius();

	const number v = vrt_values[CableEquation<TDomain>::_v_];
	number& m = m_aaMGate[vrt];
	number& h = m_aaHGate[vrt];

	const number qt = pow(m_q10, (celsius - 24.0) / 10.0);
	number a = trap0(v, m_tha, m_Ra, m_qa);
	number b = trap0(-v, -m_tha, m_Rb, m_qa);
	const number taum = std::max(m_mmin, 1.0 / (qt * (a+b)));
	const number minf = a / (a+b);

	a = trap0(v, m_thi1, m_Rd, m_qd);
	b = trap0(-v, -m_thi2, m_Rg, m_qg);
	const number tauh = std::max(m_hmin, 1.0 / (qt * (a+b)));
	const number hinf = 1.0 / (1.0 + exp((v-m_thinf) / m_qinf));

	// we give the exact solution (instead of implicit Euler)
	const number dt = newTime - m_pCE->time();
	m = minf + (m - minf) * exp(-dt/taum);
	h = hinf + (h - hinf) * exp(-dt/tauh);
}


template<typename TDomain>
void Nax_Golding01<TDomain>::current
(
	Vertex* vrt,
	const std::vector<number>& vrt_values,
	std::vector<number>& outCurrentValues
)
{
	const number v = vrt_values[CableEquation<TDomain>::_v_];
	const number ena = this->m_pCE->rev_pot_na();

	const number m = m_aaMGate[vrt];
	const number h = m_aaHGate[vrt];

	number gbar = m_gbar;
#ifdef UG_FOR_LUA
	if (m_bConductanceDependsOnCoordinates)
		gbar = m_aaGBar[vrt];
#endif

	outCurrentValues.push_back(gbar * m*m*m*h * (v - ena));
}


template<typename TDomain>
void Nax_Golding01<TDomain>::ce_obj_available()
{
	init_attachments();
}


template<typename TDomain> 
std::vector<number> Nax_Golding01<TDomain>::state_values(number x, number y, number z) const
{ 
#ifdef UG_PARALLLEL
	if (pcl::NumProcs() > 1)
		UG_THROW("The state_values method is not parallelized.");
#endif

	std::vector<number> states;

	if (!m_bLogMGate && !m_bLogHGate)
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

	if (m_bLogMGate)
		states.push_back(m_aaMGate[bestVrt]);
	if (m_bLogHGate)
		states.push_back(m_aaHGate[bestVrt]);

	return states;
}


template<typename TDomain>
void Nax_Golding01<TDomain>::specify_write_function_indices()
{
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_);
}


template<typename TDomain>
void Nax_Golding01<TDomain>::init_attachments()
{
	SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();

	if (spGrid->has_vertex_attachment(m_aMGate))
		UG_THROW("Attachment necessary (mGate) for Nax_Golding01 channel dynamics "
			"could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_aMGate);
	m_aaMGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aMGate);

	if (spGrid->has_vertex_attachment(m_aHGate))
		UG_THROW("Attachment necessary (hGate) for Nax_Golding01 channel dynamics "
			"could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_aHGate);
	m_aaHGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aHGate);

#ifdef UG_FOR_LUA
	// attach conductance attachment if necessary
	if (m_bConductanceDependsOnCoordinates)
	{
		if (spGrid->has_vertex_attachment(m_aGBar))
			UG_THROW("Attachment necessary (gBar) for " << name() << " channel dynamics "
				"could not be made, since it already exists.");
		spGrid->attach_to_vertices(m_aGBar);
		m_aaGBar = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_aGBar);
	}
#endif
}



//	explicit template instantiations
#ifdef UG_DIM_3
template class Nax_Golding01<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
