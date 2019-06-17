/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Pascal Gottmann
 * Creation date: 2014-06-13
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

#include "implicit_active_cable_disc_nernst.h"

#include "common/error.h"  // for UG_THROW
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"  // for GeometryProvider
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"  // for FV1Geometry

#include <cmath>


namespace ug {
namespace cable_neuron {


template<typename TDomain>
ImplicitActiveCableDiscNernst<TDomain>::
ImplicitActiveCableDiscNernst(const char* functions, const char* subsets)
: ImplicitActiveCableDiscBase<TDomain>(functions,subsets),
  m_diff_K(1.96e-9), m_diff_Na(2.03e-9),
  m_kOut(4.0), m_naOut(140.0)
{}


template<typename TDomain>
void ImplicitActiveCableDiscNernst<TDomain>::
set_diffusion_constants(number diffK, number diffNa)
{
	m_diff_K = diffK;
	m_diff_Na = diffNa;
}


template<typename TDomain>
void ImplicitActiveCableDiscNernst<TDomain>::
set_outside_concs(number concK, number concNa)
{
	m_kOut = concK;
	m_naOut = concNa;
}


template<typename TDomain>
void ImplicitActiveCableDiscNernst<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check number of unknown functions
	if (vLfeID.size() != 6)
		UG_THROW("ImplicitActiveCableDiscNernst: Wrong number of functions given. Need exactly "<< 6);

	// check shape function type
	if (vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("ImplicitActiveCableDiscNernst scheme only implemented for 1st order.");

	// update assemble functions
	register_all_funcs();
}


// ///////////////////////
// Assembling functions //
// ///////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	// set local positions
	static const int refDim = TElem::dim;
	TFVGeom& geo = GeomProvider<TFVGeom>::get();

	const MathVector<refDim>* vSCVip = geo.scv_local_ips();
	const size_t numSCVip = geo.num_scv_ips();
	const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
	const size_t numSCVFip = geo.num_scvf_ips();

	m_spec_cap.template set_local_ips<refDim>(vSCVip,numSCVip, false);
	m_spec_res.template set_local_ips<refDim>(vSCVFip,numSCVFip, false);
	m_gK.template set_local_ips<refDim>(vSCVip,numSCVip, false);
	m_gNa.template set_local_ips<refDim>(vSCVip,numSCVip, false);
	m_gL.template set_local_ips<refDim>(vSCVip,numSCVip, false);
	m_eK.template set_local_ips<refDim>(vSCVip,numSCVip, false);
	m_eNa.template set_local_ips<refDim>(vSCVip,numSCVip, false);
	m_eL.template set_local_ips<refDim>(vSCVip,numSCVip, false);
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("ElemDiscHHNernstFV1::prep_elem: Cannot update finite volume geometry.");

	// set global positions
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();

	m_spec_cap.set_global_ips(vSCVip, numSCVip);
	m_spec_res.set_global_ips(vSCVFip, numSCVFip);
	m_gK.set_global_ips(vSCVip, numSCVip);
	m_gNa.set_global_ips(vSCVip, numSCVip);
	m_gL.set_global_ips(vSCVip, numSCVip);
	m_eK.set_global_ips(vSCVip, numSCVip);
	m_eNa.set_global_ips(vSCVip, numSCVip);
	m_eL.set_global_ips(vSCVip, numSCVip);
}



static number vtrap(number x, number y)
{
    if (fabs(x/y) < 1e-6)
    	return y*(1 - x/y/2);
    return x/(exp(x/y) - 1);
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	number element_length = 0.0;
	number pre_resistance = 0.0;
	number volume = 0.0;

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

// channel kinetics
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// get membrane potential and gating params
		const number& vm = u(_vm_,co);
		const number& n = u(_n_,co);
		const number& m = u(_m_,co);
		const number& h = u(_h_,co);

		// get diam from attachment for element
		number diam = m_aaDiameter[pElem->vertex(co)];

		volume += scv.volume();
		// add length of scv to element length
		element_length += scv.volume();
		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*diam*diam);


		// values for m gate
		number AlphaHm = 1e5 * vtrap(-(vm+0.04), 0.010);
		number BetaHm =  4e3 * exp(-(vm+0.065)/0.018);

		// values for h gate
		number AlphaHh = 70.0 * exp(-(vm+0.065)/0.020);
		number BetaHh = 1e3 / (exp(-(vm+0.035)/0.010) + 1.0);

		// values for n gate
		number AlphaHn = 1e4*vtrap(-(vm+0.055), 0.010);
		number BetaHn = 125.0*exp(-(vm+0.065)/0.080);

		const number tmp = m_T - 273.15;
		const number tmp_factor = m_bTempDep ? std::pow(2.3, (tmp-23.0)/10.0) : 1.0;
		number rate_h = -tmp_factor*((AlphaHh * (1.0-h)) - BetaHh * h);
		number rate_m = -tmp_factor*((AlphaHm * (1.0-m)) - BetaHm * m);
		number rate_n = -tmp_factor*((AlphaHn * (1.0-n)) - BetaHn * n);


		const number helpV = (m_R*m_T)/m_F;
		const number potassium_nernst_eq = helpV*(log(m_kOut/u(_K_,co)));
		const number sodium_nernst_eq = helpV*(log(m_naOut/u(_Na_,co)));

		// single channel type fluxes
		const number potassium_part_of_flux = m_gK[ip] * pow(n,4) * (vm - potassium_nernst_eq);
		const number sodium_part_of_flux =  m_gNa[ip] * pow(m,3) * h * (vm - sodium_nernst_eq);
		const number leakage_part_of_flux = m_gL[ip] * (vm - m_eL[ip]);


		// injection flux
		number inject = 0.0;
		const number time = this->time();
		const MathVector<dim>& coco = vCornerCoords[ip];
		if (dim == 3)
			(*m_Injection)(inject, 4, time, coco[0], coco[1], coco[2]);

		if (dim == 2)
			(*m_Injection)(inject, 3, time, coco[0], coco[1]);

		if (dim == 1)
			(*m_Injection)(inject, 2, time, coco[0]);


		const number flux = potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux;

		d(_vm_, co) += scv.volume()*PI*diam*(flux - inject);
		d(_h_, co) += rate_h;
		d(_m_, co) += rate_m;
		d(_n_, co) += rate_n;

		d(_K_, co)  += potassium_part_of_flux/m_F * PI*diam*scv.volume();
		d(_Na_, co) += sodium_part_of_flux/m_F * PI*diam*scv.volume();
	}

// cable equation, axial current part; diffusion
	MathVector<dim> grad_c, grad_k, grad_na;

	for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		// compute gradient at ip
		VecSet(grad_c, 0.0);
		VecSet(grad_k, 0.0);
		VecSet(grad_na, 0.0);

		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
			VecScaleAppend(grad_c, u(_vm_,sh), scvf.global_grad(sh));
			VecScaleAppend(grad_k, u(_K_,sh), scvf.global_grad(sh));
			VecScaleAppend(grad_na, u(_Na_,sh), scvf.global_grad(sh));
			}

		// scalar product with normal
		number diff_flux = VecDot(grad_c, scvf.normal());
		number diff_fluxNa = VecDot(grad_na, scvf.normal());
		number diff_fluxK = VecDot(grad_k, scvf.normal());

		number diam_fromTo = std::min(m_aaDiameter[pElem->vertex(scvf.from())],
		                          	  m_aaDiameter[pElem->vertex(scvf.to())]);

		// scale by 1/resistance and by length of element
		diff_flux *= element_length / (m_spec_res[ip]*pre_resistance);

		// diffusion for Na and K
		diff_fluxNa *= m_diff_Na * 0.25*PI*diam_fromTo*diam_fromTo;
		diff_fluxK *= m_diff_K * 0.25*PI*diam_fromTo*diam_fromTo;

		// add to local defect
		d(_vm_, scvf.from()) -= diff_flux;
		d(_vm_, scvf.to()) += diff_flux;
		d(_K_, scvf.from()) -= diff_fluxK;
		d(_K_, scvf.to()) += diff_fluxK;
		d(_Na_, scvf.from()) -= diff_fluxNa;
		d(_Na_, scvf.to()) += diff_fluxNa;
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		//get Diameter from element
		number diam = m_aaDiameter[pElem->vertex(co)];

		// get spec capacity
		number spec_capacity = m_spec_cap[ip];

		// gating parameters time derivative
		d(_h_, co) += u(_h_, co);
		d(_m_, co) += u(_m_, co);
		d(_n_, co) += u(_n_, co);

		// Nernst Paras time derivative
		d(_Na_, co) += u(_Na_, co)*scv.volume()*0.25*PI*diam*diam;
		d(_K_, co)  += u(_K_, co)*scv.volume()*0.25*PI*diam*diam;

		// potential equation time derivative
		d(_vm_, co) += PI*diam*scv.volume()*u(_vm_, co)*spec_capacity;
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	number element_length = 0.0;
	number pre_resistance = 0.0;

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

// channel kinetics derivatives
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// get membrane potential and gating params
		const number& vm = u(_vm_,co);
		const number& n = u(_n_,co);
		const number& m = u(_m_,co);
		const number& h = u(_h_,co);

		// get diameter from element
		number diam = m_aaDiameter[pElem->vertex(co)];

		// add length of scv to element length
		element_length += scv.volume();
		
		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*diam*diam);


		// values for m gate
		number AlphaHm = 1e5 * vtrap(-(vm+0.040), 0.010);
		number BetaHm =  4e3 * exp(-(vm+0.065)/0.018);

		// values for h gate
		number AlphaHh = 70.0 * exp(-(vm+0.065)/0.020);
		number BetaHh = 1e3 / (exp(-(vm+0.035)/0.010) + 1.0);

		// values for n gate
		number AlphaHn = 1e4*vtrap(-(vm+0.055), 0.010);
		number BetaHn = 125.0*exp(-(vm+0.065)/0.080);


		// gating param m derivatives
		number help = -100.0*(vm+0.040);
		number dAlphaHm_dVm;
		if (fabs(exp(help)-1.0) > 1e-6)
			dAlphaHm_dVm = -1e5 * ((1.0-help)*exp(help)-1.0) / pow(exp(help)-1.0, 2);
		else
			dAlphaHm_dVm = 5e4;

		number dBetaHm_dVm = -4e3/0.018 * exp(-(vm+0.065)/0.018);

		// gating param h derivatives
		number dAlphaHh_dVm = -70.0/0.020 * exp(-(vm+0.065)/0.020);
		help = exp(-100.0*(vm+0.035));
		number dBetaHh_dVm = 1e5*help / pow(help+1.0, 2);

		// gating param n derivatives
		help = 100.0*(vm+0.055);
		number dAlphaHn_dVm;
		if (fabs(exp(help)-1.0) > 1e-6)
			dAlphaHn_dVm = -1e4 * ((1.0-help)*exp(help)-1.0) / pow(exp(help)-1.0, 2);
		else
			dAlphaHn_dVm = 5e4;

		number dBetaHn_dVm = 125.0/0.08 * exp((vm+0.065)/0.080);


		// Nernst potential of potassium and sodium
		const number helpV = (m_R*m_T)/m_F;	// unit must be mV!
		const number potassium_nernst_eq = helpV*(log(m_kOut/u(_K_,co)));
		const number sodium_nernst_eq = helpV*(log(m_naOut/u(_Na_,co)));

		// derivatives of nernst equations
		const number potassium_nernst_eq_dK = -helpV / u(_K_,co);
		const number sodium_nernst_eq_dNa = -helpV / u(_Na_,co);

		number factor = PI*diam*scv.volume() / m_F * m_gK[ip];
		J(_K_, co, _K_, co)   +=  -factor * pow(n,4) * potassium_nernst_eq_dK;
		J(_K_, co, _n_, co)   +=  factor * 4*pow(n,3) * (vm - potassium_nernst_eq);
		J(_K_, co, _vm_, co)  +=  factor * pow(n,4);

		factor = PI*diam*scv.volume() / m_F * m_gNa[ip];
		J(_Na_, co, _Na_, co)	+= -factor * pow(m,3) * h * sodium_nernst_eq_dNa;
		J(_Na_, co, _m_, co)	+= factor * 3*pow(m,2) * h * (vm - sodium_nernst_eq);
		J(_Na_, co, _h_, co)	+= factor * pow(m,3) * (vm - sodium_nernst_eq);
		J(_Na_, co, _vm_, co)	+= factor * pow(m,3) * h;

		J(_vm_, co, _K_,co) += -scv.volume()*PI*diam * m_gK[ip] * pow(n,4) * potassium_nernst_eq_dK;
		J(_vm_, co, _Na_,co) += -scv.volume()*PI*diam * m_gNa[ip] * pow(m,3) * h * sodium_nernst_eq_dNa;

		// derivatives of channel states
		const number tmp = m_T - 273.15;
		const number tmp_factor = m_bTempDep ? std::pow(2.3, (tmp-23.0)/10.0) : 1.0;
		J(_h_, co, _h_, co) += tmp_factor * (AlphaHh + BetaHh);
		J(_m_, co, _m_, co) += tmp_factor * (AlphaHm + BetaHm);
		J(_n_, co, _n_, co) += tmp_factor * (AlphaHn + BetaHn);
		J(_h_, co, _vm_, co) += -tmp_factor * ((dAlphaHh_dVm * (1.0-h)) - dBetaHh_dVm * h);
		J(_m_, co, _vm_, co) += -tmp_factor * ((dAlphaHm_dVm * (1.0-m)) - dBetaHm_dVm * m);
		J(_n_, co, _vm_, co) += -tmp_factor * ((dAlphaHn_dVm * (1.0-n)) - dBetaHn_dVm * n);

		// derivatives of potential from HH channels
		J(_vm_, co, _h_, co) += scv.volume()*PI*diam * m_gNa[ip]*pow(m,3) * (vm - sodium_nernst_eq);
		J(_vm_, co, _m_, co) += scv.volume()*PI*diam * 3.0*m_gNa[ip]*pow(m,2) * h * (vm - sodium_nernst_eq);
		J(_vm_, co, _n_, co) += scv.volume()*PI*diam * 4.0*m_gK[ip]*pow(n,3) * (vm - potassium_nernst_eq);
		J(_vm_, co, _vm_, co) += scv.volume()*PI*diam * (m_gK[ip]*pow(n,4) + m_gNa[ip]*pow(m,3)*h + m_gL[ip]);
	}

// axial current and diffusion derivatives

	// loop sub-control volume faces (SCVF)
	for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		// loop shape functions
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			// scalar product with normal
			number diam_fromTo = std::min(m_aaDiameter[pElem->vertex(scvf.from())],
										  m_aaDiameter[pElem->vertex(scvf.to())]);

			number d_diff_flux = VecDot(scvf.global_grad(sh), scvf.normal());
			number d_diff_fluxNa = VecDot(scvf.global_grad(sh), scvf.normal());
			number d_diff_fluxK = VecDot(scvf.global_grad(sh), scvf.normal());

			// scale by 1/resistance and by length of element
			d_diff_flux *= element_length / (m_spec_res[ip]*pre_resistance);
			d_diff_fluxNa *= m_diff_Na * 0.25*PI*diam_fromTo*diam_fromTo;
			d_diff_fluxK *= m_diff_K * 0.25*PI*diam_fromTo*diam_fromTo;

			// add flux term to local matrix
			J(_vm_, scvf.from(), _vm_, sh) -= d_diff_flux;
			J(_vm_, scvf.to(), _vm_, sh) += d_diff_flux;

			J(_K_, scvf.from(), _K_, sh) -= d_diff_fluxK;
			J(_K_, scvf.to(), _K_, sh) += d_diff_fluxK;

			J(_Na_, scvf.from(), _Na_, sh) -= d_diff_fluxNa;
			J(_Na_, scvf.to(), _Na_, sh) += d_diff_fluxNa;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		//get Diameter from element
		number diam = m_aaDiameter[pElem->vertex(co)];

		// get spec capacity
		number spec_capacity = m_spec_cap[ip];

		// gating parameters
		J(_h_, co, _h_, co) += 1.0;
		J(_m_, co, _m_, co) += 1.0;
		J(_n_, co, _n_, co) += 1.0;

		J(_Na_, co, _Na_, co) += 1.0*scv.volume()*0.25*PI*diam*diam;
		J(_K_, co, _K_, co) += 1.0*scv.volume()*0.25*PI*diam*diam;

		// potential equation
		J(_vm_, co, _vm_, co) += PI*diam*scv.volume()*spec_capacity;
	}
}


// register assemble functions
template<typename TDomain>
void ImplicitActiveCableDiscNernst<TDomain>::
register_all_funcs()
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ImplicitActiveCableDiscNernst<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef ImplicitActiveCableDiscNernst<TDomain> T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TFVGeom>);
}



// explicit template instantiations
#ifdef UG_DIM_3
template class ImplicitActiveCableDiscNernst<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug

