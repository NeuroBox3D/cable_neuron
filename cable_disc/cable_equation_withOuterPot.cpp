/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-01-09
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

#include "../cable_disc/cable_equation_withOuterPot.h"


namespace ug {
namespace cable_neuron {


template<typename TDomain>
CableEquationWithOuterPot<TDomain>::CableEquationWithOuterPot(const char* fcts, const char* subsets)
: CableEquation<TDomain>(fcts, subsets)
{}


template<typename TDomain>
void CableEquationWithOuterPot<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if (vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("CableEquation FV scheme only implemented for 1st order.");

	// remember
	this->m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	my_register_all_funcs(this->m_bNonRegularGrid);
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CableEquationWithOuterPot<TDomain>::
my_add_def_A_elem
(
	LocalVector& d,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type (in order to allow access to attachments)
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	// calculate some helper variables for cable equation
	number element_length = 0.0;
	number pre_resistance = 0.0;
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// get diam from attachment
		number diam = this->m_aaDiameter[pElem->vertex(co)];

		// add length of scv to element length
		element_length += scv.volume();

		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*diam*diam);
	}

	// diffusive parts
	MathVector<dim> grad_c;
	for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		// compute inner(!) potential gradient at ip
		VecSet(grad_c, 0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad_c, u(_VM_, sh) + u(_PHIO_, sh), scvf.global_grad(sh));

		// scalar product with normal
		number grad_normal = VecDot(grad_c, scvf.normal());

		// scale by 1/resistance and by length of element
		number diff_flux = grad_normal * element_length / (this->m_spec_res*pre_resistance);

		// add to local defect of VM
		d(_VM_, scvf.from()) -= diff_flux;
		d(_VM_, scvf.to()  ) += diff_flux;
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CableEquationWithOuterPot<TDomain>::
my_add_jac_A_elem
(
	LocalMatrix& J,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// some helper vars
	number element_length = 0.0;
	number pre_resistance = 0.0;

	// cast elem to appropriate type  (in order to allow access to attachments)
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	// only helper calculations for axial current here (membrane fluxes are purely explicit)
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// get diam from attachment
		number diam = this->m_aaDiameter[pElem->vertex(co)];

		// add length of scv to element length
		element_length += scv.volume();

		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*diam*diam);
	}


	// diffusive part
	for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		// loop shape functions
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			// scalar product with normal
			number grad_normal = VecDot(scvf.global_grad(sh), scvf.normal());

			// scale by 1/resistance and by length of element
			number d_diff_flux = grad_normal * element_length / (this->m_spec_res*pre_resistance);

			// add flux term to local matrix
			J(_VM_, scvf.from(), _VM_, sh) -= d_diff_flux;
			J(_VM_, scvf.to()  , _VM_, sh) += d_diff_flux;
			J(_VM_, scvf.from(), _PHIO_, sh) += d_diff_flux;
			J(_VM_, scvf.to()  , _PHIO_, sh) -= d_diff_flux;
		}
	}
}


// ///////////////////////////////
//	register assemble functions //
// ///////////////////////////////

template<typename TDomain>
void CableEquationWithOuterPot<TDomain>::
my_register_all_funcs(bool bHang)
{
	// call parent registration routines
	this->register_all_funcs(bHang);

	// register own assembling functionality
	my_register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CableEquationWithOuterPot<TDomain>::
my_register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef CableEquationWithOuterPot<TDomain> T;

	this->set_add_def_A_elem_fct(id, &T::template my_add_def_A_elem<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template my_add_jac_A_elem<TElem, TFVGeom>);
}



// ////////////////////////////////////
//	explicit template instantiations //
// ////////////////////////////////////
#ifdef UG_DIM_3
	template class CableEquationWithOuterPot<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
