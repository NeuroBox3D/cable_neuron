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

#ifndef UG__PLUGINS__CABLE_NEURON__CABLE_DISC__IMPLICIT_ACTIVE_CABLE_DISC_NERNST_H
#define UG__PLUGINS__CABLE_NEURON__CABLE_DISC__IMPLICIT_ACTIVE_CABLE_DISC_NERNST_H

#include "../cable_disc/implicit_active_cable_disc_base.h"

namespace ug {
namespace cable_neuron {


template <typename TDomain>
class ImplicitActiveCableDiscNernst
: public ImplicitActiveCableDiscBase<TDomain>
{
	private:
		///	base class type
		typedef ImplicitActiveCableDiscBase<TDomain> base_type;

	public:
		///	world dimension
		static const int dim = base_type::dim;

	protected:
		using base_type::_vm_;
		using base_type::_n_;
		using base_type::_m_;
		using base_type::_h_;

		static const size_t _K_ = 4;
		static const size_t _Na_ = 5;

	public:
		///	constructor
		ImplicitActiveCableDiscNernst(const char* functions, const char* subsets);

		/// sets diffusion consts
		void set_diffusion_constants(number diffK, number diffNa);

		/// set outside concentrations
		void set_outside_concs(number concK, number concNa);

	private:
		///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		///	prepares the loop over all elements
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	prepares the element for assembling
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[]);

		///	finishes the loop over all elements
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

		///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the stiffness part of the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the mass part of the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the local right hand side
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	protected:
		/// registration util
		/// \{
		void register_all_funcs();

		template <typename TElem, typename TFVGeom>
		void register_func();
		/// \}

	protected:
		/// diffusion constants
		/// \{
		number m_diff_K;
		number m_diff_Na;
		/// \}

		/// outside concentrations
		/// \{
		number m_kOut;
		number m_naOut;
		/// \}
		///
		using base_type::m_R;
		using base_type::m_T;
		using base_type::m_F;

		using base_type::m_spec_cap;
		using base_type::m_spec_res;

		using base_type::m_gK;
		using base_type::m_gNa;
		using base_type::m_gL;

		using base_type::m_eK;
		using base_type::m_eNa;
		using base_type::m_eL;

		using base_type::m_Injection;

		using base_type::m_aaDiameter;

		using base_type::m_bTempDep;
};


} // end namespace cable_neuron
} // end namespace ug


#endif // UG__PLUGINS__CABLE_NEURON__CABLE_DISC__IMPLICIT_ACTIVE_CABLE_DISC_NERNST_H
