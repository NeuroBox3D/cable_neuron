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

#ifndef __UG__PLUGINS__CABLE_NEURON__CABLE_DISC__CABLE_EQUATION_WITHOUTERPOT_H__
#define __UG__PLUGINS__CABLE_NEURON__CABLE_DISC__CABLE_EQUATION_WITHOUTERPOT_H__

// other ug4 modules
#include "cable_equation.h"

namespace ug {
namespace cable_neuron {


template <typename TDomain>
class CableEquationWithOuterPot
: public CableEquation<TDomain>
{
	public:
		// indices for unknowns
		enum {_VM_ = 0, _PHIO_};
		static const int dim = CableEquation<TDomain>::dim;

	public:
		///	constructor
		/// first function must be membrane potential (not inner!), second outer potential
	    CableEquationWithOuterPot(const char* fcts, const char* subsets);

		///	destructor
		virtual ~CableEquationWithOuterPot() {};

		///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		///	assembles stiffness part of local defect
		template <typename TElem, typename TFVGeom>
		void my_add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		/// assembles jacobian of stiffness part
		template<typename TElem, typename TFVGeom>
		void my_add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	protected:
		///	register utils
		///	\{
		void my_register_all_funcs(bool bHang);

		template <typename TElem, typename TFVGeom>
		void my_register_func();
		/// \}
};


} // namespace cable_neuron
} // namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__CABLE_DISC__CABLE_EQUATION_H__
