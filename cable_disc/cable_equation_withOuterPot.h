/*
 * cable_equation.h
 *
 *  Created on: 2019-01-09
 *  Author: mbreit
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
