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

#ifndef UG__PLUGINS__CABLE_NEURON__CABLE_DISC__IMPLICIT_ACTIVE_CABLE_DISC_BASE_H
#define UG__PLUGINS__CABLE_NEURON__CABLE_DISC__IMPLICIT_ACTIVE_CABLE_DISC_BASE_H

#include "common/types.h"  // for number
#include "common/util/smart_pointer.h"  // for SmartPtr
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"  // for IElemDisc
#include "lib_disc/spatial_disc/user_data/data_import.h"  // for DataImport
#include "lib_disc/spatial_disc/user_data/user_data.h"  // for CplUserData
#include "lib_disc/spatial_disc/user_data/user_function.h"  // for IFunction
#include "lib_grid/common_attachments.h"  // for ANumber
#include "lib_grid/grid/grid.h"  // for Grid::AttachmentAccessor
#include "lib_grid/grid/grid_base_objects.h"  // for Vertex

#include "../util/diam_attachment_handler.h"  // for DiamAttachmentHandler


namespace ug {
namespace cable_neuron {


template <typename TDomain>
class ImplicitActiveCableDiscBase
: public IElemDisc<TDomain>
{
	private:
		///	base class type
		typedef IElemDisc<TDomain> base_type;

	public:
		///	world dimension
		static const int dim = base_type::dim;

	protected:
		///	abbreviation for the local solution components
		static const size_t _vm_ = 0;
		static const size_t _h_ = 1;
		static const size_t _m_ = 2;
		static const size_t _n_ = 3;

	public:
		///	constructor
		ImplicitActiveCableDiscBase(const char* functions, const char* subsets);

		///	@brief set specific capacitance of the membrane
		///	\{
		void set_spec_cap(number val);
		void set_spec_cap(SmartPtr<CplUserData<number, dim> > fct);
#ifdef UG_FOR_LUA
		void set_spec_cap(const char* fctName);
#endif
		///	\}

		///	@brief set specific resistivity of the cytosol
		///	\{
		void set_spec_res(number val);
		void set_spec_res(SmartPtr<CplUserData<number, dim> > fct);
#ifdef UG_FOR_LUA
		void set_spec_res(const char* fctName);
#endif
		///	\}

		///	@brief set membrane conductances for K+, Na+ and leakage
		///	\{
		void set_conductances(number gK, number gNa, number gL);
		void set_conductances
		(
			SmartPtr<CplUserData<number, dim> > gKFct,
			SmartPtr<CplUserData<number, dim> > gNaFct,
			SmartPtr<CplUserData<number, dim> > gLFct
		);
#ifdef UG_FOR_LUA
		void set_conductances(const char* gKFctName, const char* gNaFctName, const char* gLFctName);
#endif
		///	\}

		///	@brief set K+, Na+ and leakage reversal potentials
		///	\{
		void set_rev_pot(number eK, number eNa, number eL);
		void set_rev_pot
		(
			SmartPtr<CplUserData<number, dim> > eKFct,
			SmartPtr<CplUserData<number, dim> > eNaFct,
			SmartPtr<CplUserData<number, dim> > eLFct
		);
#ifdef UG_FOR_LUA
		void set_rev_pot(const char* eKFctName, const char* eNaFctName, const char* eLFctName);
#endif
		///	\}


		/// set constant diameter
		void set_diameter(number d);

		/// set injecting electrode per function
		void set_injection(IFunction<number>& functor);

		/// set temperature in units of K
		void set_temperature(number k);

		/// whether or not temperature dependency is to be enabled
		void enable_temperature_dependency(bool enable);


	public:
		/// @copydoc IElemDisc::approximation_space_changed()
		virtual void approximation_space_changed();

		///	returns if hanging nodes are needed
		virtual bool use_hanging() const;


	protected:
		/// ideal gas constant
		const number m_R;

		/// temperature
		number m_T;

		/// Faraday's constant
		const number m_F;

		///	specific membrane capacitance
		DataImport<number, dim> m_spec_cap;

		/// specific resistivity of the cytosol
		DataImport<number, dim> m_spec_res;

		/// membrane conductances for K+, Na+ and leakage
		/// \{
		DataImport<number, dim> m_gK;
		DataImport<number, dim> m_gNa;
		DataImport<number, dim> m_gL;
		/// \}

		/// reversal potentials K+, Na+ and leakage
		/// \{
		DataImport<number, dim> m_eK;
		DataImport<number, dim> m_eNa;
		DataImport<number, dim> m_eL;
		/// \}

		/// injection electrode function
		IFunction<number>* m_Injection;

		ANumber m_aDiameter;									///< dendritic diameter attachment
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaDiameter;	///< dendritic diameter attachment accessor
		DiamAttachmentHandler m_dah;							///< handler for multigrid usage of diameter attachment
		number m_constDiam;										///< constant diameter (if set)
		bool m_bConstDiamSet;									///< whether const diameter is set

		// temperature dependency flag
		bool m_bTempDep;
};


} // end namespace cable_neuron
} // end namespace ug


#endif // UG__PLUGINS__CABLE_NEURON__CABLE_DISC__IMPLICIT_ACTIVE_CABLE_DISC_BASE_H
