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

#include "common/error.h"  // for UG_COND_THROW
#include "implicit_active_cable_disc_base.h"
#include "lib_disc/domain.h"  // for Domain3d
#include "lib_disc/spatial_disc/user_data/const_user_data.h"  // for ConstUserNumber
#include "lib_grid/global_attachments.h"  // for GlobalAttachment
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"  // for LuaUserDataFactory
#endif

namespace ug {
namespace cable_neuron {



template<typename TDomain>
ImplicitActiveCableDiscBase<TDomain>::
ImplicitActiveCableDiscBase(const char* functions, const char* subsets)
 : IElemDisc<TDomain>(functions,subsets),
   m_R(8.314), m_T(310.0), m_F(96485.0),
   m_Injection(NULL),
   m_aDiameter(GlobalAttachments::attachment<ANumber>("diameter")),
   m_constDiam(1e-6),
   m_bConstDiamSet(false),
   m_bTempDep(false)
{
	// check number of functions
	UG_COND_THROW(this->num_fct() < 4, "Wrong number of functions: "
		"ImplicitActiveCableDiscBase needs a minimum of 4 symbolic functions "
		"(1 for V_m and 3 for gating variables).");

	// register imports
	this->register_import(m_spec_cap);
	this->register_import(m_spec_res);
	this->register_import(m_gK);
	this->register_import(m_gNa);
	this->register_import(m_gL);
	this->register_import(m_eK);
	this->register_import(m_eNa);
	this->register_import(m_eL);
}



template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_spec_cap(number val)
{
	if (val == 0.0)
		set_spec_cap(SmartPtr<CplUserData<number, dim> >());
	else
		set_spec_cap(make_sp(new ConstUserNumber<dim>(val)));
}

template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_spec_cap(SmartPtr<CplUserData<number, dim> > fct)
{
	m_spec_cap.set_data(fct);
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_spec_cap(const char* fctName)
{
	set_spec_cap(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif



template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_spec_res(number val)
{
	if (val == 0.0)
		set_spec_res(SmartPtr<CplUserData<number, dim> >());
	else
		set_spec_res(make_sp(new ConstUserNumber<dim>(val)));
}

template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_spec_res(SmartPtr<CplUserData<number, dim> > fct)
{
	m_spec_res.set_data(fct);
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_spec_res(const char* fctName)
{
	set_spec_res(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif



template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_conductances(number gK, number gNa, number gL)
{
	if (gK == 0.0)
		m_gK.set_data(SmartPtr<CplUserData<number, dim> >());
	else
		m_gK.set_data(make_sp(new ConstUserNumber<dim>(gK)));

	if (gNa == 0.0)
		m_gNa.set_data(SmartPtr<CplUserData<number, dim> >());
	else
		m_gNa.set_data(make_sp(new ConstUserNumber<dim>(gNa)));

	if (gL == 0.0)
		m_gL.set_data(SmartPtr<CplUserData<number, dim> >());
	else
		m_gL.set_data(make_sp(new ConstUserNumber<dim>(gL)));
}

template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_conductances
(
	SmartPtr<CplUserData<number, dim> > gKFct,
	SmartPtr<CplUserData<number, dim> > gNaFct,
	SmartPtr<CplUserData<number, dim> > gLFct
)
{
	m_gK.set_data(gKFct);
	m_gNa.set_data(gNaFct);
	m_gL.set_data(gLFct);
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_conductances(const char* gKFctName, const char* gNaFctName, const char* gLFctName)
{
	m_gK.set_data(LuaUserDataFactory<number, dim>::create(gKFctName));
	m_gNa.set_data(LuaUserDataFactory<number, dim>::create(gNaFctName));
	m_gL.set_data(LuaUserDataFactory<number, dim>::create(gLFctName));
}
#endif



template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_rev_pot(number eK, number eNa, number eL)
{
	if (eK == 0.0)
		m_eK.set_data(SmartPtr<CplUserData<number, dim> >());
	else
		m_eK.set_data(make_sp(new ConstUserNumber<dim>(eK)));

	if (eNa == 0.0)
		m_eNa.set_data(SmartPtr<CplUserData<number, dim> >());
	else
		m_eNa.set_data(make_sp(new ConstUserNumber<dim>(eNa)));

	if (eL == 0.0)
		m_eL.set_data(SmartPtr<CplUserData<number, dim> >());
	else
		m_eL.set_data(make_sp(new ConstUserNumber<dim>(eL)));
}

template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_rev_pot
(
	SmartPtr<CplUserData<number, dim> > eKFct,
	SmartPtr<CplUserData<number, dim> > eNaFct,
	SmartPtr<CplUserData<number, dim> > eLFct
)
{
	m_eK.set_data(eKFct);
	m_eNa.set_data(eNaFct);
	m_eL.set_data(eLFct);
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_rev_pot(const char* eKFctName, const char* eNaFctName, const char* eLFctName)
{
	m_eK.set_data(LuaUserDataFactory<number, dim>::create(eKFctName));
	m_eNa.set_data(LuaUserDataFactory<number, dim>::create(eNaFctName));
	m_eL.set_data(LuaUserDataFactory<number, dim>::create(eLFctName));
}
#endif



template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_injection(IFunction<number>& functor)
{
	m_Injection = &functor;
}


template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::set_temperature(number k)
{
	m_T = k;
}



template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::
set_diameter(const number d)
{
	m_constDiam = d;
	m_bConstDiamSet = true;
}



template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::enable_temperature_dependency(bool enable)
{
	m_bTempDep = enable;
}



template<typename TDomain>
void ImplicitActiveCableDiscBase<TDomain>::approximation_space_changed()
{
	// handle diameter attachment
	SmartPtr<MultiGrid> grid = this->approx_space()->domain()->grid();
	if (!grid->has_attachment<Vertex>(m_aDiameter))
		grid->attach_to_vertices_dv(m_aDiameter, m_constDiam);
	else
	{
		if (m_bConstDiamSet)
		{
			UG_LOG("HINT: Even though you have explicitly set a constant diameter to the domain\n"
				   "      this discretization will use the diameter information attached to the grid\n"
				   "      you specified.\n");
		}
	}

	// this will distribute the attachment values to the whole grid
	m_dah.set_attachment(m_aDiameter);
	m_dah.set_grid(grid);

	// create accessor
	m_aaDiameter = Grid::AttachmentAccessor<Vertex, ANumber>(*grid, m_aDiameter);
}


template<typename TDomain>
bool ImplicitActiveCableDiscBase<TDomain>::
use_hanging() const
{
	// As this is basically a 1D discretization,
	// there will never be any hanging nodes.
	return false;
}




// explicit template instantiations
#ifdef UG_DIM_3
	template class ImplicitActiveCableDiscBase<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
