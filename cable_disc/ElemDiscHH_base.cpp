/*
 * Created on: 13.06.2014
 * 		Author: Pascal Gottmann
 * ElemDiscHH_base.h
 *
 * Based on
 * convection_diffusion.h
 * from andreasvogel
 */

#include "../cable_disc/ElemDiscHH_base.h"

#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug {
namespace cable_neuron {


////////////////////////////////////////////////////////////////////////////////
//	user data
////////////////////////////////////////////////////////////////////////////////

//////// Diffusion

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user)
{
	m_imDiffusion.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::set_diffusion(number val)
{
	if(val == 0.0) set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >());
	else set_diffusion(make_sp(new ConstUserMatrix<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::set_diffusion(const char* fctName)
{
	set_diffusion(LuaUserDataFactory<MathMatrix<dim,dim>, dim>::create(fctName));
}
#endif


//////// Velocity

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imVelocity.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::set_velocity(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_velocity(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_velocity(const char* fctName)
{
	set_velocity(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
#endif

//////// Flux

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imFlux.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::set_flux(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_flux(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_flux(const char* fctName)
{
	set_flux(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
#endif

//////// Reaction Rate

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_rate(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionRate.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_rate(number val)
{
	if(val == 0.0) set_reaction_rate(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_rate(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_rate(const char* fctName)
{
	set_reaction_rate(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Reaction Rate Explicit

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_rate_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionRateExpl.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_rate_explicit(number val)
{
	if(val == 0.0) set_reaction_rate_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_rate_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_rate_explicit(const char* fctName)
{
	set_reaction_rate_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Reaction

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReaction.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction(number val)
{
	if(val == 0.0) set_reaction(SmartPtr<CplUserData<number, dim> >());
	else set_reaction(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction(const char* fctName)
{
	set_reaction(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Reaction Explicit

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionExpl.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_explicit(number val)
{
	if(val == 0.0) set_reaction_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_reaction_explicit(const char* fctName)
{
	set_reaction_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif


//////// Source

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_source(SmartPtr<CplUserData<number, dim> > user)
{
	m_imSource.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_source(number val)
{
	if(val == 0.0) set_source(SmartPtr<CplUserData<number, dim> >());
	else set_source(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_source(const char* fctName)
{
	set_source(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Source explicit

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_source_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imSourceExpl.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_source_explicit(number val)
{
	if(val == 0.0) set_source_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_source_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_source_explicit(const char* fctName)
{
	set_source_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif


//////// Vector Source

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_vector_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imVectorSource.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::set_vector_source(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_vector_source(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_vector_source(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::set_vector_source(const char* fctName)
{
	set_vector_source(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
#endif

//////// Mass Scale

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_spec_capa(SmartPtr<CplUserData<number, dim> > user)
{
	m_spec_capa.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_spec_capa(number val)
{
	if(val == 0.0) set_spec_capa(SmartPtr<CplUserData<number, dim> >());
	else set_spec_capa(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_spec_capa(const char* fctName)
{
	set_spec_capa(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Mass

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_mass(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMass.set_data(user);
}

template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_mass(number val)
{
	if(val == 0.0) set_mass(SmartPtr<CplUserData<number, dim> >());
	else set_mass(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ElemDiscHH_Base<TDomain>::
set_mass(const char* fctName)
{
	set_mass(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

////////////////////////////////////////////////////////////////////////////////
//	Exports
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
typename ElemDiscHH_Base<TDomain>::NumberExport
ElemDiscHH_Base<TDomain>::
value() {return m_exValue;}


template <typename TDomain>
typename ElemDiscHH_Base<TDomain>::GradExport
ElemDiscHH_Base<TDomain>::
gradient() {return m_exGrad;}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ElemDiscHH_Base<TDomain>::
ElemDiscHH_Base(const char* functions, const char* subsets)
 : IElemDisc<TDomain>(functions,subsets),
   m_exValue(new DataExport<number, dim>(functions)),
   m_exGrad(new DataExport<MathVector<dim>, dim>(functions))
{
//	check number of functions
	if(this->num_fct() < 4)
		UG_THROW("Wrong number of functions: The ElemDiscHH 'ElemDiscHH_Base'"
					   " needs a minimum of "<<4<<" symbolic function (1 for VM and 3 for Gate-Vars)");

//	register imports
	this->register_import(m_imDiffusion);
	this->register_import(m_imVelocity);
	this->register_import(m_imFlux);
	this->register_import(m_imReactionRate);
	this->register_import(m_imReaction);
	this->register_import(m_imReactionRateExpl);
	this->register_import(m_imReactionExpl);
	this->register_import(m_imSourceExpl);
	this->register_import(m_imSource);
	this->register_import(m_imVectorSource);
	this->register_import(m_spec_capa);
	this->register_import(m_imMass);

	//m_imMassScale.set_mass_part();
	m_imMass.set_mass_part();
	m_imSource.set_rhs_part();
	m_imVectorSource.set_rhs_part();
	m_imSourceExpl.set_expl_part();
	m_imReactionExpl.set_expl_part();
	m_imReactionRateExpl.set_expl_part();


}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ElemDiscHH_Base<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ElemDiscHH_Base<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ElemDiscHH_Base<Domain3d>;
#endif

} // namespace cable_neuron
} // namespace ug