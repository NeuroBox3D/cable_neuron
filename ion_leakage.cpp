/*
 * leakage.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann, mbreit
 */

#include "ion_leakage.h"


namespace ug {
namespace cable {

template<typename TDomain>
IonLeakage<TDomain>::IonLeakage(const char* functions, const char* subsets)
try : IChannel<TDomain>(functions, subsets),
m_perm(1.0e-6), m_leaking_fct(std::string("v")), m_lfInd(0),
m_flux_at_rest(0.0), m_conc_in_rest(5e-5), m_conc_out_rest(1.5), m_vm_rest(-65) {}
UG_CATCH_THROW("Error in ChannelHH initializer list.");

template<typename TDomain>
IonLeakage<TDomain>::IonLeakage
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : IChannel<TDomain>(functions, subsets),
m_perm(1.0e-6), m_leaking_fct(std::string("v")), m_lfInd(0),
m_flux_at_rest(0.0), m_conc_in_rest(5e-5), m_conc_out_rest(1.5), m_vm_rest(-65) {}
UG_CATCH_THROW("Error in ChannelHH initializer list.");


template<typename TDomain>
std::string IonLeakage<TDomain>::
name()
{
	return std::string("ionic leakage");
}


template<typename TDomain>
void IonLeakage<TDomain>::
set_leaking_quantity(const std::string& lq)
{
	m_leaking_fct = lq;
}


template<typename TDomain>
void IonLeakage<TDomain>::
set_perm(number flux_at_rest, number conc_in_rest, number conc_out_rest, number vm_rest)
{
	// save parameters for evaluation when VMDisc avail.
	m_flux_at_rest = flux_at_rest;
	m_conc_in_rest = conc_in_rest;
	m_conc_out_rest = conc_out_rest;
	m_vm_rest = 1e-3*vm_rest;	// scale from mV to V
}


template<typename TDomain>
void IonLeakage<TDomain>::vm_disc_available()
{
	try
	{
		m_lfInd = this->m_pVMDisc->function_pattern()->fct_id_by_name(m_leaking_fct.c_str());
	}
	UG_CATCH_THROW("Leaking quantity for Leakage machanism not defined in functions of VMDisc.")

	// set permeability coeff
	const number& R = m_pVMDisc->R;
	const number& F = m_pVMDisc->F;
	const number& T = m_pVMDisc->temperature();

	if (fabs(m_vm_rest) < 1e-8) m_perm = m_flux_at_rest / ((m_conc_out_rest - m_conc_in_rest) - F/(R*T) * (m_conc_out_rest + m_conc_in_rest)*m_vm_rest);
	else m_perm = m_flux_at_rest / (2*F/(R*T) * m_vm_rest * (m_conc_out_rest - m_conc_in_rest*exp(2*F/(R*T)*m_vm_rest)) / (1.0 - exp(2*F/(R*T)*m_vm_rest)));

	// check that permeability is positive (else: modeling error!)
	UG_COND_THROW(m_perm < 0, "The permeability coefficient of your ion leakage term is negative.\n"
				  "This is not allowed since this mechanism only represents passive fluxes. You may "
				  "want to consider adding an active mechanism (pump) to your model.");
}


template<typename TDomain>
void IonLeakage<TDomain>::init_attachments()
{

}


// Methods for using gatings
template<typename TDomain>
void IonLeakage<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	// nothing to do
}

template<typename TDomain>
void IonLeakage<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{
	// nothing to do
}


template<typename TDomain>
void IonLeakage<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	const number& conc_in 	 = vrt_values[m_lfInd];
	const number& conc_out 	 = this->m_pVMDisc->conc_out(m_lfInd);
	const number& VM	 	 = 1e-3*vrt_values[VMDisc<TDomain>::_v_]; // scale to V!

	const number& R = m_pVMDisc->R;
	const number& F = m_pVMDisc->F;
	const number& T = m_pVMDisc->temperature();

	number leak;

	if (fabs(VM) < 1e-8) leak = -m_perm * ((conc_out - conc_in) - F/(R*T) * (conc_out + conc_in)*VM);
		else leak = m_perm * 2*F/(R*T) * VM * (conc_out - conc_in*exp(2*F/(R*T)*VM)) / (1.0 - exp(2*F/(R*T)*VM));

	outCurrentValues.push_back(leak);
}


template<typename TDomain>
void IonLeakage<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing VMDisc fct indices which this channel writes to
	this->m_vWFctInd.push_back(m_lfInd);
}



////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class IonLeakage<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class IonLeakage<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class IonLeakage<Domain3d>;
#endif


} // namespace cable
} // namespace ug
