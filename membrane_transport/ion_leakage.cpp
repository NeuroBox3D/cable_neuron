/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Pascal Gottmann
 * Creation date: 2015-08-24
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

#include "ion_leakage.h"


namespace ug {
namespace cable_neuron {

template<typename TDomain>
IonLeakage<TDomain>::IonLeakage(const char* functions, const char* subsets)
try : ICableMembraneTransport<TDomain>(functions, subsets),
m_perm(1.0), m_bOhmic(false), m_cond(0.0), m_revPot(0.0), m_leaking_fct(std::string("")), m_lfInd(0),
m_flux_at_rest(0.0), m_conc_in_rest(5e-5), m_conc_out_rest(1.5), m_vm_rest(-0.065), m_bPermSetByRestingConds(false),
m_valency(2)
{
    // check that there is one given function at most
    std::vector<std::string> function_tokens = TokenizeString(functions, ',');
    UG_COND_THROW(function_tokens.size() > 1, "Only one function argument is supported for ion leakage.");

    // set function as leaking function
    if (function_tokens.size())
        m_leaking_fct = function_tokens[0];
}
UG_CATCH_THROW("Error in ChannelHH initializer list.");

template<typename TDomain>
IonLeakage<TDomain>::IonLeakage
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : ICableMembraneTransport<TDomain>(functions, subsets),
m_perm(1.0), m_bOhmic(false), m_cond(0.0), m_revPot(0.0), m_leaking_fct(std::string("")), m_lfInd(0),
m_flux_at_rest(0.0), m_conc_in_rest(5e-5), m_conc_out_rest(1.5), m_vm_rest(-0.065), m_bPermSetByRestingConds(false),
m_valency(2)
{
    // check that there is one given function at most
    UG_COND_THROW(functions.size() > 1, "Only one function argument is supported for ion leakage.");

    // set function as leaking function
    if (functions.size())
        m_leaking_fct = functions[0];
}
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
set_perm(number flux_at_rest, number conc_in_rest, number conc_out_rest, number vm_rest, int valency)
{
	// save parameters for evaluation when CableEquation avail.
	m_flux_at_rest = flux_at_rest;
	m_conc_in_rest = conc_in_rest;
	m_conc_out_rest = conc_out_rest;
	m_vm_rest = vm_rest;
	m_valency = valency;
	m_bPermSetByRestingConds = true;
}


template<typename TDomain>
void IonLeakage<TDomain>::
set_perm(number perm)
{
	m_perm = perm;
	m_bPermSetByRestingConds = false;
}


template<typename TDomain>
void IonLeakage<TDomain>::set_valency(int v)
{
	m_valency = v;
}


template<typename TDomain>
void IonLeakage<TDomain>::set_ohmic(bool b)
{
	m_bOhmic = b;
}


template<typename TDomain>
void IonLeakage<TDomain>::set_cond(number g)
{
	m_cond = g;
}


template<typename TDomain>
void IonLeakage<TDomain>::set_rev_pot(number e)
{
	m_revPot = e;
}


template<typename TDomain>
void IonLeakage<TDomain>::ce_obj_available()
{
	try
	{
		m_lfInd = this->m_pCE->function_pattern()->fct_id_by_name(m_leaking_fct.c_str());
	}
	UG_CATCH_THROW("Leaking quantity for Leakage machanism not defined in functions of CableEquation.")

	if (m_bPermSetByRestingConds)
	{

		// set permeability coeff
		const number& R = m_pCE->R;
		const number& F = m_pCE->F;
		const number& T = m_pCE->temperature();
		const int z = m_valency;

		if (fabs(m_vm_rest) < 1e-8) m_perm = m_flux_at_rest / ((m_conc_out_rest - m_conc_in_rest) - z*F/(2*R*T) * (m_conc_out_rest + m_conc_in_rest)*m_vm_rest);
		else m_perm = m_flux_at_rest / (z*F/(R*T) * m_vm_rest * (m_conc_out_rest - m_conc_in_rest*exp(z*F/(R*T)*m_vm_rest)) / (1.0 - exp(z*F/(R*T)*m_vm_rest)));

		// check that permeability is positive (else: modeling error!)
		UG_COND_THROW(m_perm < 0, "The permeability coefficient of your ion leakage term is negative.\n"
					  "This is not allowed since this mechanism only represents passive fluxes. You may "
					  "want to consider adding an active mechanism (pump) to your model.");
	}
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
void IonLeakage<TDomain>::current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	const number& vm = vrt_values[CableEquation<TDomain>::_v_];

	// ohmic currents
	if (m_bOhmic)
	{
		const number leak = m_cond * (vm - m_revPot);
		outCurrentValues.push_back(leak);
		outCurrentValues.push_back(leak / (m_valency*m_pCE->F));
		return;
	}

	// GHK currents

	// getting attachments for vertex
	const number& conc_in = vrt_values[m_lfInd];
	const number& conc_out = this->m_pCE->conc_out(m_lfInd);

	const number& R = m_pCE->R;
	const number& F = m_pCE->F;
	const number& T = m_pCE->temperature();
	const int z = m_valency;

	number leak;

	if (fabs(vm) < 1e-8) leak = -m_perm * ((conc_out - conc_in) - z*F/(2*R*T) * (conc_out + conc_in)*vm);
		else leak = m_perm * z*F/(R*T) * vm * (conc_out - conc_in*exp(z*F/(R*T)*vm)) / (1.0 - exp(z*F/(R*T)*vm));

	outCurrentValues.push_back(leak * z*F);
	outCurrentValues.push_back(leak);
}


template<typename TDomain>
number IonLeakage<TDomain>::lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values)
{
	// ohmic currents
	if (m_bOhmic)
		return m_cond;

	// GHK currents
	const number& vm = vrt_values[CableEquation<TDomain>::_v_];
	const number& conc_in = vrt_values[m_lfInd];
	const number& conc_out = this->m_pCE->conc_out(m_lfInd);

	const number& R = m_pCE->R;
	const number& F = m_pCE->F;
	const number& T = m_pCE->temperature();
	const int z = m_valency;

	// flux derived from Goldman-Hodgkin-Katz equation,
	number leakDeriv;

	// near V_m == 0: approximate by first order Taylor to avoid relative errors and div-by-0
	if (fabs(vm) < 1e-8)
		leakDeriv = m_perm * z*F/(2*R*T) * (conc_out + conc_in);
	else
	{
		number in = z*F/(R*T)*vm;
		number ex = exp(in);
		leakDeriv = m_perm * 2*F/(R*T) * (conc_out*(1.0 - (1.0 - in)*ex) + conc_in*(ex - (1.0 + in))*ex)
						/ ((1.0 - ex) * (1.0 - ex));
	}

	return leakDeriv * z*F;
}


template<typename TDomain>
void IonLeakage<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing CableEquation fct indices which this channel writes to
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_);
	this->m_vWFctInd.push_back(m_lfInd);
}



////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////
#ifdef UG_DIM_3
	template class IonLeakage<Domain3d>;
#endif


} // namespace cable_neuron
} // namespace ug
