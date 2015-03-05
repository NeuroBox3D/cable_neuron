/*
 * channel_interface.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann
 */

#include "channel_interface.h"


namespace ug {

////////////////////////////////////////////////
// Methods for HH-Channel-Class
////////////////////////////////////////////////

template<typename TDomain>
void ChannelHH<TDomain>::
set_accuracy(number ac)
{
	m_accuracy = ac;
}


template<typename TDomain>
void ChannelHH<TDomain>::
set_conductivities(number Na, number K, number L)
{
	m_g_K = K;
	m_g_Na = Na;
	m_g_I = L;
}


template<typename TDomain>
void ChannelHH<TDomain>::
set_rev_pot(number R_Na, number R_K)
{
  m_rev_pot_Na = R_Na;
  m_rev_pot_K = R_K;

}

template<typename TDomain>
number ChannelHH<TDomain>::
vtrap(number x, number y)
{
	number vtrap;
    if (fabs(x/y) < 1e-6) {
           	   vtrap = y*(1 - x/y/2);
        }else{
               vtrap = x/(exp(x/y) - 1);
        }
    return vtrap;
}


template<typename TDomain>
void ChannelHH<TDomain>::vm_disc_available()
{
	init_attachments();
}


template<typename TDomain>
void ChannelHH<TDomain>::init_attachments()
{
	// attach attachments
	SmartPtr<Grid> spGrid = m_spVMDisc->approx_space()->domain()->grid();

	if (spGrid->has_vertex_attachment(m_MGate))
		UG_THROW("Attachment necessary (MGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_MGate);

	if (spGrid->has_vertex_attachment(m_HGate))
		UG_THROW("Attachment necessary (HGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_HGate);

	if (spGrid->has_vertex_attachment(m_NGate))
		UG_THROW("Attachment necessary (NGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_NGate);

	// create attachment accessors
	m_aaMGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, m_MGate);
	m_aaNGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, m_NGate);
	m_aaHGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, m_HGate);
}


// Methods for using gatings
template<typename TDomain>
void ChannelHH<TDomain>::init(const LocalVector& u, Edge* edge)
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	vrt_list vl;
	m_spVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge);
	for (size_t k = 0; k < vl.size(); ++k)
	{
		Vertex* vrt = vl[k];

		// update Vm
		number VM = u(m_spVMDisc->_v_, k);

		//TODO: is that right? not the same as in the Nernst version!
		// values for m gate
		/*
		number AlphaHm = 0.1 * vtrap(-(VM+40.0),10.0);
		number BetaHm =  4 * exp(-(VM+65.0)/18.0);

		// values for n gate
		number AlphaHn = 0.01*vtrap(-(VM+55.0),10.0);
		number BetaHn = 0.125*exp(-(VM+65.0)/80);

		// values for h gate
		number AlphaHh = 0.07 * exp(-(VM+65.0)/20.0);
		number BetaHh = 1.0 / (exp(-(VM+35.0)/10.0) + 1.0);
		*/

		// writing Gating-Params in attachment
		// gating param h
		number AlphaHh = 0.07*exp(-(VM + 65.0)/20.0);
		number BetaHh = 1.0/(exp(3.0-0.1*(VM  + 65.0))+1.0);

		// gating param m
		number AlphaHm;
		number AlphaHm_test = exp(2.5-0.1*(VM + 65.0))-1.0;
		if (fabs(AlphaHm_test) > m_accuracy)
			AlphaHm = (2.5 - 0.1*(VM + 65.0)) / AlphaHm_test;
		else
			AlphaHm = 1.0;

		number BetaHm = 4.0*exp(-(VM + 65.0)/18.0);


		// gating param n
		number AlphaHn;
		number AlphaHn_test;
		AlphaHn_test = exp(1.0-0.1*(VM + 65.0))-1.0;
		if (fabs(AlphaHn_test) > m_accuracy)
			AlphaHn = (0.1-0.01*(VM + 65.0)) / AlphaHn_test;
		else
			AlphaHn = 0.1;

		number BetaHn = 0.125*exp((VM + 65.0)/80.0);

		// setting initial gating params as equilibrium states
		this->m_aaHGate[vrt] = AlphaHh / (AlphaHh + BetaHh);
		this->m_aaMGate[vrt] = AlphaHm / (AlphaHm + BetaHm);
		this->m_aaNGate[vrt] = AlphaHn / (AlphaHn + BetaHn);
	}
}

template<typename TDomain>
void ChannelHH<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge)
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	vrt_list vl;
	m_spVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge);
	for (size_t k = 0; k < vl.size(); ++k)
	{
		Vertex* vrt = vl[k];
		number dt = newTime - m_spVMDisc->m_aaTime[vrt];

		// only update vertices that are not yet updated
		if (dt == 0.0) continue;

		// update Vm
		number VM = u(m_spVMDisc->_v_, k);


		// values for m gate
		number AlphaHm = 0.1 * vtrap(-(VM+40),10);
		number BetaHm =  4 * exp(-(VM+65)/18);

		// values for n gate
		number AlphaHn = 0.01*vtrap(-(VM+55),10);
		number BetaHn = 0.125*exp(-(VM+65)/80);

		// values for h gate
		number AlphaHh = 0.07 * exp(-(VM+65)/20.0);
		number BetaHh = 1.0 / (exp(-(VM+35.0)/10.0) + 1.0);

		/*
		// set gating params
		number AlphaHh = 0.07*exp(-(m_aaVm[*iter] + 65.0)/20.0);
		number BetaHh = 1.0/(exp(3.0-0.1*(m_aaVm[*iter]  + 65.0))+1.0);

		// gating param m
		number AlphaHm;
		number AlphaHm_test = exp(2.5-0.1*(m_aaVm[*iter] + 65.0))-1.0;
		if (fabs(AlphaHm_test) > m_accuracy)
			AlphaHm = (2.5 - 0.1*( m_aaVm[*iter] + 65.0)) / AlphaHm_test;
		else
			AlphaHm = 1.0;

		number BetaHm = 4.0*exp(-(m_aaVm[*iter] + 65.0)/18.0);

		// gating param n
		number AlphaHn;
		number AlphaHn_test;
		AlphaHn_test = exp(1.0-0.1*(m_aaVm[*iter] + 65.0))-1.0;
		if (fabs(AlphaHn_test) > m_accuracy)
			AlphaHn = (0.1-0.01*(m_aaVm[*iter] + 65.0)) / AlphaHn_test;
		else
			AlphaHn = 0.1;

		number BetaHn = 0.125*exp((m_aaVm[*iter] + 65.0)/80.0);
		*/

		number rate_h = ((AlphaHh/(AlphaHh+BetaHh)) - m_aaHGate[vrt]) / (1/(AlphaHh+BetaHh)) * dt;
		number rate_m = ((AlphaHm/(AlphaHm+BetaHm)) - m_aaMGate[vrt]) / (1/(AlphaHm+BetaHm)) * dt;
		number rate_n = ((AlphaHn/(AlphaHn+BetaHn)) - m_aaNGate[vrt]) / (1/(AlphaHn+BetaHn)) * dt;



		m_aaHGate[vrt] += rate_h;
		m_aaMGate[vrt] += rate_m;
		m_aaNGate[vrt] += rate_n;
		//std::cout << "VM: " << (m_aaVm[*iter]) << "h: "<< m_aaHGate[*iter] << "m: "<< m_aaMGate[*iter] <<  "n: "<< m_aaNGate[*iter] <<std::endl;
		//std::cout << "Rates: " << rate_h << " , " << rate_m << " , " << rate_n << std::endl;
	}
}


template<typename TDomain>
void ChannelHH<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	double NGate = m_aaNGate[vrt];
 	double MGate = m_aaMGate[vrt];
	double HGate = m_aaHGate[vrt];
	double VM 	 = vrt_values[VMDisc<TDomain>::_v_];

	// TODO Influx values needed
	// single channel type fluxes
	const number potassium_part_of_flux = m_g_K * pow(NGate,4) * (VM - m_rev_pot_K);
	const number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (VM + m_rev_pot_Na);
	const number leakage_part_of_flux = m_g_I * (VM + 54.4);

	number flux_value = (potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux);
	outCurrentValues.push_back(flux_value);
}


#if 0
template<typename TDomain>
void ChannelHH<TDomain>::Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux)
{
	double NGate = m_aaNGate[vrt];
	double MGate = m_aaMGate[vrt];
	double HGate = m_aaHGate[vrt];


	number Jac = (m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);

	outJFlux.push_back(Jac);

}
#endif

////////////////////////////////////////////////
// Methods for HH-Channel-Nernst-Class
////////////////////////////////////////////////

template<typename TDomain>
void ChannelHHNernst<TDomain>::
set_accuracy(double ac)
{
	m_accuracy = ac;
}


template<typename TDomain>
void ChannelHHNernst<TDomain>::
set_conductivities(number Na, number K, number L)
{
	m_g_K = K;
	m_g_Na = Na;
	m_g_I = L;
}


template<typename TDomain>
void ChannelHHNernst<TDomain>::vm_disc_available()
{
	init_attachments();
}


template<typename TDomain>
void ChannelHHNernst<TDomain>::init_attachments()
{
	// attach attachments
	SmartPtr<Grid> spGrid = m_spVMDisc->approx_space()->domain()->grid();

	if (spGrid->has_vertex_attachment(m_MGate))
		UG_THROW("Attachment necessary (MGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_MGate);

	if (spGrid->has_vertex_attachment(m_HGate))
		UG_THROW("Attachment necessary (HGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_HGate);

	if (spGrid->has_vertex_attachment(m_NGate))
		UG_THROW("Attachment necessary (NGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_NGate);

	// create attachment accessors
	m_aaMGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, m_MGate);
	m_aaNGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, m_NGate);
	m_aaHGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, m_HGate);
}



// Methods for using gatings
template<typename TDomain>
void ChannelHHNernst<TDomain>::init(const LocalVector& u, Edge* edge)
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	vrt_list vl;
	m_spVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge);
	for (size_t k = 0; k < vl.size(); ++k)
	{
		Vertex* vrt = vl[k];

		// update Vm
		number VM = u(m_spVMDisc->_v_, k);

		// writing gating params in attachments
		// gating param h
		number AlphaHh = 0.07*exp(-(VM + 65.0)/20.0);
		number BetaHh = 1.0/(exp(3.0-0.1*(VM + 65.0))+1.0);

		// gating param m
		number AlphaHm;
		number AlphaHm_test = exp(2.5-0.1*(VM + 65.0))-1.0;
		if (fabs(AlphaHm_test) > m_accuracy)
			AlphaHm = (2.5 - 0.1*(VM + 65.0)) / AlphaHm_test;
		else
			AlphaHm = 1.0;

		number BetaHm = 4.0*exp(-(VM + 65.0)/18.0);

		// gating param n
		number AlphaHn;
		number AlphaHn_test;
		AlphaHn_test = exp(1.0-0.1*(VM + 65.0))-1.0;
		if (fabs(AlphaHn_test) > m_accuracy)
			AlphaHn = (0.1-0.01*(VM + 65.0)) / AlphaHn_test;
		else
			AlphaHn = 0.1;

		number BetaHn = 0.125*exp((VM + 65.0)/80.0);

		// setting initial gating params as equilibrium states
		this->m_aaHGate[vrt] = AlphaHh / (AlphaHh + BetaHh);
		this->m_aaMGate[vrt] = AlphaHm / (AlphaHm + BetaHm);
		this->m_aaNGate[vrt] = AlphaHn / (AlphaHn + BetaHn);
	}
}

template<typename TDomain>
void ChannelHHNernst<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge)
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	vrt_list vl;
	m_spVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge);
	for (size_t k = 0; k < vl.size(); ++k)
	{
		Vertex* vrt = vl[k];
		number dt = newTime - m_spVMDisc->m_aaTime[vrt];

		// only update vertices that are not yet updated
		if (dt == 0.0) continue;

		// update Vm
		number VM = u(m_spVMDisc->_v_, k);

		// set new gating states
		// gating param h
		number AlphaHh = 0.07*exp(-(VM + 65.0)/20.0);
		number BetaHh = 1.0/(exp(3.0-0.1*(VM + 65.0))+1.0);

		// gating param m
		number AlphaHm;
		number AlphaHm_test = exp(2.5-0.1*(VM + 65.0))-1.0;
		if (fabs(AlphaHm_test) > m_accuracy)
			AlphaHm = (2.5 - 0.1*(VM + 65.0)) / AlphaHm_test;
		else
			AlphaHm = 1.0;

		number BetaHm = 4.0*exp(-(VM + 65.0)/18.0);

		// gating param n
		number AlphaHn;
		number AlphaHn_test;
		AlphaHn_test = exp(1.0-0.1*(VM + 65.0))-1.0;
		if (fabs(AlphaHn_test) > m_accuracy)
			AlphaHn = (0.1-0.01*(VM + 65.0)) / AlphaHn_test;
		else
			AlphaHn = 0.1;

		number BetaHn = 0.125*exp((VM + 65.0)/80.0);

		number rate_h = ((AlphaHh/(AlphaHh+BetaHh)) - m_aaHGate[vrt]) / (1.0/(AlphaHh+BetaHh)) * dt;
		number rate_m = ((AlphaHm/(AlphaHm+BetaHm)) - m_aaMGate[vrt]) / (1.0/(AlphaHm+BetaHm)) * dt;
		number rate_n = ((AlphaHn/(AlphaHn+BetaHn)) - m_aaNGate[vrt]) / (1.0/(AlphaHn+BetaHn)) * dt;

		m_aaHGate[vrt] += rate_h;
		m_aaMGate[vrt] += rate_m;
		m_aaNGate[vrt] += rate_n;
	}
}


template<typename TDomain>
void ChannelHHNernst<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments out of Vertex
	number NGate = m_aaNGate[vrt];
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];
	number v 	 = vrt_values[VMDisc<TDomain>::_v_];
	number k 	 = vrt_values[VMDisc<TDomain>::_k_];
	number na 	 = vrt_values[VMDisc<TDomain>::_na_];


	UG_ASSERT(m_spVMDisc, "Channel has not been assigned a vmDisc object yet!");
	const number helpV = 1e3*(m_R*m_T)/m_F;
	number potassium_nernst_eq 	= helpV*(std::log(m_spVMDisc->k_out/k));
	number sodium_nernst_eq	 	= helpV*(std::log(m_spVMDisc->na_out/na));

	// single channel ion fluxes
	number potassium_part_of_flux = m_g_K * pow(NGate,4) * (v - potassium_nernst_eq);
	number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (v - sodium_nernst_eq);
	number leakage_part_of_flux = m_g_I * (v - (-54.4));

	//std::cout << "potassium_part_of_flux: " << potassium_part_of_flux << std::endl;
	//std::cout << "sodium_part_of_flux: " << sodium_part_of_flux << std::endl;

	//std::cout << " n: " << NGate << " m: "<< MGate << " h: " << HGate << std::endl;

	outCurrentValues.push_back(potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux);
	outCurrentValues.push_back(potassium_part_of_flux/m_F);
	outCurrentValues.push_back(sodium_part_of_flux/m_F);

	//std::cout << "pot: " << potassium_part_of_flux << " sod : " << sodium_part_of_flux << " leak: " << leakage_part_of_flux << std::endl;
	//std::cout << "leakeage term: " << leakage_part_of_flux << std::endl;
	//std::cout << "outCurrentValues: " << outCurrentValues[0] << ", " << outCurrentValues[1] << ", " << outCurrentValues[2] << ", " << std::endl;
	//std::cout << "end ionic current" << std::endl;
}


#if 0
template<typename TDomain>
void ChannelHHNernst<TDomain>::Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux)
{
	number NGate = m_aaNGate[vrt];
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];
	number k 	 = vrt_values[VMDisc<TDomain>::_k_];
	number na 	 = vrt_values[VMDisc<TDomain>::_na_];

	UG_ASSERT(m_spVMDisc, "Channel has not been assigned a vmDisc object yet!");
	const number helpV = (m_R*m_T)/m_F;
	const number potassium_nernst_eq_dK 	=  helpV * (-m_spVMDisc->k_out/k)*0.18; //helpV * (-K_out/pow(u(_K_,co),2));
	const number sodium_nernst_eq_dNa		=  -helpV * (-m_spVMDisc->na_out/na)*0.003;

	outJFlux.push_back(m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);
	outJFlux.push_back(sodium_nernst_eq_dNa);
	outJFlux.push_back(potassium_nernst_eq_dK);

	//std::cout << "outJFlux: " << outJFlux[0] << ", " << outJFlux[1] << ", " << outJFlux[2] << ", " << std::endl;
}
#endif

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class IChannel<Domain1d>;
	template class ChannelHH<Domain1d>;
	template class ChannelHHNernst<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class IChannel<Domain2d>;
	template class ChannelHH<Domain2d>;
	template class ChannelHHNernst<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class IChannel<Domain3d>;
	template class ChannelHH<Domain3d>;
	template class ChannelHHNernst<Domain3d>;
#endif


} // namespace ug
