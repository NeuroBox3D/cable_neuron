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

template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::
set_accuracy(number ac)
{
	m_accuracy = ac;
}


template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::
set_conductivities(number Na, number K, number L)
{
	m_g_K = K;
	m_g_Na = Na;
	m_g_I = L;
}


template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::
set_rev_pot(number R_Na, number R_K)
{
  m_rev_pot_Na = R_Na;
  m_rev_pot_K = R_K;

}

template<typename TDomain, typename TAlgebra>
number ChannelHH<TDomain, TAlgebra>::
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


// Methods for using gatings
template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::init(number time, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct)
{
	// attach attachments
	if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->m_MGate))
		UG_THROW("Attachment necessary (MGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->m_MGate);

	if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->m_HGate))
		UG_THROW("Attachment necessary (HGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
		spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->m_HGate);

	if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->m_NGate))
		UG_THROW("Attachment necessary (NGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->m_NGate);

	// create attachment accessors
	this->m_aaMGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGridFct->domain()->grid(), this->m_MGate);
	this->m_aaNGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGridFct->domain()->grid(), this->m_NGate);
	this->m_aaHGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGridFct->domain()->grid(), this->m_HGate);


	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(spGridFct->domain()->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");


	std::vector<DoFIndex> multInd;

	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// update Vm
			const size_t VM_fct_ind = spGridFct->fct_id_by_name("v");

			spGridFct->dof_indices(*iter, VM_fct_ind, multInd);
			UG_ASSERT(multInd.size() == 1, "multi-index has size != 1");

			number VM = DoFRef(*spGridFct, multInd[0]);

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
			this->m_aaHGate[*iter] = AlphaHh / (AlphaHh + BetaHh);
			this->m_aaMGate[*iter] = AlphaHm / (AlphaHm + BetaHm);
			this->m_aaNGate[*iter] = AlphaHn / (AlphaHn + BetaHn);
		}
	}
}

template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::update_gating(number newTime, ConstSmartPtr<typename TAlgebra::vector_type> uOld)
{
	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(m_vmDisc->approx_space()->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");


	ConstSmartPtr<DoFDistribution> dd = m_vmDisc->approx_space()->dof_distribution(GridLevel::TOP);
	std::vector<DoFIndex> multInd;

	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// update Vm
			const size_t VM_fct_ind = dd->fct_id_by_name("VM");
			dd->dof_indices(*iter, VM_fct_ind, multInd);
			UG_ASSERT(multInd.size() == 1, "multi-index has size != 1");
			number VM = DoFRef(*uOld, multInd[0]);

	        // values for m gate
	        number AlphaHm = 0.1 * vtrap(-(VM+40),10);
	        number BetaHm =  4 * exp(-(VM+65)/18);

	        // values for n gate
	        number AlphaHn = 0.01*vtrap(-(VM+55),10);
	        number BetaHn = 0.125*exp(-(VM+65)/80);

	        // values for h gate
	        number AlphaHh = 0.07 * exp(-(VM+65)/20.0);
	        number BetaHh = 1.0 / (exp(-(VM+35.0)/10.0) + 1.0);

			/*//set gatting Params
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

			number BetaHn = 0.125*exp((m_aaVm[*iter] + 65.0)/80.0);*/

			// needs import later
			number dt = 0.01;

			number rate_h = ((AlphaHh/(AlphaHh+BetaHh)) - m_aaHGate[*iter]) / (1/(AlphaHh+BetaHh)) *dt;
			number rate_m = ((AlphaHm/(AlphaHm+BetaHm)) - m_aaMGate[*iter]) / (1/(AlphaHm+BetaHm)) *dt;
			number rate_n = ((AlphaHn/(AlphaHn+BetaHn)) - m_aaNGate[*iter]) / (1/(AlphaHn+BetaHn)) *dt;



			m_aaHGate[*iter] += rate_h;
			m_aaMGate[*iter] += rate_m;
			m_aaNGate[*iter] += rate_n;
			//std::cout << "VM: " << (m_aaVm[*iter]) << "h: "<< m_aaHGate[*iter] << "m: "<< m_aaMGate[*iter] <<  "n: "<< m_aaNGate[*iter] <<std::endl;
			//std::cout << "Rates: " << rate_h << " , " << rate_m << " , " << rate_n << std::endl;
		}
	}
}


template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	double NGate = m_aaNGate[ver];
 	double MGate = m_aaMGate[ver];
	double HGate = m_aaHGate[ver];
	double VM 	 = vrt_values[VMDisc<TDomain,TAlgebra>::_v_];

	// TODO Influx values needed
	// single channel type fluxes
	const number potassium_part_of_flux = m_g_K * pow(NGate,4) * (VM - m_rev_pot_K);
	const number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (VM + m_rev_pot_Na);
	const number leakage_part_of_flux = m_g_I * (VM + 54.4);

	number flux_value = (potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux);
	outCurrentValues.push_back(flux_value);
}



template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::Jacobi_sets(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outJFlux)
{
	double NGate = m_aaNGate[v];
	double MGate = m_aaMGate[v];
	double HGate = m_aaHGate[v];


	number Jac = (m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);

	outJFlux.push_back(Jac);

}


////////////////////////////////////////////////
// Methods for HH-Channel-Nernst-Class
////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::
set_accuracy(double ac)
{
	m_accuracy = ac;
}


template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::
set_conductivities(number Na, number K, number L)
{
	m_g_K = K;
	m_g_Na = Na;
	m_g_I = L;
}


// Methods for using gatings
template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::init(number time, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct)
{
	//UG_LOG("init starts");
	// attach attachments

	if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->m_MGate))
		UG_THROW("Attachment necessary (MGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->m_MGate);

	if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->m_HGate))
		UG_THROW("Attachment necessary (HGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->m_HGate);

	if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->m_NGate))
		UG_THROW("Attachment necessary (NGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->m_NGate);

	// create attachment accessors
	this->m_aaMGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->m_MGate);
	this->m_aaNGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->m_NGate);
	this->m_aaHGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->m_HGate);

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try {ssGrp = SubsetGroup(spGridFct->domain()->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");


	std::vector<DoFIndex> multInd;

	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = spGridFct->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = spGridFct->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get Vm
			const size_t VM_fct_ind = spGridFct->fct_id_by_name("VM");
			spGridFct->dof_indices(*iter, VM_fct_ind, multInd);
			UG_ASSERT(multInd.size() == 1, "multi-index has size != 1");
			number VM = DoFRef(*spGridFct, multInd[0]);

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
			this->m_aaHGate[*iter] = AlphaHh / (AlphaHh + BetaHh);
			this->m_aaMGate[*iter] = AlphaHm / (AlphaHm + BetaHm);
			this->m_aaNGate[*iter] = AlphaHn / (AlphaHn + BetaHn);
		}
	}
}

template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::update_gating(number newTime, ConstSmartPtr<typename TAlgebra::vector_type> uOld)
{
	//std::cout << "update gatting start" << std::endl;
	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(m_vmDisc->approx_space()->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	ConstSmartPtr<DoFDistribution> dd = m_vmDisc->approx_space()->dof_distribution(GridLevel::TOP);

	std::vector<DoFIndex> multInd;

	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get Vm
			const size_t VM_fct_ind = dd->fct_id_by_name("VM");
			dd->dof_indices(*iter, VM_fct_ind, multInd);
			UG_ASSERT(multInd.size() == 1, "multi-index has size != 1");
			number VM = DoFRef(*uOld, multInd[0]);

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

			// TODO: check that this is correct!
			number dt = newTime;

			number rate_h = ((AlphaHh/(AlphaHh+BetaHh)) - m_aaHGate[*iter]) / (1.0/(AlphaHh+BetaHh)) * dt;
			number rate_m = ((AlphaHm/(AlphaHm+BetaHm)) - m_aaMGate[*iter]) / (1.0/(AlphaHm+BetaHm)) * dt;
			number rate_n = ((AlphaHn/(AlphaHn+BetaHn)) - m_aaNGate[*iter]) / (1.0/(AlphaHn+BetaHn)) * dt;

			m_aaHGate[*iter] += rate_h;
			m_aaMGate[*iter] += rate_m;
			m_aaNGate[*iter] += rate_n;
		}
	}
}


template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments out of Vertex
	number NGate = m_aaNGate[ver];
	number MGate = m_aaMGate[ver];
	number HGate = m_aaHGate[ver];
	number v 	 = vrt_values[VMDisc<TDomain,TAlgebra>::_v_];
	number k 	 = vrt_values[VMDisc<TDomain,TAlgebra>::_k_];
	number na 	 = vrt_values[VMDisc<TDomain,TAlgebra>::_na_];


	UG_ASSERT(m_vmDisc, "Channel has not been assigned a vmDisc object yet!");
	const number helpV = 1e3*(m_R*m_T)/m_F;
	number potassium_nernst_eq 	= helpV*(std::log(m_vmDisc->k_out/k));
	number sodium_nernst_eq	 	= helpV*(std::log(m_vmDisc->na_out/na));

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



template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::Jacobi_sets(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outJFlux)
{
	number NGate = m_aaNGate[v];
	number MGate = m_aaMGate[v];
	number HGate = m_aaHGate[v];
	number k 	 = vrt_values[VMDisc<TDomain,TAlgebra>::_k_];
	number na 	 = vrt_values[VMDisc<TDomain,TAlgebra>::_na_];

	UG_ASSERT(m_vmDisc, "Channel has not been assigned a vmDisc object yet!");
	const number helpV = (m_R*m_T)/m_F;
	const number potassium_nernst_eq_dK 	=  helpV * (-m_vmDisc->k_out/k)*0.18; //helpV * (-K_out/pow(u(_K_,co),2));
	const number sodium_nernst_eq_dNa		=  -helpV * (-m_vmDisc->na_out/na)*0.003;

	outJFlux.push_back(m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);
	outJFlux.push_back(sodium_nernst_eq_dNa);
	outJFlux.push_back(potassium_nernst_eq_dK);

	//std::cout << "outJFlux: " << outJFlux[0] << ", " << outJFlux[1] << ", " << outJFlux[2] << ", " << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////


#ifdef UG_DIM_1
#ifdef UG_CPU_1
template class IChannel<Domain1d, CPUAlgebra>;
template class ChannelHH<Domain1d, CPUAlgebra>;
template class ChannelHHNernst<Domain1d, CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class IChannel<Domain1d, CPUBlockAlgebra<2> >;
template class ChannelHH<Domain1d, CPUBlockAlgebra<2> >;
template class ChannelHHNernst<Domain1d, CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class IChannel<Domain1d, CPUBlockAlgebra<3> >;
template class ChannelHH<Domain1d, CPUBlockAlgebra<3> >;
template class ChannelHHNernst<Domain1d, CPUBlockAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class IChannel<Domain1d, CPUBlockAlgebra<4> >;
template class ChannelHH<Domain1d, CPUBlockAlgebra<4> >;
template class ChannelHHNernst<Domain1d, CPUBlockAlgebra<4> >;
#endif
#ifdef UG_CPU_VAR
template class IChannel<Domain1d, CPUVariableBlockAlgebra >;
template class ChannelHH<Domain1d, CPUVariableBlockAlgebra >;
template class ChannelHHNernst<Domain1d, CPUVariableBlockAlgebra >;
#endif
#endif



#ifdef UG_DIM_2
#ifdef UG_CPU_1
template class IChannel<Domain2d, CPUAlgebra>;
template class ChannelHH<Domain2d, CPUAlgebra>;
template class ChannelHHNernst<Domain2d, CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class IChannel<Domain2d, CPUBlockAlgebra<2> >;
template class ChannelHH<Domain2d, CPUBlockAlgebra<2> >;
template class ChannelHHNernst<Domain2d, CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class IChannel<Domain2d, CPUBlockAlgebra<3> >;
template class ChannelHH<Domain2d, CPUBlockAlgebra<3> >;
template class ChannelHHNernst<Domain2d, CPUBlockAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class IChannel<Domain2d, CPUBlockAlgebra<4> >;
template class ChannelHH<Domain2d, CPUBlockAlgebra<4> >;
template class ChannelHHNernst<Domain2d, CPUBlockAlgebra<4> >;
#endif
#ifdef UG_CPU_VAR
template class IChannel<Domain2d, CPUVariableBlockAlgebra >;
template class ChannelHH<Domain2d, CPUVariableBlockAlgebra >;
template class ChannelHHNernst<Domain2d, CPUVariableBlockAlgebra >;
#endif
#endif


#ifdef UG_DIM_3
#ifdef UG_CPU_1
template class IChannel<Domain3d, CPUAlgebra>;
template class ChannelHH<Domain3d, CPUAlgebra>;
template class ChannelHHNernst<Domain3d, CPUAlgebra >;
#endif
#ifdef UG_CPU_2
template class IChannel<Domain3d, CPUBlockAlgebra<2> >;
template class ChannelHH<Domain3d, CPUBlockAlgebra<2> >;
template class ChannelHHNernst<Domain3d, CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class IChannel<Domain3d, CPUBlockAlgebra<3> >;
template class ChannelHH<Domain3d, CPUBlockAlgebra<3> >;
template class ChannelHHNernst<Domain3d, CPUBlockAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class IChannel<Domain3d, CPUBlockAlgebra<4> >;
template class ChannelHH<Domain3d, CPUBlockAlgebra<4> >;
template class ChannelHHNernst<Domain3d, CPUBlockAlgebra<4> >;
#endif
#ifdef UG_CPU_VAR
template class IChannel<Domain3d, CPUVariableBlockAlgebra >;
template class ChannelHH<Domain3d, CPUVariableBlockAlgebra >;
template class ChannelHHNernst<Domain3d, CPUVariableBlockAlgebra >;
#endif
#endif




} /* namespace ug */
