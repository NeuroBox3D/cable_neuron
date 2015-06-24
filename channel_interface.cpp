/*
 * channel_interface.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann
 */

#include "channel_interface.h"


namespace ug {
namespace cable {


////////////////////////////////////////////////
// Methods for HH-Channel-Class
////////////////////////////////////////////////


template<typename TDomain> void ChannelHH<TDomain>::set_log_nGate(bool bLogNGate) { m_log_nGate = bLogNGate; }
template<typename TDomain> void ChannelHH<TDomain>::set_log_hGate(bool bLogHGate) { m_log_hGate = bLogHGate; }
template<typename TDomain> void ChannelHH<TDomain>::set_log_mGate(bool bLogMGate) { m_log_mGate = bLogMGate; }


template<typename TDomain>
void ChannelHH<TDomain>::
set_conductivities(number Na, number K, number L)
{
	m_g_K = K;
	m_g_Na = Na;
	m_g_I = L;
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
std::vector<number> ChannelHH<TDomain>::allGatingAccesors(number x, number y, number z)
{
	//var for output
	std::vector<number> GatingAccesors;
	typedef ug::MathVector<TDomain::dim> position_type;

	position_type coord;

	if (coord.size()==1)
	{
		coord[0] = x;
	}
	if (coord.size()==2)
	{
		coord[0] = x;
		coord[1] = y;
	}
	if (coord.size()==3)
	{
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}

	// accessors
	typedef Attachment<position_type> position_attachment_type;
	typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

	// Definitions for Iterating over all Elements
	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try{ ssGrp = SubsetGroup(m_pVMDisc->approx_space()->domain()->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	//UG_LOG("Channel: Before iteration" << std::endl);

	itType iter;
	number bestDistSq, distSq;
	Vertex* bestVrt;


	// Iterate only if there is one Gatting needed
	if (m_log_mGate==true || m_log_hGate==true || m_log_nGate==true)
	{
		// iterating over all elements
		for (size_t si=0; si < ssGrp.size(); si++)
		{
			itType iterBegin = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);
			itType iterEnd = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]);

			const position_accessor_type& aaPos = m_pVMDisc->approx_space()->domain()->position_accessor();
			// if the right vertex of needed Position is found write out values
			if (si==0)
			{
				bestVrt = *iterBegin;
				bestDistSq = VecDistanceSq(coord, aaPos[bestVrt]);
			}
				iter = iterBegin;
				iter++;
				while(iter != iterEnd)
				{
					distSq = VecDistanceSq(coord, aaPos[*iter]);
					if(distSq < bestDistSq)
					{
						bestDistSq = distSq;
						bestVrt = *iter;
					}
					++iter;
				}
		}

		if (m_log_mGate == true)
			GatingAccesors.push_back(this->m_aaMGate[bestVrt]);
		if (m_log_hGate == true)
			GatingAccesors.push_back(this->m_aaHGate[bestVrt]);
		if (m_log_nGate == true)
			GatingAccesors.push_back(this->m_aaNGate[bestVrt]);
	}


	return GatingAccesors;
}


template<typename TDomain>
void ChannelHH<TDomain>::init_attachments()
{
	// attach attachments
	SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid();

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
	m_aaMGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_MGate);
	m_aaNGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_NGate);
	m_aaHGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_HGate);
}


// Methods for using gatings
template<typename TDomain>
void ChannelHH<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	number VM = vrt_values[VMDisc<TDomain>::_v_];

	// values for m gate
	number AlphaHm = 0.1 * vtrap(-(VM+40.0),10.0);
	number BetaHm =  4.0 * exp(-(VM+65.0)/18.0);

	// values for n gate
	number AlphaHn = 0.01*vtrap(-(VM+55.0),10.0);
	number BetaHn = 0.125*exp(-(VM+65.0)/80.0);

	// values for h gate
	number AlphaHh = 0.07 * exp(-(VM+65.0)/20.0);
	number BetaHh = 1.0 / (exp(-(VM+35.0)/10.0) + 1.0);

	// setting initial gating params as equilibrium states
	this->m_aaHGate[vrt] = AlphaHh / (AlphaHh + BetaHh);
	this->m_aaMGate[vrt] = AlphaHm / (AlphaHm + BetaHm);
	this->m_aaNGate[vrt] = AlphaHn / (AlphaHn + BetaHn);
}

template<typename TDomain>
void ChannelHH<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{
	number dt = newTime - m_pVMDisc->time();
	number VM = vrt_values[VMDisc<TDomain>::_v_];

	// values for m gate
	number AlphaHm = 0.1 * vtrap(-(VM+40.0),10.0);
	number BetaHm =  4.0 * exp(-(VM+65.0)/18.0);
	number sumABHm = 4.5 * (AlphaHm+BetaHm);

	// values for n gate
	number AlphaHn = 0.01 * vtrap(-(VM+55.0),10.0);
	number BetaHn = 0.125 * exp(-(VM+65.0)/80.0);
	number sumABHn = 4.5 * (AlphaHn+BetaHn);

	// values for h gate
	number AlphaHh = 0.07 * exp(-(VM+65.0)/20.0);
	number BetaHh = 1.0 / (exp(-(VM+35.0)/10.0) + 1.0);
	number sumABHh = 4.5 * (AlphaHh+BetaHh);

	number rate_h = ((AlphaHh/((AlphaHh+BetaHh))) - m_aaHGate[vrt]) * sumABHh * dt;
	number rate_m = ((AlphaHm/((AlphaHm+BetaHm))) - m_aaMGate[vrt]) * sumABHm * dt;
	number rate_n = ((AlphaHn/((AlphaHn+BetaHn))) - m_aaNGate[vrt]) * sumABHn * dt;

	m_aaHGate[vrt] += rate_h;
	m_aaMGate[vrt] += rate_m;
	m_aaNGate[vrt] += rate_n;
	//std::cout << "VM: " << VM << "   h: "<< m_aaHGate[vrt] << "   m: "<< m_aaMGate[vrt] <<  "   n: "<< m_aaNGate[vrt] <<std::endl;
	//std::cout << "Rates: " << rate_h << " , " << rate_m << " , " << rate_n << std::endl;
}


template<typename TDomain>
void ChannelHH<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	number NGate = m_aaNGate[vrt];
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];
	number VM 	 = vrt_values[VMDisc<TDomain>::_v_];

	number rev_pot_K = m_pVMDisc->ek();
	number rev_pot_Na = m_pVMDisc->ena();
	number rev_pot_leak = m_pVMDisc->eleak();

	// single channel type fluxes
	const number potassium_part_of_flux = m_g_K * pow(NGate,4) * (VM - rev_pot_K);
	const number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (VM - rev_pot_Na);
	const number leakage_part_of_flux = m_g_I * (VM - rev_pot_leak);

	/*
	std::cout << "VM: " << VM << std::endl;
	std::cout << "NGate: " << NGate << std::endl;
	std::cout << "MGate: " << MGate << std::endl;
	std::cout << "HGate: " << HGate << std::endl;
	std::cout << "potassium_part_of_flux: " << potassium_part_of_flux << std::endl;
	std::cout << "sodium_part_of_flux: " << sodium_part_of_flux << std::endl;
	std::cout << "leakage_part_of_flux: " << leakage_part_of_flux << std::endl;
	std::cout << "flux: " << potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux << std::endl;
	*/

	number flux_value = (potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux);
	outCurrentValues.push_back(flux_value);
}


#if 0
template<typename TDomain>
void ChannelHH<TDomain>::Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux)
{
	number NGate = m_aaNGate[vrt];
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];


	number Jac = (m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);

	outJFlux.push_back(Jac);

}
#endif

////////////////////////////////////////////////
// Methods for HH-Channel-Nernst-Class
////////////////////////////////////////////////


template<typename TDomain> void ChannelHHNernst<TDomain>::set_log_nGate(bool bLogNGate) { m_log_nGate = bLogNGate; }
template<typename TDomain> void ChannelHHNernst<TDomain>::set_log_hGate(bool bLogHGate) { m_log_hGate = bLogHGate; }
template<typename TDomain> void ChannelHHNernst<TDomain>::set_log_mGate(bool bLogMGate) { m_log_mGate = bLogMGate; }



template<typename TDomain>
void ChannelHHNernst<TDomain>::
set_conductivities(number Na, number K, number L)
{
	m_g_K = K;
	m_g_Na = Na;
	m_g_I = L;
}


template<typename TDomain>
number ChannelHHNernst<TDomain>::
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
void ChannelHHNernst<TDomain>::vm_disc_available()
{
	init_attachments();
}

template<typename TDomain>
std::vector<number> ChannelHHNernst<TDomain>::allGatingAccesors(number x, number y, number z)
{
	//var for output
	std::vector<number> GatingAccesors;
	typedef ug::MathVector<TDomain::dim> position_type;

	position_type coord;
	if (coord.size()==1)
	{
		coord[0] = x;
	}
	if (coord.size()==2)
	{
		coord[0] = x;
		coord[1] = y;
	}
	if (coord.size()==3)
	{
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}



	// accessors
	typedef Attachment<position_type> position_attachment_type;
	typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

	// Definitions for Iterating over all Elements
	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try{ ssGrp = SubsetGroup(m_pVMDisc->approx_space()->domain()->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	//UG_LOG("Channel: Before iteration" << std::endl);

	itType iter;
	number bestDistSq, distSq;
	Vertex* bestVrt;


	// Iterate only if there is one Gatting needed
	if (m_log_mGate==true || m_log_hGate==true || m_log_nGate==true)
	{
		// iterating over all elements
		for (size_t si=0; si < ssGrp.size(); si++)
		{
			itType iterBegin = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);
			itType iterEnd = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]);

			const position_accessor_type& aaPos = m_pVMDisc->approx_space()->domain()->position_accessor();

			if (si==0)
			{
				bestVrt = *iterBegin;
				bestDistSq = VecDistanceSq(coord, aaPos[bestVrt]);
			}
			iter = iterBegin;
			iter++;
			while(iter != iterEnd)
			{
				distSq = VecDistanceSq(coord, aaPos[*iter]);
				if(distSq < bestDistSq)
				{
					bestDistSq = distSq;
					bestVrt = *iter;
				}
				++iter;
			}
		}

		if (m_log_mGate == true)
			GatingAccesors.push_back(this->m_aaMGate[bestVrt]);
		if (m_log_hGate == true)
			GatingAccesors.push_back(this->m_aaHGate[bestVrt]);
		if (m_log_nGate == true)
			GatingAccesors.push_back(this->m_aaNGate[bestVrt]);
	}


	return GatingAccesors;
}

template<typename TDomain>
void ChannelHHNernst<TDomain>::init_attachments()
{
	// attach attachments
	SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid();

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
	m_aaMGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_MGate);
	m_aaNGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_NGate);
	m_aaHGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_HGate);
}



// Methods for using gatings
template<typename TDomain>
void ChannelHHNernst<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	number VM = vrt_values[VMDisc<TDomain>::_v_];

	// values for m gate
	number AlphaHm = 0.1 * vtrap(-(VM+40.0),10.0);
	number BetaHm =  4.0 * exp(-(VM+65.0)/18.0);

	// values for n gate
	number AlphaHn = 0.01*vtrap(-(VM+55.0),10.0);
	number BetaHn = 0.125*exp(-(VM+65.0)/80.0);

	// values for h gate
	number AlphaHh = 0.07 * exp(-(VM+65.0)/20.0);
	number BetaHh = 1.0 / (exp(-(VM+35.0)/10.0) + 1.0);

	// setting initial gating params as equilibrium states
	this->m_aaHGate[vrt] = AlphaHh / (AlphaHh + BetaHh);
	this->m_aaMGate[vrt] = AlphaHm / (AlphaHm + BetaHm);
	this->m_aaNGate[vrt] = AlphaHn / (AlphaHn + BetaHn);
}

template<typename TDomain>
void ChannelHHNernst<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{
	number dt = newTime - m_pVMDisc->time();
	number VM = vrt_values[VMDisc<TDomain>::_v_];

	// set new gating states
	// values for m gate
	number AlphaHm = 0.1 * vtrap(-(VM+40.0),10.0);
	number BetaHm =  4.0 * exp(-(VM+65.0)/18.0);

	// values for n gate
	number AlphaHn = 0.01 * vtrap(-(VM+55.0),10.0);
	number BetaHn = 0.125 * exp(-(VM+65.0)/80.0);

	// values for h gate
	number AlphaHh = 0.07 * exp(-(VM+65.0)/20.0);
	number BetaHh = 1.0 / (exp(-(VM+35.0)/10.0) + 1.0);

	number rate_h = ((AlphaHh/(AlphaHh+BetaHh)) - m_aaHGate[vrt]) / (1.0/(AlphaHh+BetaHh)) * dt;
	number rate_m = ((AlphaHm/(AlphaHm+BetaHm)) - m_aaMGate[vrt]) / (1.0/(AlphaHm+BetaHm)) * dt;
	number rate_n = ((AlphaHn/(AlphaHn+BetaHn)) - m_aaNGate[vrt]) / (1.0/(AlphaHn+BetaHn)) * dt;

	m_aaHGate[vrt] += rate_h;
	m_aaMGate[vrt] += rate_m;
	m_aaNGate[vrt] += rate_n;
}


template<typename TDomain>
void ChannelHHNernst<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments out of Vertex
	number NGate = m_aaNGate[vrt];
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];
	number v 	 = vrt_values[m_pVMDisc->_v_];
	number k 	 = vrt_values[m_pVMDisc->_k_];
	number na 	 = vrt_values[m_pVMDisc->_na_];

	//UG_ASSERT(m_pVMDisc->valid(), "Channel has not been assigned a vmDisc object yet!");
	const number helpV = 1e3*(m_R*m_T)/m_F;
	number potassium_nernst_eq 	= helpV*(std::log(m_pVMDisc->k_out()/k));
	number sodium_nernst_eq	 	= helpV*(std::log(m_pVMDisc->na_out()/na));
	number rev_pot_leak = m_pVMDisc->eleak();

	// single channel ion fluxes
	number potassium_part_of_flux = m_g_K * pow(NGate,4) * (v - potassium_nernst_eq);
	number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (v - sodium_nernst_eq);
	number leakage_part_of_flux = m_g_I * (v - rev_pot_leak);

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

	UG_ASSERT(m_pVMDisc, "Channel has not been assigned a vmDisc object yet!");
	const number helpV = (m_R*m_T)/m_F;
	const number potassium_nernst_eq_dK 	=  helpV * (-m_pVMDisc->k_out()/k)*0.18; //helpV * (-K_out/pow(u(_K_,co),2));
	const number sodium_nernst_eq_dNa		=  -helpV * (-m_pVMDisc->na_out()/na)*0.003;

	outJFlux.push_back(m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);
	outJFlux.push_back(sodium_nernst_eq_dNa);
	outJFlux.push_back(potassium_nernst_eq_dK);

	//std::cout << "outJFlux: " << outJFlux[0] << ", " << outJFlux[1] << ", " << outJFlux[2] << ", " << std::endl;
}
#endif







////////////////////////////////////////////////
// Methods for Leakeage-Channel-Class
////////////////////////////////////////////////
template<typename TDomain>
void ChannelLeak<TDomain>::
set_leak_cond(number L)
{
	m_g_I = L;
}


template<typename TDomain>
void ChannelLeak<TDomain>::vm_disc_available()
{

}


template<typename TDomain>
void ChannelLeak<TDomain>::init_attachments()
{

}

template<typename TDomain>
std::vector<number> ChannelLeak<TDomain>::allGatingAccesors(number x, number y, number z)
{
	std::vector<number> GatingAccesors;


	return GatingAccesors;
}

// Methods for using gatings
template<typename TDomain>
void ChannelLeak<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	// nothing to do
}

template<typename TDomain>
void ChannelLeak<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{
	// nothing to do
}


template<typename TDomain>
void ChannelLeak<TDomain>::ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	number VM 	 = vrt_values[m_pVMDisc->_v_];
	number leak_equilibrium = m_pVMDisc->eleak();

	const number leakage_part_of_flux = m_g_I * (VM - leak_equilibrium);

	number flux_value = (leakage_part_of_flux);
	outCurrentValues.push_back(flux_value);
}


#if 0
template<typename TDomain>
void ChannelLeak<TDomain>::Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux)
{
	number NGate = m_aaNGate[vrt];
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];


	number Jac = (m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);

	outJFlux.push_back(Jac);

}

#endif




////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class IChannel<Domain1d>;
	template class ChannelHH<Domain1d>;
	template class ChannelLeak<Domain1d>;
	template class ChannelHHNernst<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class IChannel<Domain2d>;
	template class ChannelHH<Domain2d>;
	template class ChannelLeak<Domain2d>;
	template class ChannelHHNernst<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class IChannel<Domain3d>;
	template class ChannelHH<Domain3d>;
	template class ChannelLeak<Domain3d>;
	template class ChannelHHNernst<Domain3d>;
#endif

} // namespace cable
} // namespace ug
