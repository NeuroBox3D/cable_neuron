/*
 * channel_hh.cpp
 *
 *  Created on: 29.10.2014
 *      Author: pgottmann, mbreit
 */

#include "channel_hh.h"
#include <limits> // numeric_limits

namespace ug {
namespace cable {


////////////////////////////////////////////////
// Methods for HH-Channel-Class
////////////////////////////////////////////////

template<typename TDomain>
ChannelHH<TDomain>::ChannelHH(const char* functions, const char* subsets)
try : IChannel<TDomain>(functions, subsets),
m_g_K(3.6e-4), m_g_Na(1.2e-3),
m_log_nGate(false), m_log_hGate(false), m_log_mGate(false) {}
UG_CATCH_THROW("Error in ChannelHH initializer list.");

template<typename TDomain>
ChannelHH<TDomain>::ChannelHH
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : IChannel<TDomain>(functions, subsets),
m_g_K(3.6e-4), m_g_Na(1.2e-3),
m_log_nGate(false), m_log_hGate(false), m_log_mGate(false) {}
UG_CATCH_THROW("Error in ChannelHH initializer list.");


template<typename TDomain>
std::string ChannelHH<TDomain>::
name()
{
	return std::string("ChannelHH");
}


template<typename TDomain> void ChannelHH<TDomain>::set_log_nGate(bool bLogNGate) { m_log_nGate = bLogNGate; }
template<typename TDomain> void ChannelHH<TDomain>::set_log_hGate(bool bLogHGate) { m_log_hGate = bLogHGate; }
template<typename TDomain> void ChannelHH<TDomain>::set_log_mGate(bool bLogMGate) { m_log_mGate = bLogMGate; }


template<typename TDomain>
void ChannelHH<TDomain>::
set_conductances(number gK, number gNa)
{
	m_g_K = gK;
	m_g_Na = gNa;
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
std::vector<number> ChannelHH<TDomain>::state_values(number x, number y, number z)
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
			itType iterBegin = m_pVMDisc->approx_space()->dof_distribution(GridLevel(), false)->template begin<Vertex>(ssGrp[si]);
			itType iterEnd = m_pVMDisc->approx_space()->dof_distribution(GridLevel(), false)->template end<Vertex>(ssGrp[si]);

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
	number tmp = m_pVMDisc->temperature_celsius();
	number tmp_factor = std::pow(2.3, (tmp-23.0)/10.0);

	// values for m gate
	number AlphaHm = 0.1 * vtrap(-(VM+40.0),10.0);
	number BetaHm =  4.0 * exp(-(VM+65.0)/18.0);

	// values for n gate
	number AlphaHn = 0.01 * vtrap(-(VM+55.0),10.0);
	number BetaHn = 0.125 * exp(-(VM+65.0)/80.0);

	// values for h gate
	number AlphaHh = 0.07 * exp(-(VM+65.0)/20.0);
	number BetaHh = 1.0 / (exp(-(VM+35.0)/10.0) + 1.0);

	number rate_h = tmp_factor * (AlphaHh - m_aaHGate[vrt] * (AlphaHh+BetaHh));
	number rate_m = tmp_factor * (AlphaHm - m_aaMGate[vrt] * (AlphaHm+BetaHm));
	number rate_n = tmp_factor * (AlphaHn - m_aaNGate[vrt] * (AlphaHn+BetaHn));

	m_aaHGate[vrt] += rate_h * dt;
	m_aaMGate[vrt] += rate_m * dt;
	m_aaNGate[vrt] += rate_n * dt;
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

	number tmp = m_pVMDisc->temperature_celsius();
	number tmp_factor = std::pow(2.3, (tmp-23.0)/10.0);

	// single channel type fluxes
	const number potassium_part_of_flux = tmp_factor * m_g_K * pow(NGate,4) * (VM - rev_pot_K);
	const number sodium_part_of_flux =    tmp_factor * m_g_Na * pow(MGate,3) * HGate * (VM - rev_pot_Na);

	/*
	std::cout << "VM: " << VM << std::endl;
	std::cout << "NGate: " << NGate << std::endl;
	std::cout << "MGate: " << MGate << std::endl;
	std::cout << "HGate: " << HGate << std::endl;
	std::cout << "potassium_part_of_flux: " << potassium_part_of_flux << std::endl;
	std::cout << "sodium_part_of_flux: " << sodium_part_of_flux << std::endl;
	std::cout << "flux: " << potassium_part_of_flux + sodium_part_of_flux << std::endl;
	*/

	number flux_value = (potassium_part_of_flux + sodium_part_of_flux);
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




template<typename TDomain>
number ChannelHH<TDomain>::
lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values)
{
	// getting attachments for vertex
	number NGate = m_aaNGate[vrt];
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];

	number tmp = m_pVMDisc->temperature_celsius();
	number tmp_factor = std::pow(2.3, (tmp-23.0)/10.0);

	// single channel type fluxes
	const number potassium_part_of_flux = tmp_factor * m_g_K * pow(NGate,4);
	const number sodium_part_of_flux =    tmp_factor * m_g_Na * pow(MGate,3) * HGate;

	return potassium_part_of_flux + sodium_part_of_flux;
}


template<typename TDomain>
void ChannelHH<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing VMDisc fct indices which this channel writes to
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_);
}

////////////////////////////////////////////////
// Methods for HH-Channel-Nernst-Class
////////////////////////////////////////////////

template <typename TDomain>
ChannelHHNernst<TDomain>::ChannelHHNernst(const char* functions, const char* subsets)
try : IChannel<TDomain>(functions, subsets),
m_g_K(3.6e-4), m_g_Na(1.2e-3),
m_log_nGate(false), m_log_hGate(false), m_log_mGate(false) {}
UG_CATCH_THROW("Error in ChannelHHNernst initializer list.");

template <typename TDomain>
ChannelHHNernst<TDomain>::ChannelHHNernst
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : IChannel<TDomain>(functions, subsets),
m_g_K(3.6e-4), m_g_Na(1.2e-3),
m_log_nGate(false), m_log_hGate(false), m_log_mGate(false) {}
UG_CATCH_THROW("Error in ChannelHHNernst initializer list.");


template<typename TDomain>
std::string ChannelHHNernst<TDomain>::
name()
{
	return std::string("ChannelHHNernst");
}


template<typename TDomain> void ChannelHHNernst<TDomain>::set_log_nGate(bool bLogNGate) { m_log_nGate = bLogNGate; }
template<typename TDomain> void ChannelHHNernst<TDomain>::set_log_hGate(bool bLogHGate) { m_log_hGate = bLogHGate; }
template<typename TDomain> void ChannelHHNernst<TDomain>::set_log_mGate(bool bLogMGate) { m_log_mGate = bLogMGate; }


template<typename TDomain>
void ChannelHHNernst<TDomain>::
set_conductances(number gK, number gNa)
{
	m_g_K = gK;
	m_g_Na = gNa;
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
void ChannelHHNernst<TDomain>::vmDisc_available()
{
	init_attachments();
}

template<typename TDomain>
std::vector<number> ChannelHHNernst<TDomain>::state_values(number x, number y, number z)
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
	number bestDistSq = std::numeric_limits<number>::max();
	number distSq;
	Vertex* bestVrt = NULL;


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

		if (!bestVrt) UG_THROW("No vertex found for coords.");

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

	number rate_h = AlphaHh - m_aaHGate[vrt]*(AlphaHh+BetaHh);
	number rate_m = AlphaHm - m_aaMGate[vrt]*(AlphaHm+BetaHm);
	number rate_n = AlphaHn - m_aaNGate[vrt]*(AlphaHn+BetaHn);

	m_aaHGate[vrt] += rate_h * dt;
	m_aaMGate[vrt] += rate_m * dt;
	m_aaNGate[vrt] += rate_n * dt;
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

	const number R = m_pVMDisc->R;
	const number F = m_pVMDisc->F;
	const number T = m_pVMDisc->temperature();

	//UG_ASSERT(m_pVMDisc->valid(), "Channel has not been assigned a vmDisc object yet!");
	const number helpV = 1e3 * R*T/F;
	number potassium_nernst_eq 	= helpV*(std::log(m_pVMDisc->k_out()/k));
	number sodium_nernst_eq	 	= helpV*(std::log(m_pVMDisc->na_out()/na));

	// single channel ion fluxes
	number potassium_part_of_flux = m_g_K * pow(NGate,4) * (v - potassium_nernst_eq);
	number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (v - sodium_nernst_eq);

	outCurrentValues.push_back(potassium_part_of_flux + sodium_part_of_flux);
	outCurrentValues.push_back(potassium_part_of_flux / F);
	outCurrentValues.push_back(sodium_part_of_flux / F);

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
	outJFlux.push_back(potassium_nernst_eq_dK);
	outJFlux.push_back(sodium_nernst_eq_dNa);

	//std::cout << "outJFlux: " << outJFlux[0] << ", " << outJFlux[1] << ", " << outJFlux[2] << ", " << std::endl;
}
#endif


template<typename TDomain>
void ChannelHHNernst<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing VMDisc fct indices which this channel writes to
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_);
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_k_);
	this->m_vWFctInd.push_back(VMDisc<TDomain>::_na_);
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class ChannelHH<Domain1d>;
	template class ChannelHHNernst<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class ChannelHH<Domain2d>;
	template class ChannelHHNernst<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class ChannelHH<Domain3d>;
	template class ChannelHHNernst<Domain3d>;
#endif

} // namespace cable
} // namespace ug
