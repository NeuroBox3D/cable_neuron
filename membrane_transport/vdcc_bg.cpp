/*
 * vdcc_bg.cpp
 *
 *  Created on: 24.08.2015
 *      Author: pgottmann, mbreit
 */

#include "vdcc_bg.h"
#include <limits> // numeric_limits

namespace ug {
namespace cable_neuron {


////////////////////////////////////////////////
// Methods for HH-Channel-Class
////////////////////////////////////////////////

template<typename TDomain>
VDCC_BG_cable<TDomain>::VDCC_BG_cable(const char* functions, const char* subsets)
try : ICableMembraneTransport<TDomain>(functions, subsets),
	m_mp(2), m_z_m(3.4), m_v12_m(-21.0), m_tau0_m(1.5),
	m_hp(1), m_z_h(-2.0), m_v12_h(-40.0), m_tau0_h(75.0),
	m_perm(3.8e-7),
	m_log_hGate(false), m_log_mGate(false) {}
UG_CATCH_THROW("Error in VDCC_BG_cable initializer list.");

template<typename TDomain>
VDCC_BG_cable<TDomain>::VDCC_BG_cable
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
try : ICableMembraneTransport<TDomain>(functions, subsets),
	m_mp(2), m_z_m(3.4), m_v12_m(-21.0), m_tau0_m(1.5),
	m_hp(1), m_z_h(-2.0), m_v12_h(-40.0), m_tau0_h(75.0),
	m_perm(3.8e-7),
	m_log_hGate(false), m_log_mGate(false) {}
UG_CATCH_THROW("Error in VDCC_BG_cable initializer list.");


template<typename TDomain>
std::string VDCC_BG_cable<TDomain>::
name()
{
	return std::string("VDCC_BG_cable");
}


template<typename TDomain> void VDCC_BG_cable<TDomain>::set_log_hGate(bool bLogHGate) { m_log_hGate = bLogHGate; }
template<typename TDomain> void VDCC_BG_cable<TDomain>::set_log_mGate(bool bLogMGate) { m_log_mGate = bLogMGate; }



template<typename TDomain>
void VDCC_BG_cable<TDomain>::ce_obj_available()
{
	init_attachments();
}

template<typename TDomain>
std::vector<number> VDCC_BG_cable<TDomain>::state_values(number x, number y, number z)
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
	try{ ssGrp = SubsetGroup(m_pCE->approx_space()->domain()->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	//UG_LOG("Channel: Before iteration" << std::endl);

	itType iter;
	number bestDistSq = std::numeric_limits<number>::max();
	number distSq;
	Vertex* bestVrt = NULL;


	// iterate only if any state is needed
	if (m_log_mGate || m_log_hGate)
	{
		// iterating over all elements
		for (size_t si=0; si < ssGrp.size(); si++)
		{
			itType iterBegin = m_pCE->approx_space()->dof_distribution(GridLevel(), false)->template begin<Vertex>(ssGrp[si]);
			itType iterEnd = m_pCE->approx_space()->dof_distribution(GridLevel(), false)->template end<Vertex>(ssGrp[si]);

			const position_accessor_type& aaPos = m_pCE->approx_space()->domain()->position_accessor();
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

		UG_COND_THROW(!bestVrt, "No vertex found.");

		if (m_log_mGate)
			GatingAccesors.push_back(this->m_aaMGate[bestVrt]);
		if (m_log_hGate)
			GatingAccesors.push_back(this->m_aaHGate[bestVrt]);
	}


	return GatingAccesors;
}


template<typename TDomain>
void VDCC_BG_cable<TDomain>::init_attachments()
{
	// attach attachments
	SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid();

	if (spGrid->has_vertex_attachment(m_MGate))
		UG_THROW("Attachment necessary (MGate) for Borg-Graham type VDCC dynamics "
				 "could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_MGate);

	if (spGrid->has_vertex_attachment(m_HGate))
		UG_THROW("Attachment necessary (HGate) for Borg-Graham type VDCC dynamics "
				 "could not be made, since it already exists.");
	spGrid->attach_to_vertices(m_HGate);

	// create attachment accessors
	m_aaMGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_MGate);
	m_aaHGate = Grid::AttachmentAccessor<Vertex, ANumber>(*spGrid, m_HGate);
}


// Methods for using gatings
template<typename TDomain>
void VDCC_BG_cable<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values)
{
	number VM = 1e3*vrt_values[CableEquation<TDomain>::_v_];	// potential in units of mV

	const number& R = m_pCE->R;
	const number& F = m_pCE->F;
	const number& T = m_pCE->temperature();

	m_aaMGate[vrt] =  1.0 / (1.0 + exp(-m_z_m * (VM - m_v12_m) * 1e-3*F/(R*T)));
	m_aaHGate[vrt] = 1.0 / (1.0 + exp(-m_z_h * (VM - m_v12_h) * 1e-3*F/(R*T)));
}

template<typename TDomain>
void VDCC_BG_cable<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values)
{
	number dt = 1e3*(newTime - m_pCE->time());					// time in units of ms
	number VM = 1e3*vrt_values[CableEquation<TDomain>::_v_];	// potential in units of mV

	const number& R = m_pCE->R;
	const number& F = m_pCE->F;
	const number& T = m_pCE->temperature();

	number& mGate = m_aaMGate[vrt];
	number& hGate = m_aaHGate[vrt];

	number m_inf = 1.0 / (1.0 + exp(-m_z_m * (VM - m_v12_m) * 1e-3*F/(R*T)));
	number h_inf = 1.0 / (1.0 + exp(-m_z_h * (VM - m_v12_h) * 1e-3*F/(R*T)));

	// forward step: implicit
	if (dt >= 0)
	{
		// For calculating the next gating step it is recommended to use a time step
		// size no larger than 1e-5s = 1e-2ms in order to meet a sufficient accuracy
		if (dt > 1e-2)
		{
			number vdcc_dt = 1e-2;
			number t0 = 0.0;

			// loop intermediate time steps until the current intermediate time point
			// is the final time point dt
			while (t0 < dt)
			{
				// compute next time point
				number t = t0 + vdcc_dt;

				// check if out of bounds, if yes:
				// set to final time point and adjust step size accordingly
				if (t > dt)
				{
					t = dt;
					vdcc_dt = dt - t0;
				}

				// compute next gating value
				mGate = (mGate + vdcc_dt/m_tau0_m * m_inf) / (1.0 + vdcc_dt/m_tau0_m);
				hGate = (hGate + vdcc_dt/m_tau0_h * h_inf) / (1.0 + vdcc_dt/m_tau0_h);

				// save new time as current time
				t0 = t;
			}
		}
		// sufficiently small time step size already specified
		else
		{
			mGate = (mGate + dt/m_tau0_m * m_inf) / (1.0 + dt/m_tau0_m);
			hGate = (hGate + dt/m_tau0_h * h_inf) / (1.0 + dt/m_tau0_h);
		}
	}

	// backward step: explicit // TODO: do the above backwards!
	else
	{
		if (dt < -1e-2) UG_THROW("time step too large; not implemented yet!");

		mGate += dt/m_tau0_m * (m_inf - mGate);
		hGate += dt/m_tau0_h * (h_inf - hGate);
	}
}


template<typename TDomain>
void VDCC_BG_cable<TDomain>::current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues)
{
	// getting attachments for vertex
	number MGate = m_aaMGate[vrt];
	number HGate = m_aaHGate[vrt];
	number VM 	 = vrt_values[CableEquation<TDomain>::_v_];
	number caCyt = vrt_values[CableEquation<TDomain>::_ca_];
	number caExt = m_pCE->ca_out();

	const number& R = m_pCE->R;
	const number& F = m_pCE->F;
	const number& T = m_pCE->temperature();

	number gating = pow(MGate, m_mp);
	gating *= pow(HGate, m_hp);

	// flux derived from Goldman-Hodgkin-Katz equation,
	number maxFluxDensity;

	// near V_m == 0: approximate by first order Taylor to avoid relative errors and div-by-0
	if (fabs(VM) < 1e-8) maxFluxDensity = m_perm * ((caExt - caCyt) - F/(R*T) * (caExt + caCyt)*VM);
	else maxFluxDensity = -m_perm * 2*F/(R*T) * VM * (caExt - caCyt*exp(2*F/(R*T)*VM)) / (1.0 - exp(2*F/(R*T)*VM));

	outCurrentValues.push_back(-gating * maxFluxDensity);
}


template<typename TDomain>
void VDCC_BG_cable<TDomain>::
specify_write_function_indices()
{
	// prepare vector containing CableEquation fct indices which this channel writes to
	this->m_vWFctInd.push_back(CableEquation<TDomain>::_ca_);
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
	template class VDCC_BG_cable<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class VDCC_BG_cable<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class VDCC_BG_cable<Domain3d>;
#endif

} // namespace cable_neuron
} // namespace ug
