/*
 * channel_interface.cpp
 *
 *  Created on: 29.10.2014
 *      Author: mbreit
 */

#include "channel_interface.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"


namespace ug {

// Methods for Interface class



template<typename TDomain>
template<typename TElem, typename TFVGeom>
void IChannel<TDomain>::prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("Cannot update Finite Volume Geometry.\n");

}

template<typename TDomain>
void IChannel<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_funcs(m_bNonRegularGrid);
}


/* IChannel:: use_hanging()
{
	return false;
}*/


// Methods for HH-Channel-Class


template<typename TDomain>
void ChannelHH<TDomain>::init(number time)
{
	// attach attachments
	if (this->m_mg->has_vertex_attachment(this->m_MGate))
		UG_THROW("Attachment necessary (MGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_MGate);

	if (this->m_mg->has_vertex_attachment(this->m_HGate))
	{
		if (this->m_mg->has_vertex_attachment(this->m_HGate))
			UG_THROW("Attachment necessary (HGate) for Hodgkin and Huxley channel dynamics "
					 "could not be made, since it already exists.");
		this->m_mg->attach_to_vertices(this->m_HGate);
	}

	if (this->m_mg->has_vertex_attachment(this->m_NGate))
	{
		if (this->m_mg->has_vertex_attachment(this->m_NGate))
			UG_THROW("Attachment necessary (NGate) for Hodgkin and Huxley channel dynamics "
					 "could not be made, since it already exists.");
		this->m_mg->attach_to_vertices(this->m_NGate);
	}


	if (this->m_mg->has_vertex_attachment(this->m_Vm))
		UG_THROW("Attachment necessary (Vm) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_Vm);


	// create attachment accessors
	this->m_aaMGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_MGate);
	this->m_aaNGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_NGate);
	this->m_aaHGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_HGate);
	this->m_aaVm = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_Vm);

	//TODO set accuracy later on another point
	number m_accuracy = 1e-6;


	// creates Multi index
	std::vector<DoFIndex> multInd;

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	const typename TDomain::position_accessor_type& aaPos = this->m_dom->position_accessor();
	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// TODO getting potential out of approxSpace into attachement
			// TODO later the number of attachements needs to be variable, also some outer concentrations have to be given
			const size_t fct = 0;
			const LFEID& lfeID = m_dd->local_finite_element_id(fct);
			//const
			//m_dd->constrained_vertex_dof_indices(fct, multInd, lfeID);
			m_dd->inner_dof_indices(*iter, fct, multInd);
			// we need here some better settings
			m_aaVm[*iter] = multInd[0];

			// Writting Gatting-Params in attachement
			// gating param h
			number AlphaHh = 0.07*exp(-(m_aaVm[*iter] + 65.0)/20.0);
			number BetaHh = 1.0/(exp(3.0-0.1*(m_aaVm[*iter]  + 65.0))+1.0);

			// gating param m
			number AlphaHm;
			number AlphaHm_test = exp(2.5-0.1*(m_aaVm[*iter] + 65.0))-1.0;
			if (fabs(AlphaHm_test) > m_accuracy)
				AlphaHm = (2.5 - 0.1*( m_aaVm[*iter] + 65.0)) / AlphaHm_test;
			else
				AlphaHm = 1.0;

			number BetaHm = 4.0*exp((-m_aaVm[*iter] + 65.0)/18.0);

			// gating param n
			number AlphaHn;
			number AlphaHn_test;
			AlphaHn_test = exp(1.0-0.1*(m_aaVm[*iter] + 65.0))-1.0;
			if (fabs(AlphaHn_test) > m_accuracy)
				AlphaHn = (0.1-0.01*(m_aaVm[*iter] + 65.0)) / AlphaHn_test;
			else
				AlphaHn = 0.1;

			number BetaHn = 0.125*exp((m_aaVm[*iter] + 65.0)/80.0);

			m_aaHGate[*iter] = AlphaHh / (AlphaHh + BetaHh);
			m_aaMGate[*iter] = AlphaHm / (AlphaHm + BetaHm);
			m_aaNGate[*iter] = AlphaHn / (AlphaHn + BetaHn);

		}
	}

}

template<typename TDomain>
void ChannelHH<TDomain>::update_gating(number newTime)
{
	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	//TODO set accuracy later on another point
	number m_accuracy = 1e-6;


	const typename TDomain::position_accessor_type& aaPos = this->m_dom->position_accessor();
	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			number AlphaHh = 0.07*exp(-(m_aaVm[*iter] + 65.0)/20.0);
			number BetaHh = 1.0/(exp(3.0-0.1*(m_aaVm[*iter]  + 65.0))+1.0);

			// gating param m
			number AlphaHm;
			number AlphaHm_test = exp(2.5-0.1*(m_aaVm[*iter] + 65.0))-1.0;
			if (fabs(AlphaHm_test) > m_accuracy)
				AlphaHm = (2.5 - 0.1*( m_aaVm[*iter] + 65.0)) / AlphaHm_test;
			else
				AlphaHm = 1.0;

			number BetaHm = 4.0*exp((-m_aaVm[*iter] + 65.0)/18.0);

			// gating param n
			number AlphaHn;
			number AlphaHn_test;
			AlphaHn_test = exp(1.0-0.1*(m_aaVm[*iter] + 65.0))-1.0;
			if (fabs(AlphaHn_test) > m_accuracy)
				AlphaHn = (0.1-0.01*(m_aaVm[*iter] + 65.0)) / AlphaHn_test;
			else
				AlphaHn = 0.1;

			number BetaHn = 0.125*exp((m_aaVm[*iter] + 65.0)/80.0);

			number rate_h = -((AlphaHh * (1.0-m_aaHGate[*iter])) - BetaHh * m_aaHGate[*iter]);
			number rate_m = -((AlphaHm * (1.0-m_aaMGate[*iter])) - BetaHm * m_aaMGate[*iter]);
			number rate_n = -((AlphaHn * (1.0-m_aaNGate[*iter])) - BetaHn * m_aaNGate[*iter]);

			// needs import later
			number dt = 0.01;

			m_aaHGate[*iter] -= dt * rate_h;
			m_aaMGate[*iter] -= dt * rate_m;
			m_aaNGate[*iter] -= dt * rate_n;

		}
	}


}

template<typename TDomain>
void ChannelHH<TDomain>::ionic_current(Vertex* v, std::vector<number>& outCurrentValues)
{

	const number m_g_K  = 0.36e-3;
	const number m_g_Na = 1.2e-3;
	const number m_g_I  = 0.003e-3;





	// single channel type fluxes
	//const number potassium_part_of_flux = m_g_K * pow(u(_n_,co),4) * (u(_VM_,co) - -77);
	//const number sodium_part_of_flux =  m_g_Na * pow(u(_m_,co),3) * u(_h_,co) * (u(_VM_, co) + -50);
	//const number leakage_part_of_flux = m_g_I * (u(_VM_,co) + 54.4);

	//number flux_value = potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux;

}
////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
void IChannel<TDomain>::
register_all_funcs(bool bHang)
{
	//register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void IChannel<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	//this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	//this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	//this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	//this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);

}



template<typename TDomain>
void ChannelHH<TDomain>::
register_all_funcs(bool bHang)
{
	//register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ChannelHH<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	//this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	//this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	//this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	//this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);

}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class IChannel<Domain1d>;
template class ChannelHH<Domain1d>;
#endif
#ifdef UG_DIM_2
template class IChannel<Domain2d>;
template class ChannelHH<Domain2d>;
#endif
#ifdef UG_DIM_3
template class IChannel<Domain3d>;
template class ChannelHH<Domain3d>;
#endif




} /* namespace ug */
