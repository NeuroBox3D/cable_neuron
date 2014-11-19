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
void IChannel<TDomain>::add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void IChannel<TDomain>::add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void IChannel<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void IChannel<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}

template<typename TDomain>
IChannel<TDomain>::~IChannel()
{

}

template<typename TDomain>
template <typename TElem, typename TFVGeom>
void IChannel<TDomain>::fsh_elem_loop()
{


}


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
template <typename TElem, typename TFVGeom>
void IChannel<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{



}


template<typename TDomain>
void IChannel<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_funcs(m_bNonRegularGrid);

	//TODO here we need some option to get the number of unknowns
}

template<typename TDomain>
template <typename TElem, typename TFVGeom>
void IChannel<TDomain>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// some help constants
	number element_length = 0.0;
	number pre_resistance = 0.0;
	number volume = 0.0;

	//need to set later in another way also given from elsewhere
	number Na_out = 140;
	number K_out = 2.5;

	// helper for saving defekts
	std::vector<double> outCurrentValues;

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		//TODO get Diam from somewhere
		number Diam = 0.25; //m_aaDiameter[pElem->vertex(co)];


		volume += scv.volume();
		// add length of scv to element length
		element_length += scv.volume();
		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*Diam*Diam);



		// values we are getting from ionic_flux function

		ionic_current(pElem->vertex(co), outCurrentValues);

		//need to add all fluxes for every single element in hh case only VM-Flux
		d(_VM_, co) += scv.volume()*PI*Diam*(outCurrentValues[co]);


		// add defekt of gatting Params
		//d(_h_, co) += outCurrentValues[0];
		//d(_m_, co) += outCurrentValues[1];
		//d(_n_, co) += outCurrentValues[2];
	}



}


/* IChannel:: use_hanging()
{
	return false;
}*/


// Methods for HH-Channel-Class


template<typename TDomain>
ChannelHH<TDomain>::~ChannelHH()
{

}


template<typename TDomain>
void ChannelHH<TDomain>::init(number time, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<LocalVector> u)
{
	// attach attachments
	if (this->m_mg->has_vertex_attachment(this->m_MGate))
		UG_THROW("Attachment necessary (MGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_MGate);

	if (this->m_mg->has_vertex_attachment(this->m_HGate))
		UG_THROW("Attachment necessary (HGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
		this->m_mg->attach_to_vertices(this->m_HGate);

	if (this->m_mg->has_vertex_attachment(this->m_NGate))
		UG_THROW("Attachment necessary (NGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_NGate);


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
	SmartPtr<DoFDistribution> dd;

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
			dd = (approx->dof_distribution(GridLevel::TOP));
			const LFEID& lfeID = dd->local_finite_element_id(fct);
			//const
			//getting exakt index of value
			dd->inner_dof_indices(*iter, fct, multInd);
			// has to be the solvung values
			number test = multInd[0][0];
			this->m_aaVm[*iter] = ;

			std::cout<< "getting VM is working in false way: " << m_aaVm[*iter] << std::endl;

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

			// Setting Starting gatting params independent on VM
			this->m_aaHGate[*iter] = AlphaHh / (AlphaHh + BetaHh);
			this->m_aaMGate[*iter] = AlphaHm / (AlphaHm + BetaHm);
			this->m_aaNGate[*iter] = AlphaHn / (AlphaHn + BetaHn);


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
			std::cout << "Rates: " << rate_h << " , " << rate_m << " , " << rate_n << std::endl;

		}
	}


}

template<typename TDomain>
void ChannelHH<TDomain>::ionic_current(Vertex* v, std::vector<number>& outCurrentValues)
{

	// later need to be set somewhere else
	const number m_g_K  = 0.36e-3;
	const number m_g_Na = 1.2e-3;
	const number m_g_I  = 0.003e-3;


	// getting attachments out of Vertex
	double NGate = m_aaNGate[v];
	double MGate = m_aaMGate[v];
	double HGate = m_aaHGate[v];
	double VM 	 = m_aaVm[v];


	// TODO Influx values needed
	// single channel type fluxes
	const number potassium_part_of_flux = m_g_K * pow(NGate,4) * (VM + 77);
	const number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (VM - 50);
	const number leakage_part_of_flux = m_g_I * (VM + 54.4);

	number flux_value = potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux;

	outCurrentValues.push_back(flux_value);

}
////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
void IChannel<TDomain>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void IChannel<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);

}



template<typename TDomain>
void ChannelHH<TDomain>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
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
