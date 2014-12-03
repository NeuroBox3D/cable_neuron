/*
 * channel_interface.cpp
 *
 *  Created on: 29.10.2014
 *      Author: mbreit
 */

#include "channel_interface.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/local_transfer_interface.h"




namespace ug {

//////////////////////////////////
// Methods for Interface class
//////////////////////////////////

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void IChannel<TDomain, TAlgebra>::prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("Cannot update Finite Volume Geometry.\n");

}


template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void IChannel<TDomain, TAlgebra>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{



}

template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void IChannel<TDomain, TAlgebra>::fsh_elem_loop()
{

}


template<typename TDomain, typename TAlgebra>
void IChannel<TDomain, TAlgebra>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_funcs(m_bNonRegularGrid);


}

template<typename TDomain, typename TAlgebra>
bool IChannel<TDomain, TAlgebra>::
use_hanging() const
{
	// As this is basically a 1D discretization,
	// there will never be any hanging nodes.
	return true;
}



template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void IChannel<TDomain, TAlgebra>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}



////////////////////////////////////////////////
// Methods for HH-Channel-Class
////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::
set_accuracy(double ac)
{
	m_accuracy = ac;
}


template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::
set_consts(number Na, number K, number L)
{
	m_g_K = K;
	m_g_Na = Na;
	m_g_I = L;
}

template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::
set_rev_pot(number R_Na, number R_K)
{
  m_sodium = R_Na;
  m_potassium = R_K;

}




// Methods for using gatings
template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::init(number time, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct)
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


	// creates Multi index
	std::vector<DoFIndex> multInd;
	SmartPtr<DoFDistribution> dd;
	//SmartPtr<vector_type> spU;



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

			// needed konzentration has to be set
			size_t VM = spGridFct->fct_id_by_name("VM");
			//std::cout << "VM: " << VM << std::endl;


			dd = (approx->dof_distribution(GridLevel::TOP));
			const LFEID& lfeID = dd->local_finite_element_id(VM);
			//const
			//getting exakt index of value
			//dd->inner_dof_indices(*iter, fct, multInd);
			// has to be the solvung values

			number test = multInd[0][0];
			number sub;// = spGridFct[test];

			Vertex* vrt = *iter;


			spGridFct->dof_indices(vrt, VM, multInd);


			for (int i=0; i<multInd.size(); i++)
			{
				m_aaVm[*iter] = DoFRef(*spGridFct, multInd[i]);
			}




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

			// Setting Starting gatting params independent on VM
			this->m_aaHGate[*iter] = AlphaHh / (AlphaHh + BetaHh);
			this->m_aaMGate[*iter] = AlphaHm / (AlphaHm + BetaHm);
			this->m_aaNGate[*iter] = AlphaHn / (AlphaHn + BetaHn);

		}
	}

}

template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::update_gating(number newTime, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct)
{
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

			//getting act vertex
			Vertex* vrt = *iter;

			// needed konzentration has to be set
			size_t VM = spGridFct->fct_id_by_name("VM");
			//std::cout << "VM: upg" << VM << std::endl;


			std::vector<DoFIndex> multInd;
			spGridFct->dof_indices(vrt, VM, multInd);

			//update Vm
			for (int i=0; i<multInd.size(); i++)
			{
				m_aaVm[*iter] = DoFRef(*spGridFct, multInd[i]);
			}

			//set gatting Params
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
void ChannelHH<TDomain, TAlgebra>::ionic_current(Vertex* v, std::vector<number>& outCurrentValues)
{
	// getting attachments out of Vertex

	double NGate = m_aaNGate[v];
	double MGate = m_aaMGate[v];
	double HGate = m_aaHGate[v];
	double VM 	 = m_aaVm[v];



	// TODO Influx values needed
	// single channel type fluxes
	const number potassium_part_of_flux = m_g_K * pow(NGate,4) * (VM + m_potassium);
	const number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (VM - m_sodium);
	const number leakage_part_of_flux = m_g_I * (VM + 54.4);

	number flux_value = (potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux);
	//std::cout << "Flux_Value: " << flux_value << std::endl;
	outCurrentValues.push_back(flux_value);

}



template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::Jacobi_sets(Vertex* v, std::vector<number>& outJFlux)
{
	double NGate = m_aaNGate[v];
	double MGate = m_aaMGate[v];
	double HGate = m_aaHGate[v];


	number Jac = (m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);

	outJFlux.push_back(Jac);

}

template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void ChannelHH<TDomain, TAlgebra>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{



}


////////////////////////////////////////////////
// Methods for HH-Channel-Nernst- Class
////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::
set_accuracy(double ac)
{
	m_accuracy = ac;
}


template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::
set_consts(number Na, number K, number L)
{
	m_g_K = K;
	m_g_Na = Na;
	m_g_I = L;
}

template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::
set_nernst_consts(number R, number T, number F)
{
	//std::cout << "using nernst consts" << std::endl;
	m_R = R;
	m_T = T;
	m_F = F;

}

template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::
set_out_conc(number Na, number K)
{
	m_Na_out = Na;
	m_K_out  = K;
}

// Methods for using gatings
template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::init(number time, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct)
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

	if (this->m_mg->has_vertex_attachment(this->m_Na))
		UG_THROW("Attachment necessary (HGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
		this->m_mg->attach_to_vertices(this->m_Na);

	if (this->m_mg->has_vertex_attachment(this->m_K))
		UG_THROW("Attachment necessary (NGate) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_K);

	if (this->m_mg->has_vertex_attachment(this->m_Vm))
		UG_THROW("Attachment necessary (Vm) for Hodgkin and Huxley channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_Vm);


	// create attachment accessors
	this->m_aaMGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_MGate);
	this->m_aaNGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_NGate);
	this->m_aaHGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_HGate);
	this->m_aaVm = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_Vm);
	this->m_aaNa = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_Na);
	this->m_aaK = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_K);


	// creates Multi index
	std::vector<DoFIndex> multInd, multIndK, multIndNa;
	SmartPtr<DoFDistribution> dd;
	//SmartPtr<vector_type> spU;



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

			// needed konzentration has to be set
			size_t VM = spGridFct->fct_id_by_name("VM");
			size_t K = spGridFct->fct_id_by_name("K");
			size_t Na = spGridFct->fct_id_by_name("Na");
			//std::cout << "VM: " << VM << std::endl;


			dd = (approx->dof_distribution(GridLevel::TOP));
			const LFEID& lfeID = dd->local_finite_element_id(VM);
			//const
			//getting exakt index of value
			//dd->inner_dof_indices(*iter, fct, multInd);
			// has to be the solvung values

			//number test = multInd[0][0];
			number sub;// = spGridFct[test];

			Vertex* vrt = *iter;


			spGridFct->dof_indices(vrt, VM, multInd);
			spGridFct->dof_indices(vrt, K, multIndK);
			spGridFct->dof_indices(vrt, Na, multIndNa);


			//for (int i=0; i<multInd.size(); i++)
			//{
				m_aaVm[*iter] = DoFRef(*spGridFct, multInd[0]);
				m_aaK[*iter] = DoFRef(*spGridFct, multIndK[0]);
				m_aaNa[*iter] = DoFRef(*spGridFct, multIndNa[0]);
			//}




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

			// Setting Starting gatting params independent on VM
			this->m_aaHGate[*iter] = AlphaHh / (AlphaHh + BetaHh);
			this->m_aaMGate[*iter] = AlphaHm / (AlphaHm + BetaHm);
			this->m_aaNGate[*iter] = AlphaHn / (AlphaHn + BetaHn);

		}
	}

}

template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::update_gating(number newTime, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct)
{
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

			//getting act vertex
			Vertex* vrt = *iter;

			// needed konzentration has to be set
			size_t VM = spGridFct->fct_id_by_name("VM");
			size_t Na = spGridFct->fct_id_by_name("Na");
			size_t K = spGridFct->fct_id_by_name("K");
			//std::cout << "VM: upg" << VM << std::endl;


			std::vector<DoFIndex> multInd, multIndK, multIndNa;

			spGridFct->dof_indices(vrt, VM, multInd);
			spGridFct->dof_indices(vrt, K, multIndK);
			spGridFct->dof_indices(vrt, Na, multIndNa);

			// updating all Vars
			//for (int i=0; i<multInd.size(); i++)
			//{
				m_aaVm[*iter] = DoFRef(*spGridFct, multInd[0]);
				m_aaK[*iter] = DoFRef(*spGridFct, multIndK[0]);
				m_aaNa[*iter] = DoFRef(*spGridFct, multIndNa[0]);
			//}

			//set gatting Params
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

			// needs import later
			number dt = newTime - this->time();//0.01;

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
void ChannelHHNernst<TDomain, TAlgebra>::ionic_current(Vertex* v, std::vector<number>& outCurrentValues)
{
	// getting attachments out of Vertex
	double NGate = m_aaNGate[v];
	double MGate = m_aaMGate[v];
	double HGate = m_aaHGate[v];
	double VM 	 = m_aaVm[v];
	double K 	 = m_aaK[v];
	double Na	 = m_aaNa[v];

	const number helpV = (m_R*m_T)/m_F;
	// nernst potential of potassium and sodium
	//std::cout << "K: " << K << std::endl;

	const number potassium_nernst_eq 	= helpV*(log(m_K_out/K));
	const number sodium_nernst_eq	 	= -helpV*(log(m_Na_out/Na));


	//std::cout << "potassium_nernst_eq: " << potassium_nernst_eq << std::endl;

	// single channel type fluxes
	const number potassium_part_of_flux = m_g_K * pow(NGate,4) * (VM - potassium_nernst_eq);
	const number sodium_part_of_flux =  m_g_Na * pow(MGate,3) * HGate * (VM + sodium_nernst_eq);
	const number leakage_part_of_flux = m_g_I * (VM + 54.4);

	//std::cout << "potassium_part_of_flux: " << potassium_part_of_flux << std::endl;

	outCurrentValues.push_back(potassium_part_of_flux + sodium_part_of_flux + leakage_part_of_flux);
	outCurrentValues.push_back(sodium_part_of_flux/m_F);
	outCurrentValues.push_back(potassium_part_of_flux/m_F);

	//std::cout << "outCurrentValues: " << outCurrentValues[0] << ", " << outCurrentValues[1] << ", " << outCurrentValues[2] << ", " << std::endl;

}



template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::Jacobi_sets(Vertex* v, std::vector<number>& outJFlux)
{
	double NGate = m_aaNGate[v];
	double MGate = m_aaMGate[v];
	double HGate = m_aaHGate[v];
	double K 	 = m_aaK[v];
	double Na	 = m_aaNa[v];

	const number helpV = (m_R*m_T)/m_F;

	const number potassium_nernst_eq_dK 	=  helpV * (-m_K_out/K)*0.18; //helpV * (-K_out/pow(u(_K_,co),2));
	const number sodium_nernst_eq_dNa		=  -helpV * (-m_Na_out/Na)*0.003;

	outJFlux.push_back(m_g_K*pow(NGate,4) + m_g_Na*pow(MGate,3)*HGate + m_g_I);
	outJFlux.push_back(sodium_nernst_eq_dNa);
	outJFlux.push_back(potassium_nernst_eq_dK);

	//std::cout << "outJFlux: " << outJFlux[0] << ", " << outJFlux[1] << ", " << outJFlux[2] << ", " << std::endl;
}

template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void ChannelHHNernst<TDomain, TAlgebra>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{



}











////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////
// Interface Reg
////////////////////////////
template<typename TDomain, typename TAlgebra>
void IChannel<TDomain, TAlgebra>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void IChannel<TDomain, TAlgebra>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);
	//this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	//this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	//this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	//this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);

}


////////////////////////////////////
// Hodkin and Huxley Reg
////////////////////////////////////
template<typename TDomain, typename TAlgebra>
void ChannelHH<TDomain, TAlgebra>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void ChannelHH<TDomain, TAlgebra>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	//this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	//this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	//this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);

}

////////////////////////////////////
// Hodkin and Huxley Reg with Nernst Reg
////////////////////////////////////
template<typename TDomain, typename TAlgebra>
void ChannelHHNernst<TDomain, TAlgebra>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void ChannelHHNernst<TDomain, TAlgebra>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	//this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	//this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	//this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);

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
template class ChannelHHNernst<Domain1d, CPUAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class IChannel<Domain1d, CPUBlockAlgebra<3> >;
template class ChannelHH<Domain1d, CPUBlockAlgebra<3> >;
template class ChannelHHNernst<Domain1d, CPUAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class IChannel<Domain1d, CPUBlockAlgebra<4> >;
template class ChannelHH<Domain1d, CPUBlockAlgebra<4> >;
template class ChannelHHNernst<Domain1d, CPUAlgebra<4> >;
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
