/*
 * VM_Disc.cpp
 *
 *  Created on: 26.11.2014
 *      Author: Pgottmann
 */

#include "VM_Disc.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/local_transfer_interface.h"





namespace ug {



// Standard Methods for VM-Equatations
template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_diameter(const number d)
{
	// handle the attachment
	if (m_spApproxSpace->domain()->grid()->has_vertex_attachment(m_aDiameter))
		UG_THROW("Radius attachment necessary for HH elem disc "
				 "could not be made, since it already exists.");
	m_spApproxSpace->domain()->grid()->attach_to_vertices_dv(m_aDiameter, d);

	m_aaDiameter = Grid::AttachmentAccessor<Vertex, ANumber>(*m_spApproxSpace->domain()->grid(), m_aDiameter);
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_spec_res(number val)
{
	m_spec_res = val;
}



template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_spec_cap(number val)
{
	m_spec_cap = val;
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_influx(number Flux, number x, number y, number z, number dur, number beg)
{

}


template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
add_channel(SmartPtr<IChannel<TDomain, TAlgebra> > Channel, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct)
{
	std::cout << "Channel_size: " << m_channel.size() << std::endl;
	Channel->init(0.01, approx, spGridFct);
	int Channel_size = m_channel.size();
	m_channel.push_back(Channel);
	std::cout << "Channel added" << std::endl;
	// add Channel function from IElem
}



// Methods for Interface class

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
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
	std::vector<double> AlloutCurrentValues;

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
		number Diam = m_aaDiameter[pElem->vertex(co)];


		volume += scv.volume();
		// add length of scv to element length
		element_length += scv.volume();
		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*Diam*Diam);


		//TODO Influx handling

		//const typename TDomain::position_type& vert_coords = aaPos[pElem->vertex(co)];
		// set flux variable
		double flux = 0;
		/*if 	(vert_coords[0]==flux_coords[0] and vert_coords[1]==flux_coords[1] and vert_coords[2]==flux_coords[2])
			flux = influx;*/


		// loop over all channels
		for (int i = 0; i < m_channel.size(); i++)
		{

			// values we are getting from ionic_flux function in channels
			m_channel[i].get()->ionic_current(pElem->vertex(co), outCurrentValues);
			AlloutCurrentValues[co] += outCurrentValues[co];
		}

		d(_VM_, co) += scv.volume()*PI*Diam *(AlloutCurrentValues[co]);



		//std::cout << "defekt changes: "<< d(_VM_, co) << std::endl;
	}

	// cable equation, "diffusion" part
		/*MathVector<dim> grad_c;

		for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
			// get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			// compute gradient at ip
			VecSet(grad_c, 0.0);
			for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(grad_c, u(_VM_,sh), scvf.global_grad(sh));

			// scalar product with normal
			number diff_flux = VecDot(grad_c, scvf.normal());

			number Diam_FromTo = 0.5 * (m_aaDiameter[pElem->vertex(scvf.from())]
			                          + m_aaDiameter[pElem->vertex(scvf.to())]);

			//calculates pre_resistance
			pre_resistance = volume / (0.25*PI*Diam_FromTo*Diam_FromTo);

			// scale by 1/resistance and by length of element
			diff_flux *= element_length / (m_spec_res*pre_resistance);

			// add to local defect
			d(_VM_, scvf.from()) -= diff_flux;
			d(_VM_, scvf.to()  ) += diff_flux;
		}*/




}

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		//get Diameter from element
		number Diam = 3.18e-6;

		// get spec capacity
		number spec_capacity = 1e-5;


		// potential equation time derivative
		/*std::cout << "u: " << u(_VM_, co) << std::endl;
		std::cout << "co: " << co << std::endl;
		std::cout << "time derivative: " << (PI*Diam*scv.volume()*u(_VM_, co)*spec_capacity) << std::endl;*/

		d(_VM_, co) += PI*Diam*scv.volume()*u(_VM_, co)*spec_capacity;
	}



}

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{



	// get finite volume geometry
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

		number element_length = 0.0;
		number pre_resistance = 0.0;
		number volume = 0;

		// cast elem to appropriate type
		TElem* pElem = dynamic_cast<TElem*>(elem);
		if (!pElem) {UG_THROW("Wrong element type.");}


		std::vector<number> jacFlux;
		std::vector<number> AlljacFlux;

	// channel kinetics derivatives
		for (size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
			// get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

			// get associated node
			const int co = scv.node_id();

			//get Diameter from element later in attachment
			number Diam = 3.18e-6;



			// add length of scv to element length
			element_length += scv.volume();

			// add "pre_resistance" parts
			pre_resistance += scv.volume() / (0.25*PI*Diam*Diam);

			// calculates volume for later use
			volume += scv.volume();

			for (int i = 0;  i < m_channel.size() ; i++)
			{
				// getting jacobian depending on vertex-attachments
				m_channel[i].get()->Jacobi_sets(pElem->vertex(co), jacFlux);
				AlljacFlux[co] += jacFlux[co];
			}


			// derivatives of potential from HH channels
			//J(_VM_, co, _h_, co) += scv.volume()*PI*Diam * m_g_Na*pow(u(_m_,co),3) * (u(_VM_, co) - m_sodium);
			//J(_VM_, co, _m_, co) += scv.volume()*PI*Diam * 3.0*m_g_Na*pow(u(_m_,co),2) * u(_h_,co) * (u(_VM_, co) - m_sodium);
			//J(_VM_, co, _n_, co) += scv.volume()*PI*Diam * 4.0*m_g_K*pow(u(_n_,co),3) * (u(_VM_,co) + m_potassium);
			//TODO Has to be settet elsewhere
			J(_VM_, co, _VM_, co) += scv.volume()*PI*Diam*(AlljacFlux[co]); //* (m_g_K*pow(u(_n_,co),4) + m_g_Na*pow(u(_m_,co),3)*u(_h_,co) + m_g_I);
		}


		//diffusive part
		/*for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
			{
				// get current SCVF
				const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

				// loop shape functions
				for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					// scalar product with normal
					number Diam_FromTo = 0.5 * (m_aaDiameter[pElem->vertex(scvf.from())]
					                          + m_aaDiameter[pElem->vertex(scvf.to())]);

					pre_resistance = volume / (0.25*PI*Diam_FromTo*Diam_FromTo);


					number d_diff_flux = VecDot(scvf.global_grad(sh), scvf.normal());

					// scale by 1/resistance and by length of element
					d_diff_flux *= element_length / (m_spec_res*pre_resistance);

					// add flux term to local matrix
					J(_VM_, scvf.from(), _VM_, sh) -= d_diff_flux;
					J(_VM_, scvf.to()  , _VM_, sh) += d_diff_flux;
				}
			}*/
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		//get Diameter from element later in attachment
		number Diam = 3.18e-6;

		//spec_capa has to be set later on in an varialbe

		// get spec capacity
		number spec_capacity = 1e-5;



		// potential equation
		J(_VM_, co, _VM_, co) += PI*Diam*scv.volume()*spec_capacity;
	}
}


template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::fsh_elem_loop()
{


}


template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("Cannot update Finite Volume Geometry.\n");

	//we coud update vm here...

}


template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{



}


template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_funcs(m_bNonRegularGrid);

	//TODO here we need some option to get the number of unknowns
}



template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}





////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::
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



////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////


#ifdef UG_DIM_1
#ifdef UG_CPU_1
template class VMDisc<Domain1d, CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class VMDisc<Domain1d, CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class VMDisc<Domain1d, CPUBlockAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class VMDisc<Domain1d, CPUBlockAlgebra<4> >;
#endif
#ifdef UG_CPU_VAR
template class VMDisc<Domain1d, CPUVariableBlockAlgebra >;
#endif
#endif



#ifdef UG_DIM_2
#ifdef UG_CPU_1
template class VMDisc<Domain2d, CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class VMDisc<Domain2d, CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class VMDisc<Domain2d, CPUBlockAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class VMDisc<Domain2d, CPUBlockAlgebra<4> >;
#endif
#ifdef UG_CPU_VAR
template class VMDisc<Domain2d, CPUVariableBlockAlgebra >;
#endif
#endif


#ifdef UG_DIM_3
#ifdef UG_CPU_1
template class VMDisc<Domain3d, CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class VMDisc<Domain3d, CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class VMDisc<Domain3d, CPUBlockAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class VMDisc<Domain3d, CPUBlockAlgebra<4> >;
#endif
#ifdef UG_CPU_VAR
template class VMDisc<Domain3d, CPUVariableBlockAlgebra >;
#endif
#endif




} /* namespace ug */
