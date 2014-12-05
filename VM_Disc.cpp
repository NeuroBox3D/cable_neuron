/*
 * VM_Disc.cpp
 *
 *  Created on: 26.11.2014
 *      Author: Pgottmann
 *
 *
 *      Discretization of Kabelequatation depending on function called _VM_ needed
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
	if (m_spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(m_aDiameter))
		UG_THROW("Radius attachment necessary for HH elem disc "
				 "could not be made, since it already exists.");
	m_spGridFct->approx_space()->domain()->grid()->attach_to_vertices_dv(m_aDiameter, d);

	m_aaDiameter = Grid::AttachmentAccessor<Vertex, ADouble>(*m_spGridFct->approx_space()->domain()->grid(), m_aDiameter);
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_spec_res(number val)
{
	m_spec_res = val;
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_influx_ac(number influx_ac)
{
	m_influx_ac = influx_ac;
}


template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_spec_cap(number val)
{
	m_spec_cap = val;
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_influx(number Flux, number x, number y, number z, number beg, number dur)
{

	int flux_number = m_flux_value.size();

	m_flux_value.push_back(Flux);

	m_beg_flux.push_back(beg);
	m_dur_flux.push_back(dur);

	m_coords.push_back((x, y, z));

}


template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
add_channel(SmartPtr<IChannel<TDomain, TAlgebra> > Channel, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct)
{
	std::cout << "Channel_size: " << m_channel.size() << std::endl;
	Channel->init(0.01, spGridFct);
	int Channel_size = m_channel.size();
	m_channel.push_back(Channel);
	std::cout << "Channel added" << std::endl;
	// add Channel function from IElem
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
add_func(const char* func)
{
	m_funcs.push_back(func);
	m_numb_funcs += 1;
}



// Methods for Interface class

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

	std::cout << "deff a elem starts" << std::endl;

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

		// get Diam from attachment
		number Diam = m_aaDiameter[pElem->vertex(co)];

		// calculate volume for later use
		volume += scv.volume();
		// add length of scv to element length
		element_length += scv.volume();
		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*Diam*Diam);

		// Helpers for Influx-Handling
		number time = this->time();
		double influx = 0;


	///Influx handling
		for (int i = 0; i < m_flux_value.size(); i++)
		{
			/*std::cout << "coords: " << m_coords[i][0] << " - " << vCornerCoords[0][1] << std::endl;
			std::cout << "times: " << time << " - " << m_beg_flux[i] << std::endl;
			std::cout << "echtes erg: " << (vCornerCoords[0][2] - m_coords[i][2]) << std::endl;
			std::cout << "abs erg: " << ((fabs((vCornerCoords[0][2] - m_coords[i][2])))) << std::endl;*/
			// Time depending vars
			if ((m_beg_flux[i] <= time)
			and ((m_dur_flux[i] + m_beg_flux[i]) >= time)
			// place depending vars
			and ((fabs(vCornerCoords[0][0] - m_coords[i][0])) < m_influx_ac)
			and ((fabs(vCornerCoords[0][1] - m_coords[i][1])) < m_influx_ac)
			and ((fabs(vCornerCoords[0][2] - m_coords[i][2])) < m_influx_ac))
			{
				//std::cout << "influx should happening" << std::endl;
				influx = m_flux_value[i];
				std::cout << "influx value: " << influx << std::endl;
			}

		}

		// for all functions space is needed
		for (int i=1; i<m_numb_funcs; i++)
		{
			AlloutCurrentValues.push_back(0);
		}

	/// Channel defekt adding
		for (int i = 0; i < m_channel.size(); i++)
		{
			// values we are getting from ionic_flux function in channels
			m_channel[i].get()->ionic_current(pElem->vertex(co), outCurrentValues);
			// adding defekt from every channel
			// in 0 all Vms
			for (int j=0; j<m_numb_funcs; j++)
				{
					AlloutCurrentValues[j] += (outCurrentValues[j]);
				}
		}




	/// Writing all into defekts
		d(_VM_, co) += scv.volume()*PI*Diam *(AlloutCurrentValues[0]-influx);

	/// Now all other defektes needed to be added
		for (int k=1; k < m_numb_funcs; k++)
		{
			d(k, co) += scv.volume()*PI*Diam*AlloutCurrentValues[k];
		}

		//std::cout << "defekt changes: "<< d(_VM_, co) << std::endl;
	}

	// cable equation, "diffusion" part
		MathVector<dim> grad_c;

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
		}
		std::cout << "def a elem ends" << std::endl;




}

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	std::cout << "add def m elem start" << std::endl;
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
		number Diam = m_aaDiameter[pElem->vertex(co)];

		// get spec capacity
		number spec_capacity = m_spec_cap;


		// potential equation time derivative
		/*std::cout << "u: " << u(_VM_, co) << std::endl;
		std::cout << "co: " << co << std::endl;
		std::cout << "time derivative: " << (PI*Diam*scv.volume()*u(_VM_, co)*spec_capacity) << std::endl;*/

		// Nernst Paras time derivative
		d(_Na_, co) += u(_Na_, co)*scv.volume()*0.25*PI*Diam*Diam;
		d(_K_, co)  += u(_K_, co)*scv.volume()*0.25*PI*Diam*Diam;

		d(_VM_, co) += PI*Diam*scv.volume()*u(_VM_, co)*spec_capacity;
	}


	std::cout << "add def m elem end" << std::endl;
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
		number Diam = m_aaDiameter[pElem->vertex(co)];



		// add length of scv to element length
		element_length += scv.volume();

		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*Diam*Diam);

		// calculates volume for later use
		volume += scv.volume();


		/*	ionic_current is constant in unknowns (new values for V_m, Na, K)
		// for all functions
		for (int i=1; i<m_numb_funcs; i++)
		{
			AlljacFlux.push_back(0);
		}


	/// Jacobi Matrix-Channel handling
		for (int i = 0;  i < m_channel.size() ; i++)
		{
			// getting jacobian depending on vertex-attachments
			m_channel[i].get()->Jacobi_sets(pElem->vertex(co), jacFlux);
			//adding jacobian from every channel
			for (int i=0; i<m_numb_funcs; i++)
			{
				AlljacFlux[i] += (jacFlux[i]);
			}
		}

		for (int i=0; i<m_numb_funcs; i++)
		{
			J(i, co, i, co) += scv.volume()*PI*Diam*(AlljacFlux[i]);
		}
		 */
	}


	//diffusive part
	for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
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
	}
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	std::cout << "jac m elem starts" << std::endl;
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
		number Diam = m_aaDiameter[pElem->vertex(co)];

		//spec_capa has to be set later on in an varialbe

		// get spec capacity
		number spec_capacity = m_spec_cap;

		J(_Na_, co, _Na_, co) += 1.0*scv.volume()*0.25*PI*Diam*Diam;
		J(_K_, co, _K_, co) += 1.0*scv.volume()*0.25*PI*Diam*Diam;

		// potential equation
		J(_VM_, co, _VM_, co) += PI*Diam*scv.volume()*spec_capacity;
	}
	std::cout << "jac m elem ends" << std::endl;
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

	//std::cout << "unumber of all functions " << u.num_all_fct() << std::endl;
	//number function_number = u.num_all_fct();
	//_VM_ = u.fct_id_by_name(m_funcs[0]);

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

	// check number
	if (vLfeID.size() < 1)
		UG_THROW("VMDisc: Wrong number of functions given. Need min 1 function "<< 1);

	if (vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("VMDISC FV Scheme only implemented for 1st order.");

	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_funcs(m_bNonRegularGrid);



	std::cout << "before prepare" << std::endl;
	//VM always needed in this diskretication so it is easy only getting index of VM
	for (int i = 0; i < m_numb_funcs; i++)
	{
		if (m_funcs[i] = "VM")
			_VM_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
		if (m_funcs[i] = "K")
			_K_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
		if (m_funcs[i] = "Na")
			_Na_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
	}
	std::cout << "after prepare" << std::endl;

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
