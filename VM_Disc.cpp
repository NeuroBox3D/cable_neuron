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
	if (m_spApproxSpace->domain()->grid()->has_vertex_attachment(m_aDiameter))
		UG_THROW("Radius attachment necessary for Vm disc "
				 "could not be made, since it already exists.");
	m_spApproxSpace->domain()->grid()->attach_to_vertices_dv(m_aDiameter, d);

	m_aaDiameter = Grid::AttachmentAccessor<Vertex, ADouble>(*m_spApproxSpace->domain()->grid(), m_aDiameter);
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_spec_res(number val)
{
	m_spec_res = val;
}


template<typename TDomain, typename TAlgebra>
SmartPtr<ApproximationSpace<TDomain> > VMDisc<TDomain, TAlgebra>::approx_space()
{
	return m_spApproxSpace;
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
set_diff_coeffs(const std::vector<number>& diff_coeffs)
{
	m_diff = diff_coeffs;
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
set_influx(number Flux, number x, number y, number z, number beg, number dur)
{
	m_flux_value.push_back(Flux);

	m_beg_flux.push_back(beg);
	m_dur_flux.push_back(dur);

	// TODO: a little bit ugly; maybe save influx positions differently?
	m_coords.resize(m_coords.size()+1);
	m_coords[m_coords.size()-1][0] = x;
	if (dim >= 2) m_coords[m_coords.size()-1][1] = y;
	if (dim >= 3) m_coords[m_coords.size()-1][2] = z;
}



#if 0
template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
add_channel(SmartPtr<IChannel<TDomain, TAlgebra> > Channel)
{
	m_channel.push_back(Channel);
}

template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
add_func(std::string func)
{
	m_funcs.push_back(func);
	m_numb_funcs += 1;
}
#endif


template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::
update(number dt, ConstSmartPtr<typename TAlgebra::vector_type> uOld)
{
	// save old solution for explicit assembling later
	m_uOld = uOld;

	// update channels
	for (size_t i = 0; i < m_channel.size(); ++i)
		m_channel[i]->update_gating(dt, uOld);
}


// Methods for Interface class

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// some helper vars
	number element_length = 0.0;
	number pre_resistance = 0.0;

	// cast elem to appropriate type (in order to allow access to attachments)
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	// membrane transport mechanisms and forced influx
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// get diam from attachment
		number Diam = m_aaDiameter[pElem->vertex(co)];

		// add length of scv to element length
		element_length += scv.volume();

		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*Diam*Diam);


		// influx handling
		number time = this->time();
		double influx = 0;
		for (size_t i = 0; i < m_flux_value.size(); i++)
		{
			/*std::cout << "coords: " << m_coords[i][0] << " - " << vCornerCoords[0][1] << std::endl;
			std::cout << "times: " << time << " - " << m_beg_flux[i] << std::endl;
			std::cout << "echtes erg: " << (vCornerCoords[0][2] - m_coords[i][2]) << std::endl;
			std::cout << "abs erg: " << ((fabs((vCornerCoords[0][2] - m_coords[i][2])))) << std::endl;*/
			// Time depending vars
			if (m_beg_flux[i] <= time && m_dur_flux[i] + m_beg_flux[i] >= time
				&& fabs(vCornerCoords[co][0] - m_coords[i][0]) < m_influx_ac
				&& fabs(vCornerCoords[co][1] - m_coords[i][1]) < m_influx_ac
				&& fabs(vCornerCoords[co][2] - m_coords[i][2]) < m_influx_ac
			   )
			{
				influx += m_flux_value[i];
			}
		}

		// membrane transport mechanisms
		std::vector<number> allOutCurrentValues;
		for (size_t i = 0; i < m_numb_funcs+1; ++i)
			allOutCurrentValues.push_back(0.0);

		for (size_t i = 0; i < m_channel.size(); i++)
		{
			std::vector<number> outCurrentValues;

			// values we are getting from ionic_flux function in channels
			std::vector<number> vrt_values(m_numb_funcs+1);
			//for (size_t j = 0; j < m_numb_funcs+1; ++j) vrt_values[j] = u(j, co); <-- NO! this would be implicit!
			std::vector<DoFIndex> multInd;
			for (size_t j = 0; j < m_numb_funcs+1; ++j)
			{
				const size_t fct_ind = m_spDD->fct_id_by_name(this->m_vFct[j].c_str());
				m_spDD->dof_indices(pElem->vertex(co), fct_ind, multInd);
				UG_ASSERT(multInd.size() == 1, "multi-index has size != 1");
				vrt_values[j] = DoFRef(*m_uOld, multInd[0]);
			}

			m_channel[i]->ionic_current(pElem->vertex(co), vrt_values, outCurrentValues);

			const std::vector<std::string>& functions = m_channel[i]->write_fcts();

			// adding defect for every ion species involved
			for (size_t j = 0; j < outCurrentValues.size(); j++)
				allOutCurrentValues[get_index(functions[j])] += (outCurrentValues[j]);
		}

		// writing potential defects
		d(_VM_, co) += scv.volume()*PI*Diam * (allOutCurrentValues[0]-influx);

		// writing ion species defects
		for (size_t k = 1; k < m_numb_funcs+1; k++)
			d(k, co) += scv.volume()*PI*Diam * allOutCurrentValues[k];
	}

	// diffusive parts
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
		number grad_normal = VecDot(grad_c, scvf.normal());

		// scale by 1/resistance and by length of element
		number diff_flux = grad_normal * element_length / (m_spec_res*pre_resistance);

		if (diff_flux > 1.0 || diff_flux < -1.0)
		{
			UG_LOG("m_spec_res: " << m_spec_res << "   pre_res: " << pre_resistance << std::endl);
		}

		// add to local defect of VM
		d(_VM_, scvf.from()) -= diff_flux;
		d(_VM_, scvf.to()  ) += diff_flux;

		// diameter of axial flux cross-section
		number diam_fromTo = std::min(m_aaDiameter[pElem->vertex(scvf.from())],
							   m_aaDiameter[pElem->vertex(scvf.to())]);

		for (size_t k = 1; k < m_numb_funcs+1; k++)
		{
			// compute gradient at ip
			VecSet(grad_c, 0.0);
			for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(grad_c, u(k,sh), scvf.global_grad(sh));

			// scalar product with normal
			grad_normal = VecDot(grad_c, scvf.normal());

			// scale by cross section and diff const
			diff_flux = grad_normal * m_diff[k-1] * 0.25*PI * diam_fromTo*diam_fromTo;


			if (diff_flux > 1.0 || diff_flux < -1.0)
			{
				UG_LOG("m_diff[" << k-1 << "]: " << m_diff[k-1] << std::endl);
			}

			d(k, scvf.from()) -= diff_flux;
			d(k, scvf.to()  ) += diff_flux;
		}

	}
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
		number diam = m_aaDiameter[pElem->vertex(co)];

		// get spec capacity
		number spec_capacity = m_spec_cap;


		// potential equation time derivative
		d(_VM_, co) += PI*diam*scv.volume()*u(_VM_, co)*spec_capacity;

		// ion species time derivative
		for (size_t k = 1; k < m_numb_funcs+1; k++)
			d(k, co) += u(k, co)*scv.volume()*0.25*PI*diam*diam;
	}
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// some helper vars
	number element_length = 0.0;
	number pre_resistance = 0.0;

	// cast elem to appropriate type  (in order to allow access to attachments)
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	// only helper calculations for axial current here (membrane fluxes are purely explicit)
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// get diam from attachment
		number Diam = m_aaDiameter[pElem->vertex(co)];

		// add length of scv to element length
		element_length += scv.volume();

		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*Diam*Diam);
	}


	// diffusive part
	for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		// loop shape functions
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			// scalar product with normal
			number grad_normal = VecDot(scvf.global_grad(sh), scvf.normal());

			// scale by 1/resistance and by length of element
			number d_diff_flux = grad_normal * element_length / (m_spec_res*pre_resistance);

			// add flux term to local matrix
			J(_VM_, scvf.from(), _VM_, sh) -= d_diff_flux;
			J(_VM_, scvf.to()  , _VM_, sh) += d_diff_flux;

			// diameter of axial flux cross-section
			number diam_fromTo = std::min(m_aaDiameter[pElem->vertex(scvf.from())],
								   m_aaDiameter[pElem->vertex(scvf.to())]);

			for (size_t k = 1; k < m_numb_funcs+1; k++)
			{
				// scale by cross section and diff const
				d_diff_flux = grad_normal * m_diff[k-1] * 0.25*PI * diam_fromTo*diam_fromTo;
				J(k, scvf.from(), k, sh) -= d_diff_flux;
				J(k, scvf.to(), k, sh) += d_diff_flux;
			}
		}
	}
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//std::cout << "jac m elem starts" << std::endl;
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
		for (size_t k = 1; k < m_numb_funcs+1; k++)
		{
			J(k, co, k, co) += scv.volume()*0.25*PI*Diam*Diam;
		}
		// potential equation
		J(_VM_, co, _VM_, co) += PI*Diam*scv.volume()*spec_capacity;
	}
	//std::cout << "jac m elem ends" << std::endl;
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
	try
	{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("Cannot update Finite Volume Geometry.\n");
}


template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{



}



template<typename TDomain, typename TAlgebra>
void VMDisc<TDomain, TAlgebra>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check number
	if (vLfeID.size() != 4)
		UG_THROW("VMDisc: Wrong number of functions given. Need exactly 4 functions ");

	if (vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("VMDisc FV scheme only implemented for 1st order.");

	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_funcs(m_bNonRegularGrid);

	//std::cout << "before prepare" << std::endl;
	//VM always needed in this discretization so it is easy only getting index of VM
	/*
	for (int i = 0; i < m_numb_funcs; i++)
	{
		if (m_funcs[i] == "VM")
			_VM_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
		if (m_funcs[i] == "K")
			_K_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
		if (m_funcs[i] == "Na")
			_Na_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
	}
	//std::cout << "after prepare" << std::endl;
	*/
}

template<typename TDomain, typename TAlgebra>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain, TAlgebra>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// nothing to do
}


////////////////////////////////////////////////////////////////////////////////
// Functions for ion handling
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain, typename TAlgebra>
size_t VMDisc<TDomain, TAlgebra>::get_index(std::string s)
{
	return m_spApproxSpace->fct_id_by_name(s.c_str());
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
	typedef VMDisc<TDomain, TAlgebra> T;

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
