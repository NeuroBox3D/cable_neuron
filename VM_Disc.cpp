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



/*set_celsius(number cels)
		celsius = cels;*/


template<typename TDomain>
VMDisc<TDomain>* VMDisc<TDomain>::get_VmDisc()
{
	return this;
}


template<typename TDomain>
double VMDisc<TDomain>::
get_eca()
{
	return m_eca;
}

template<typename TDomain>
double VMDisc<TDomain>::
get_ena()
{
	return m_ena;
}

template<typename TDomain>
double VMDisc<TDomain>::
get_ek()
{
	return m_ek;
}

template<typename TDomain>
void VMDisc<TDomain>::
set_eca(double value)
{
	m_eca = value;
}

template<typename TDomain>
void VMDisc<TDomain>::
set_ek(double value)
{
	m_ek = value;
}

template<typename TDomain>
void VMDisc<TDomain>::
set_ena(double value)
{
	m_ena = value;
}


template<typename TDomain>
double VMDisc<TDomain>::
get_flux_ca()
{
	return m_ca;
}

template<typename TDomain>
double VMDisc<TDomain>::
get_flux_na()
{
	return m_na;
}

template<typename TDomain>
double VMDisc<TDomain>::
get_flux_k()
{
	return m_k;
}

template<typename TDomain>
double VMDisc<TDomain>::
get_flux_v()
{
	return m_v;
}


template<typename TDomain>
void VMDisc<TDomain>::
set_diameter(const number d)
{
	// handle the attachments
	if (m_spApproxSpace->domain()->grid()->has_vertex_attachment(m_aDiameter))
		UG_THROW("Radius attachment necessary for Vm disc "
				 "could not be created, since it already exists.");
	m_spApproxSpace->domain()->grid()->attach_to_vertices_dv(m_aDiameter, d);

	m_aaDiameter = Grid::AttachmentAccessor<Vertex, ANumber>(*m_spApproxSpace->domain()->grid(), m_aDiameter);
}

template <typename TDomain>
void VMDisc<TDomain>::
set_diameter_attachment(Attachment<number> diameter) {
	m_aDiameter = diameter;
}

template<typename TDomain>
void VMDisc<TDomain>::
set_spec_res(number val)
{
	m_spec_res = val;
}


template<typename TDomain>
SmartPtr<ApproximationSpace<TDomain> > VMDisc<TDomain>::approx_space()
{
	return m_spApproxSpace;
}


template<typename TDomain>
void VMDisc<TDomain>::
set_influx_ac(number influx_ac)
{
	m_influx_ac = influx_ac;
}


template<typename TDomain>
void VMDisc<TDomain>::
set_spec_cap(number val)
{
	m_spec_cap = val;
}


template<typename TDomain>
void VMDisc<TDomain>::
set_diff_coeffs(const std::vector<number>& diff_coeffs)
{
	m_diff = diff_coeffs;
}

template<typename TDomain>
void VMDisc<TDomain>::
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


#ifdef PLUGIN_SYNAPSE_PROVIDER_ENABLED
template<typename TDomain>
void VMDisc<TDomain>::
set_synapse_provider_factory(ConstSmartPtr<SynapseProviderFactory<TDomain> > spf)
{
	m_spSPF = spf;
}


template<typename TDomain>
void VMDisc<TDomain>::
set_provider_type(const std::string& providerName) {
	m_spSP = make_sp(m_spSPF->CreateProvider(providerName));
}

template <typename TDomain>
void VMDisc<TDomain>::
set_synapse_provider(SmartPtr<SynapseProvider<TDomain> sp) {
	this->m_spSP = sp;
}
#endif




template<typename TDomain>
void VMDisc<TDomain>::
add_channel(SmartPtr<IChannel<TDomain> > Channel)
{
	m_channel.push_back(Channel);

	// set this vm disc to the newly added channel
	Channel->set_vm_disc(this);
}

#if 0
template<typename TDomain>
void VMDisc<TDomain>::
add_func(std::string func)
{
	m_funcs.push_back(func);
	m_numb_funcs += 1;
}
#endif


template<typename TDomain>
void VMDisc<TDomain>::update_time(const number newTime, Edge* edge)
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	vrt_list vl;
	m_spApproxSpace->domain()->grid()->associated_elements(vl, edge);
	for (size_t vrt = 0; vrt < vl.size(); ++vrt)
		m_aaTime[vl[vrt]] = newTime;
}


template<typename TDomain>
void VMDisc<TDomain>::save_old_sol(const LocalVector& u, Edge* edge)
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	vrt_list vl;
	m_spApproxSpace->domain()->grid()->associated_elements_sorted(vl, edge);
	for (size_t vrt = 0; vrt < vl.size(); ++vrt)
		for (size_t i = 0; i < m_numb_funcs+1; ++i)
			m_aaUold[vl[vrt]][i] = u(i, vrt);
}


template<typename TDomain>
void VMDisc<TDomain>::prep_timestep_elem
(
	const number time,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	Edge* edge = dynamic_cast<Edge*>(elem);
	if (!edge) UG_THROW("VMDisc::prep_timestep_elem() called with improper element type.");

	// update old solution
	save_old_sol(u, edge);

	if (time == m_init_time)
	{
		// init channels
		for (size_t i = 0; i < m_channel.size(); ++i)
			m_channel[i]->init(u, edge);
	}
	else
	{
		// update channels
		//std::cout << "update" << std::endl;
		for (size_t i = 0; i < m_channel.size(); ++i)
			m_channel[i]->update_gating(time, u, edge);

	}

	// update time in attachments
	update_time(time, edge);
}


// Methods for Interface class

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain>::add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// get subset handler
	MGSubsetHandler& ssh = *m_spApproxSpace->domain()->subset_handler();


	// some helper vars
	number element_length2 = 0.0;
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

		element_length2 += scv.volume()/2;

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
				&& fabs(0.5*(vCornerCoords[co][0]+vCornerCoords[(co+1)%2][0]) - m_coords[i][0]) < m_influx_ac
				&& fabs(0.5*(vCornerCoords[co][1]+vCornerCoords[(co+1)%2][1]) - m_coords[i][1]) < m_influx_ac
				&& fabs(0.5*(vCornerCoords[co][2]+vCornerCoords[(co+1)%2][2]) - m_coords[i][2]) < m_influx_ac
			   )
			{
				d(_v_, co) += -m_flux_value[i];
				if (co==1)
				{
					d(_v_, 0)*=(element_length2)*PI*m_aaDiameter[pElem->vertex(0)];
					d(_v_, 1)*=(element_length2)*PI*m_aaDiameter[pElem->vertex(1)];
				}
			}
		}

#ifdef PLUGIN_SYNAPSE_PROVIDER_ENABLED
		// influxes from synapse provider factory
		if (m_spSPF.valid())
		{
			/// if a synapse provider is available
			if	(m_spSF.valid()) 
			{
				// ... and assemble to defect
				if (m_spSP->synapse_at_location(pElem, co, time, current))
				{
					//TODO add get current from synapse location and delete current from above
					d(_v_, co) += current;
				}
			}
		}
#endif

		// membrane transport mechanisms
		std::vector<number> allOutCurrentValues;
		for (size_t i = 0; i < m_numb_funcs+1; ++i)
			allOutCurrentValues.push_back(0.0);

		for (size_t i = 0; i < m_channel.size(); i++)
		{
			// if channel working on right subset
			const std::vector<std::string> Subsets = m_channel[i]->write_subsets();
			for (size_t j = 0; j<Subsets.size(); j++)
			{
				// getting subset of vertex
				size_t indexVert = ssh.get_subset_index(pElem->vertex(co));
				std::string SName = ssh.get_subset_name(indexVert);
				// if channel works on provided subset
				if (Subsets[j]==SName)
				{
					std::vector<number> outCurrentValues;

					// values we are getting from ionic_flux function in channels
					std::vector<number> vrt_values(m_numb_funcs+1);
					//for (size_t j = 0; j < m_numb_funcs+1; ++j) vrt_values[j] = u(j, co); <-- NO! this would be implicit!
					std::vector<DoFIndex> multInd;
					for (size_t j = 0; j < m_numb_funcs+1; ++j)
						vrt_values[j] = m_aaUold[pElem->vertex(co)][j];

					m_channel[i]->ionic_current(pElem->vertex(co), vrt_values, outCurrentValues);

					const std::vector<std::string>& functions = m_channel[i]->write_fcts();

					// adding defect for every ion species involved
					for (size_t k = 0; k < outCurrentValues.size(); k++)
						allOutCurrentValues[get_index(functions[k])] += (outCurrentValues[k]);
				}
			}
		}

		// writing potential defects
		d(_v_, co) += scv.volume()*PI*Diam * (allOutCurrentValues[0]-influx);
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
			VecScaleAppend(grad_c, u(_v_,sh), scvf.global_grad(sh));

		// scalar product with normal
		number grad_normal = VecDot(grad_c, scvf.normal());

		// scale by 1/resistance and by length of element
		number diff_flux = grad_normal * element_length / (m_spec_res*pre_resistance);

		if (diff_flux > 1.0 || diff_flux < -1.0)
		{
			UG_LOG("m_spec_res: " << m_spec_res << "   pre_res: " << pre_resistance << std::endl);
		}

		// add to local defect of VM
		d(_v_, scvf.from()) -= diff_flux;
		d(_v_, scvf.to()  ) += diff_flux;

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

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain>::add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
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
		d(_v_, co) += PI*diam*scv.volume()*u(_v_, co)*spec_capacity;

		// ion species time derivative
		for (size_t k = 1; k < m_numb_funcs+1; k++)
			d(k, co) += u(k, co)*scv.volume()*0.25*PI*diam*diam;
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain>::
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
			J(_v_, scvf.from(), _v_, sh) -= d_diff_flux;
			J(_v_, scvf.to()  , _v_, sh) += d_diff_flux;

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


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain>::
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
		J(_v_, co, _v_, co) += PI*Diam*scv.volume()*spec_capacity;
	}
	//std::cout << "jac m elem ends" << std::endl;
}


template<typename TDomain>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain>::fsh_elem_loop()
{


}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain>::prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try
	{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("Cannot update Finite Volume Geometry.\n");
}


template<typename TDomain>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{



}



template<typename TDomain>
void VMDisc<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
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
			_v_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
		if (m_funcs[i] == "K")
			_K_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
		if (m_funcs[i] == "Na")
			_Na_ = m_spGridFct->fct_id_by_name(m_funcs[i]);
	}
	//std::cout << "after prepare" << std::endl;
	*/
}

template<typename TDomain>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// nothing to do
}


////////////////////////////////////////////////////////////////////////////////
// Functions for ion handling
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
size_t VMDisc<TDomain>::get_index(std::string s)
{
	return m_spApproxSpace->fct_id_by_name(s.c_str());
}





////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
void VMDisc<TDomain>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef VMDisc<TDomain> T;

	this->set_prep_timestep_elem_fct(id, &T::prep_timestep_elem);
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
	template class VMDisc<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class VMDisc<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class VMDisc<Domain3d>;
#endif


} /* namespace ug */
