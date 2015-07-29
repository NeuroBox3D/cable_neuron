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

#include "lib_grid/global_attachments.h"
#include "../neuronal_topology_importer/neuronal_topology_importer.h"


namespace ug {
namespace cable {


// TODO: rework this?
// We generally only have the following functions: v, k, na, ca.
// We should begin with generating the necessary functions in a hard-coded way
// and might then check whether some of them are not in fact needed.
template<typename TDomain>
VMDisc<TDomain>::VMDisc(const char* subsets, const number init_time)
: 	IElemDisc<TDomain>("v, k, na, ca", subsets),
	R(8.314), F(96485.0),
	m_aDiameter(GlobalAttachments::attachment<ANumber>("diameter")),
	m_constDiam(1e-6), m_bConstDiamSet(false),
	m_spec_res(1.0e6), m_spec_cap(1.0e-5),
	m_k_out(4.0), m_na_out(150.0), m_ca_out(1.5),
	m_ek(-90.0), m_ena(60.0), m_eca(140.0),
	m_temperature(310.0),
	m_influx_ac(1e-9),
	m_output(false),
	m_gating_x(0), m_gating_y(0), m_gating_z(0),
	m_gating_pfad(""),
#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
	m_spSH(SPNULL),
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
	m_spSD(SPNULL),
#endif
	m_v(0), m_na(0), m_k(0), m_ca(0),
	m_init_time(init_time), m_time(init_time),
	m_bNonRegularGrid(false),
	m_bLocked(false)
{
	// set diff constants
	m_diff.resize(3);
	m_diff[0] = 1.0e-12;
	m_diff[1] = 1.0e-12;
	m_diff[2] = 2.2e-13;
}


// /////////////////////
// setting parameters //
// /////////////////////

template<typename TDomain>
void VMDisc<TDomain>::
set_diameter(const number d)
{
	UG_COND_THROW(m_bLocked, "Diameter cannot be (re)set after VMDisc "
				  "has been added to domain discretization.");

	m_constDiam = d;
	m_bConstDiamSet = true;
}

template<typename TDomain> void VMDisc<TDomain>::set_spec_res(number val) {m_spec_res = val;}
template<typename TDomain> void VMDisc<TDomain>::set_spec_cap(number val) {	m_spec_cap = val;}

template<typename TDomain> void VMDisc<TDomain>::set_k_out(number value) {m_k_out = value;}
template<typename TDomain> void VMDisc<TDomain>::set_na_out(number value) {m_na_out = value;}
template<typename TDomain> void VMDisc<TDomain>::set_ca_out(number value) {m_ca_out = value;}

template<typename TDomain> void VMDisc<TDomain>::set_diff_coeffs
	(const std::vector<number>& diff_coeffs) {m_diff = diff_coeffs;}

template<typename TDomain> void VMDisc<TDomain>::set_ek(number value) {m_ek = value;}
template<typename TDomain> void VMDisc<TDomain>::set_ena(number value) {m_ena = value;}
template<typename TDomain> void VMDisc<TDomain>::set_eca(number value) {m_eca = value;}

template<typename TDomain> void VMDisc<TDomain>::set_temperature(number kelvin) {m_temperature = kelvin;}
template<typename TDomain> void VMDisc<TDomain>::set_temperature_celsius
	(number cels){m_temperature = cels + 273.15;}


// /////////////////////
// getting parameters //
// /////////////////////

template<typename TDomain> number VMDisc<TDomain>::diameter() {return m_constDiam;}

template<typename TDomain> number VMDisc<TDomain>::spec_res() {return m_spec_res;}
template<typename TDomain> number VMDisc<TDomain>::spec_cap() {return m_spec_cap;}

template<typename TDomain> number VMDisc<TDomain>::k_out() {return m_k_out;}
template<typename TDomain> number VMDisc<TDomain>::na_out() {return m_na_out;}
template<typename TDomain> number VMDisc<TDomain>::ca_out() {return m_ca_out;}

template<typename TDomain> const std::vector<number>& VMDisc<TDomain>::diff_coeffs() {return m_diff;}

template<typename TDomain> number VMDisc<TDomain>::ek() {return m_ek;}
template<typename TDomain> number VMDisc<TDomain>::ena() {return m_ena;}
template<typename TDomain> number VMDisc<TDomain>::eca() {return m_eca;}

template<typename TDomain> number VMDisc<TDomain>::temperature() {return m_temperature;}
template<typename TDomain> number VMDisc<TDomain>::temperature_celsius() {return m_temperature - 273.15;}


// ////////////////////////////
// setters for functionality //
// ////////////////////////////

template<typename TDomain> void VMDisc<TDomain>::set_influx_ac(number influx_ac) {m_influx_ac = influx_ac;}

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

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
template <typename TDomain>
void VMDisc<TDomain>::
set_synapse_handler(SmartPtr<synapse_handler::NETISynapseHandler<TDomain> > sh)
{
	m_spSH = sh;
}
#endif


#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
template <typename TDomain>
void VMDisc<TDomain>::
set_synapse_distributor(SmartPtr<SynapseDistributor> sd)
{
	m_spSD = sd;
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


template <typename TDomain>
void VMDisc<TDomain>::set_output(bool output, number gating_x, number gating_y, number gating_z, std::string gating_pfad)
{
	m_output = output;
	m_gating_x = gating_x;
	m_gating_y = gating_y;
	m_gating_z = gating_z;
	m_gating_pfad = gating_pfad;
}


// ////////////////////////////////
// getters for functional values //
// ////////////////////////////////

/*template<typename TDomain> number VMDisc<TDomain>::flux_k() {return m_k;}
template<typename TDomain> number VMDisc<TDomain>::flux_na() {return m_na;}
template<typename TDomain> number VMDisc<TDomain>::flux_ca() {return m_ca;}
template<typename TDomain> number VMDisc<TDomain>::flux_v() {return m_v;}*/

template <typename TDomain> number VMDisc<TDomain>::time() {return m_time;}

template<typename TDomain>
number VMDisc<TDomain>::get_vm(Vertex* vrt) const
{
	ConstSmartPtr<DoFDistribution> dd = this->approx_space()->dof_distribution(GridLevel(), false);
	std::vector<DoFIndex> dofIndex;
	dd->dof_indices(vrt, _v_, dofIndex, false, false);
	UG_COND_THROW(dofIndex.size() != 1, "Not exactly one DoF index found for vertex.");

	return DoFRef(*m_spUOld, dofIndex[0]);
}


template<typename TDomain>
void VMDisc<TDomain>::write_gatings_for_position(number x, number y, number z, std::string pfad)
{
	// Vector with all needed Filename as ofstreams
	std::vector<std::vector<SmartPtr<std::ofstream> > > Vec_ofstreams;

	// Vector with all Gating Vectors of all Channels
	std::vector<std::vector<number> > ChannelGate;



	for (size_t i=0; i < m_channel.size(); i++)
	{
		// temp vector for initialisation of every channel
		std::vector<number> temp;
		ChannelGate.push_back(temp);
		std::vector<SmartPtr<std::ofstream> > vec;
		Vec_ofstreams.push_back(vec);

		// writing all Accesors of one Channel
		ChannelGate[i] = m_channel[i]->state_values(x, y, z);
		// getting all States from channel i
		for (size_t j=0; j < ChannelGate[i].size(); j++)
		{
			// building char stream for gate
			std::stringstream ssoStreamName;
			ssoStreamName << pfad << "ChannelNumber_" << i << "_GateNumber_" << j << ".txt";


			std::string soStream = ssoStreamName.str();
			const char* CharStream = soStream.c_str();

			//creates SmartPtr ofstream for every channel Gate
			SmartPtr<std::ofstream> NewStreamm;
			try { NewStreamm = make_sp(new std::ofstream(CharStream, std::ios::app)); }
			UG_CATCH_THROW("Can't create ofstream for State output. Perhaps Pfad is missing");

			Vec_ofstreams[i].push_back(NewStreamm);

			*Vec_ofstreams[i][j] << this->time() <<" \t";
			*Vec_ofstreams[i][j] << (ChannelGate[i][j]) << "\n";
			//std::cout << (ChannelGate[i][j]) << std::endl;

		}
	}

}


template <typename TDomain>
template <typename TVector>
number VMDisc<TDomain>::
estimate_cfl_cond(ConstSmartPtr<TVector> u)
{
	ConstSmartPtr<DoFDistribution> dd = this->approx_space()->dof_distribution(GridLevel(), false);
	std::vector<DoFIndex> dofIndex;
	MGSubsetHandler& ssh = *this->approx_space()->domain()->subset_handler();

	std::vector<number> vrt_values(m_numb_ion_funcs+1);

	// iterate over surface level
	number maxLinDep = 0.0;
	typedef DoFDistribution::traits<Vertex>::const_iterator it_type;
	it_type it = dd->begin<Vertex>();
	it_type it_end = dd->end<Vertex>();
	for (; it != it_end; ++it)
	{
		Vertex* vrt = *it;

		// fill vector with solution at vertex
		for (size_t j = 0; j < m_numb_ion_funcs+1; ++j)
		{
			dd->dof_indices(vrt, j, dofIndex, false, true);
			UG_COND_THROW(dofIndex.size() != 1, "Not exactly one DoF index found for vertex.");
			vrt_values[j] = DoFRef(*m_spUOld, dofIndex[0]);
		}

		// iterate over channels
		number linDep = 0.0;
		for (size_t ch = 0; ch < m_channel.size(); ++ch)
		{
			// check that vertex belongs to an edge of a subset that this channel works on
			const std::vector<std::string>& ch_subsets = m_channel[ch]->write_subsets();

			typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
			edge_list el;
			this->approx_space()->domain()->grid()->associated_elements(el, vrt);
			for (size_t k = 0; k < el.size(); ++k)
			{
				Edge* edge = el[k];
				size_t siEdge = ssh.get_subset_index(edge);
				std::string sName = ssh.get_subset_name(siEdge);

				for (size_t l = 0; l < ch_subsets.size(); ++l)
				{
					if (sName == ch_subsets[l])
						goto compute_linDep;
				}
			}
			continue;

		compute_linDep:
			linDep += m_channel[ch]->lin_dep_on_pot(vrt, vrt_values);
		}

		maxLinDep = std::max(maxLinDep, linDep);
	}


	double cfl = 2.0 * m_spec_cap / maxLinDep;
	// communicate
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		com.allreduce(&cfl, &cfl, 1, PCL_DT_DOUBLE, PCL_RO_MIN);
	}
#endif
	return (number) cfl;
}


// ///////////////////////////
// inherited from IElemDisc //
// ///////////////////////////

template<typename TDomain>
void VMDisc<TDomain>::approximation_space_changed()
{
	// only do this the first time the approx changes (when it is initially set)
	if (m_bLocked) return;

	SmartPtr<MultiGrid> grid = this->approx_space()->domain()->grid();

	// handle diameter attachment
	if (!grid->has_attachment<Vertex>(m_aDiameter))
		grid->attach_to_vertices_dv(m_aDiameter, m_constDiam);
	else
	{
		if (m_bConstDiamSet)
		{
			UG_LOG("HINT: Even though you have explicitly set a constant diameter to the domain\n"
				   "      this discretization will use the diameter information attached to the grid\n"
				   "      you specified.\n");

		}
	}

	// this will distribute the attachment values to the whole grid
	m_dah.set_attachment(m_aDiameter);
	m_dah.set_grid(grid);

	// create accessor
	m_aaDiameter = Grid::AttachmentAccessor<Vertex, ANumber>(*grid, m_aDiameter);

	// call channel init functions
	for (size_t i = 0; i < m_channel.size(); i++)
		m_channel[i]->vm_disc_available();


#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
	// call init method for synapse handler
	if (m_spSH.valid())
		m_spSH->grid_first_available();
#endif

	// lock discretization
	m_bLocked = true;
}


template<typename TDomain>
void VMDisc<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check number
	if (vLfeID.size() != 4)
		UG_THROW("VMDisc: Wrong number of functions given. Needs exactly 4 functions ");

	if (vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("VMDisc FV scheme only implemented for 1st order.");

	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_funcs(m_bNonRegularGrid);
}


template <typename TDomain>
void VMDisc<TDomain>::
prep_timestep(number time, VectorProxyBase* upb)
{

	//Write Gattings out
	if (m_output==true)
		write_gatings_for_position(m_gating_x, m_gating_y, m_gating_z, m_gating_pfad);


	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;

	m_spUOld = u.clone();

	ConstSmartPtr<DoFDistribution> dd = this->approx_space()->dof_distribution(GridLevel(), false);
	std::vector<DoFIndex> dofIndex;
	MGSubsetHandler& ssh = *this->approx_space()->domain()->subset_handler();

	std::vector<number> vrt_values(m_numb_ion_funcs+1);

	// iterate over surface level
	typedef DoFDistribution::traits<Vertex>::const_iterator it_type;
	it_type it = dd->begin<Vertex>();
	it_type it_end = dd->end<Vertex>();
	for (; it != it_end; ++it)
	{
		Vertex* vrt = *it;

		// fill vector with solution at vertex
		for (size_t j = 0; j < m_numb_ion_funcs+1; ++j)
		{
			dd->dof_indices(vrt, j, dofIndex, false, true);
			UG_COND_THROW(dofIndex.size() != 1, "Not exactly one DoF index found for vertex.");
			vrt_values[j] = DoFRef(*m_spUOld, dofIndex[0]);
		}

		// iterate over channels
		for (size_t ch = 0; ch < m_channel.size(); ++ch)
		{
			// check that vertex belongs to an edge of a subset that this channel works on
			const std::vector<std::string>& ch_subsets = m_channel[ch]->write_subsets();

			typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
			edge_list el;
			this->approx_space()->domain()->grid()->associated_elements(el, vrt);
			for (size_t k = 0; k < el.size(); ++k)
			{
				Edge* edge = el[k];
				size_t siEdge = ssh.get_subset_index(edge);
				std::string sName = ssh.get_subset_name(siEdge);

				for (size_t l = 0; l < ch_subsets.size(); ++l)
				{
					if (sName == ch_subsets[l])
						goto update;
				}
			}
			continue;


		update:
			if (time == m_init_time)
			{
				// init channel
				m_channel[ch]->init(vrt, vrt_values);
			}
			else
			{
				// update channel
				m_channel[ch]->update_gating(time, vrt, vrt_values);
			}
		}
	}

	// update time in attachments (must be done AFTER update_gating of channels)
	m_time = time;

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
	// call update_presyn() method for synapse handler
	if (m_spSH.valid())
		m_spSH->update_presyn();
#endif
}


template<typename TDomain>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	// nothing to do
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
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain>::add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type (in order to allow access to attachments)
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	// calculate some helper variables for cable equation
	number element_length = 0.0;
	number pre_resistance = 0.0;
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// get diam from attachment
		number diam = m_aaDiameter[pElem->vertex(co)];

		// add length of scv to element length
		element_length += scv.volume();

		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*diam*diam);
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

		// debug
		UG_ASSERT(std::fabs(diff_flux) < 1.0,
				  "m_spec_res: " << m_spec_res << "   pre_res: " << pre_resistance << std::endl);

		// add to local defect of VM
		d(_v_, scvf.from()) -= diff_flux;
		d(_v_, scvf.to()  ) += diff_flux;

		// diameter of axial flux cross-section
		number diam_fromTo = std::min(m_aaDiameter[pElem->vertex(scvf.from())],
							   m_aaDiameter[pElem->vertex(scvf.to())]);

		for (size_t k = 1; k < m_numb_ion_funcs+1; k++)
		{
			// compute gradient at ip
			VecSet(grad_c, 0.0);
			for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(grad_c, u(k,sh), scvf.global_grad(sh));

			// scalar product with normal
			grad_normal = VecDot(grad_c, scvf.normal());

			// scale by cross section and diff const
			diff_flux = grad_normal * m_diff[k-1] * 0.25*PI * diam_fromTo*diam_fromTo;

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
		for (size_t k = 1; k < m_numb_ion_funcs+1; k++)
			d(k, co) += u(k, co)*scv.volume()*0.25*PI*diam*diam;
	}
}


template<typename TDomain>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// dof distribution and DoFIndex vector (for accessing old solution)
	ConstSmartPtr<DoFDistribution> dd = this->approx_space()->dof_distribution(GridLevel(), false);
	std::vector<DoFIndex> dofIndex;

	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// get subset handler
	MGSubsetHandler& ssh = *this->approx_space()->domain()->subset_handler();

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
		number diam = m_aaDiameter[pElem->vertex(co)];


		// influx handling
		number time = this->time();
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
				// use real current here, thus the influx is independent from the geometry
				d(_v_, co) += m_flux_value[i];
			}
		}

		// synapses handled by synapse_handler
#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		if	(m_spSH.valid())
		{
			number current = 0;
			//UG_LOG_ALL_PROCS("In synapse Provider!!!" << "!"<<std::endl);
			// ... and assemble to defect if synapse present
			if (m_spSH->synapse_on_edge(pElem, co, time, current))
			{
				//UG_LOG_ALL_PROCS("Setting Current" << "!"<<std::endl);
				//UG_LOG_ALL_PROCS("Current: " << current << std::endl);
				d(_v_, co) -= 1e-12*current; // scaling from nA to C/ms
			}
		}
#endif

		// synapses handled by synapse_distributor (to be removed)
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
		if	(m_spSD.valid())
		{
			number current = 0;

			// ... and assemble to defect if synapse present
			if (m_spSD->has_active_synapses(pElem, co, time, current))
			{
				//UG_LOG_ALL_PROCS("Setting Current" << "!"<<std::endl);
				//UG_LOG_ALL_PROCS("Current: " << current << std::endl);
				d(_v_, co) -= current;
			}
		}
#endif

		// membrane transport mechanisms (IChannels)
		std::vector<number> allOutCurrentValues;
		for (size_t i = 0; i < m_numb_ion_funcs+1; ++i)
			allOutCurrentValues.push_back(0.0);

		for (size_t i = 0; i < m_channel.size(); i++)
		{
			// if channel working on right subset
			const std::vector<std::string> Subsets = m_channel[i]->write_subsets();

			// getting subset of vertex
			size_t siElem = ssh.get_subset_index(pElem);
			std::string sName = ssh.get_subset_name(siElem);

			for (size_t j = 0; j<Subsets.size(); j++)
			{
				// if channel works on provided subset
				if (Subsets[j]==sName)
				{
					std::vector<number> outCurrentValues;

					// values we are getting from ionic_flux function in channels
					std::vector<number> vrt_values(m_numb_ion_funcs+1);
					for (size_t j = 0; j < m_numb_ion_funcs+1; ++j)
					{
						dd->dof_indices(pElem->vertex(co), j, dofIndex, false, true);
						UG_COND_THROW(dofIndex.size() != 1, "Not exactly one DoF index found for vertex.");
						vrt_values[j] = DoFRef(*m_spUOld, dofIndex[0]);
					}
					m_channel[i]->ionic_current(pElem->vertex(co), vrt_values, outCurrentValues);

					const std::vector<std::string>& functions = m_channel[i]->write_fcts();

					// adding defect for every ion species involved
					for (size_t k = 0; k < outCurrentValues.size(); k++)
					{
						size_t fct_index = this->approx_space()->fct_id_by_name(functions[k].c_str());
						allOutCurrentValues[fct_index] += (outCurrentValues[k]);
					}
				}
			}
		}

		// writing potential defects
		d(_v_, co) -= scv.volume()*PI*diam * allOutCurrentValues[0];

		// writing ion species defects
		for (size_t k = 1; k < m_numb_ion_funcs+1; k++)
			d(k, co) -= scv.volume()*PI*diam * allOutCurrentValues[k];
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
		number diam = m_aaDiameter[pElem->vertex(co)];

		// add length of scv to element length
		element_length += scv.volume();

		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*diam*diam);
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

			for (size_t k = 1; k < m_numb_ion_funcs+1; k++)
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
		number diam = m_aaDiameter[pElem->vertex(co)];
		//UG_COND_THROW(fabs(Diam) <= 1e-12, "Diam zero!\n");


		//spec_capa has to be set later on in an varialbe

		// get spec capacity
		number spec_capacity = m_spec_cap;
		for (size_t k = 1; k < m_numb_ion_funcs+1; k++)
		{
			J(k, co, k, co) += scv.volume()*0.25*PI*diam*diam;
		}
		// potential equation
		J(_v_, co, _v_, co) += PI*diam*scv.volume()*spec_capacity;
	}
	//std::cout << "jac m elem ends" << std::endl;
}


template<typename TDomain>
template <typename TElem, typename TFVGeom>
void VMDisc<TDomain>::fsh_elem_loop()
{
	// nothing to do
}


// ///////////////////////////////
//	register assemble functions //
// ///////////////////////////////

template<typename TDomain>
void VMDisc<TDomain>::
register_all_funcs(bool bHang)
{
	// register prepare_timestep functionality separately, only for CPU1
#ifdef UG_CPU_1
	size_t aid = bridge::AlgebraTypeIDProvider::instance().id<CPUAlgebra>();
	this->set_prep_timestep_fct(aid, &VMDisc<TDomain>::prep_timestep);
#else
	UG_THROW("CPUAlgebra type not present. Please make sure that your UG4 compile options include it.");
#endif

	// register assembling functionality
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void VMDisc<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef VMDisc<TDomain> T;

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
}



// ////////////////////////////////////
//	explicit template instantiations //
// ////////////////////////////////////

#ifdef UG_DIM_1
	template class VMDisc<Domain1d>;
	#ifdef UG_CPU_1
		template number VMDisc<Domain1d>::estimate_cfl_cond(ConstSmartPtr<CPUAlgebra::vector_type> u);
	#endif
#endif

#ifdef UG_DIM_2
	template class VMDisc<Domain2d>;
	#ifdef UG_CPU_1
		template number VMDisc<Domain2d>::estimate_cfl_cond(ConstSmartPtr<CPUAlgebra::vector_type> u);
	#endif
#endif

#ifdef UG_DIM_3
	template class VMDisc<Domain3d>;
	#ifdef UG_CPU_1
		template number VMDisc<Domain3d>::estimate_cfl_cond(ConstSmartPtr<CPUAlgebra::vector_type> u);
	#endif
#endif


} // namespace cable
} // namespace ug
