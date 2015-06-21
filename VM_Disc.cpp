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


template<typename TDomain>
VMDisc<TDomain>::VMDisc(const char* subsets, const number init_time)
: 	IElemDisc<TDomain>("v, k, na, ca", subsets),
	m_k_out(2.5), m_na_out(140.0), m_ca_out(1.5), m_celsius(37.0),
	m_v(0), m_na(0), m_k(0), m_ca(0),
	m_ena(50.0), m_ek(-77.0), m_eca(138.0), m_eleak(-54.4),
	m_spec_res(1.0e6), m_spec_cap(1.0e-5), m_influx_ac(1e-9),
	m_aDiameter(GlobalAttachments::attachment<ANumber>("diameter")),
	m_constDiam(1e-6), m_bConstDiamSet(false),
#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
	m_spSH(SPNULL),
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
	m_spSD(SPNULL),
#endif
	m_bNonRegularGrid(false),
	m_init_time(init_time), m_ass_time(init_time-1.0),
	m_bLocked(false)
{
	// set diff constants
	m_diff.resize(3);
	m_diff[0] = 1.0e-12;
	m_diff[1] = 1.0e-12;
	m_diff[2] = 2.2e-13;
}


template<typename TDomain>
size_t VMDisc<TDomain>::_v_()
{
	return m_v_;
}

template<typename TDomain>
size_t VMDisc<TDomain>::_k_()
{
	return m_k_;
}

template<typename TDomain>
size_t VMDisc<TDomain>::_na_()
{
	return m_na_;
}

template<typename TDomain>
size_t VMDisc<TDomain>::_ca_()
{
	return m_ca_;
}



template<typename TDomain>
number VMDisc<TDomain>::ca_out()
{
        return m_ca_out;
}


template<typename TDomain>
number VMDisc<TDomain>::na_out()
{
        return m_na_out;
}


template<typename TDomain>
number VMDisc<TDomain>::k_out()
{
	return m_k_out;
}

template<typename TDomain>
void VMDisc<TDomain>::set_celsius(number cels)
{
	m_celsius = cels;
}

template<typename TDomain>
number VMDisc<TDomain>::celsius()
{
	return m_celsius;
}



template<typename TDomain>
VMDisc<TDomain>* VMDisc<TDomain>::get_VmDisc()
{
	return this;
}


template<typename TDomain>
void VMDisc<TDomain>::
set_eca(number value)
{
	m_eca = value;
}

template<typename TDomain>
number VMDisc<TDomain>::
eca()
{
	return m_eca;
}

template<typename TDomain>
void VMDisc<TDomain>::
set_ek(number value)
{
	m_ek = value;
}


template<typename TDomain>
number VMDisc<TDomain>::ek()
{
	return m_ek;
}

template<typename TDomain>
void VMDisc<TDomain>::
set_ena(number value)
{
	m_ena = value;
}

template<typename TDomain>
number VMDisc<TDomain>::ena()
{
	return m_ena;
}

template<typename TDomain>
number VMDisc<TDomain>::eleak()
{
	return m_eleak;
}

template<typename TDomain>
void VMDisc<TDomain>::
set_eleak(number value)
{
	m_eleak = value;
}

template<typename TDomain>
number VMDisc<TDomain>::
flux_ca()
{
	return m_ca;
}

template<typename TDomain>
number VMDisc<TDomain>::
flux_na()
{
	return m_na;
}

template<typename TDomain>
number VMDisc<TDomain>::
flux_k()
{
	return m_k;
}

template<typename TDomain>
number VMDisc<TDomain>::
flux_v()
{
	return m_v;
}


template<typename TDomain>
void VMDisc<TDomain>::
set_diameter(const number d)
{
	UG_COND_THROW(m_bLocked, "Diameter cannot be (re)set after VMDisc "
				  "has been added to domain discretization.");

	m_constDiam = d;
	m_bConstDiamSet = true;
}


template<typename TDomain>
void VMDisc<TDomain>::
set_spec_res(number val)
{
	m_spec_res = val;
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
void VMDisc<TDomain>::write_AllGattings_on_position(number x, number y, number z)
{

	// Vector with all needed Filename as ofstreams
	std::vector<std::vector<SmartPtr<std::ofstream> > > Vec_ofstreams;

	// Vector with all Gating Vectors of all Channels
	std::vector<number> ChannelGate;

	for (size_t i=0; i < m_channel.size(); i++)
	{
		UG_LOG("first all gating" << std::endl);
		ChannelGate = m_channel[i]->allGatingAccesors(x, y, z);
		UG_LOG("after first gating" << std::endl);
		/*std::vector<SmartPtr<std::ofstream> > vec;
		Vec_ofstreams.push_back(vec);
		// getting all States from channel i
		for (size_t j=0; j < ChannelGate[i].size(); j++)
		{
			// building char stream for gate
			std::stringstream ssoStreamName;
			ssoStreamName << "ChannelNumber_" << i << "_GateNumber_" << j << ".txt";
			std::string soStream = ssoStreamName.str();
			const char* CharStream = soStream.c_str();

			//creates ofstream for every channel Gate
			//std::ofstream NewStream(CharStream);
			SmartPtr<std::ofstream> NewStreamm;
			NewStreamm = make_sp(new std::ofstream(CharStream));

			Vec_ofstreams[i].push_back(NewStreamm);


			UG_LOG("before output" << std::endl);
			*Vec_ofstreams[i][j] << (ChannelGate[i][j]) << "/n";
			std::cout << (ChannelGate[i][j]) << std::endl;


		}*/
	}

	// testing if it is working
	/*for (size_t i=0; i < Vec_ofstreams.size(); i++)
	{
		for (size_t j=0; j < Vec_ofstreams[i].size(); j++)
		{
			std::cout << Vec_ofstreams[i][j] << std::endl;
		}
	}*/


}




template<typename TDomain>
void VMDisc<TDomain>::update_time(const number newTime, Edge* edge)
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	vrt_list vl;
	this->approx_space()->domain()->grid()->associated_elements(vl, edge);
	for (size_t vrt = 0; vrt < vl.size(); ++vrt)
		m_aaTime[vl[vrt]] = newTime;
}


template<typename TDomain>
void VMDisc<TDomain>::save_old_sol(const LocalVector& u, Edge* edge)
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	vrt_list vl;
	this->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge);
	for (size_t vrt = 0; vrt < vl.size(); ++vrt)
		for (size_t i = 0; i < m_numb_funcs+1; ++i)
			m_aaUold[vl[vrt]][i] = u(i, vrt);
}

template<typename TDomain>
void VMDisc<TDomain>::get_vm(std::vector<number>& outValues, Edge* edge) const
{
    for (size_t vrt = 0; vrt < edge->num_vertices(); ++vrt)
    	outValues.push_back(m_aaUold[edge->vertex(vrt)][m_v_]);
}


template<typename TDomain>
number VMDisc<TDomain>::get_vm(Vertex* vrt) const
{
	return m_aaUold[vrt][m_v_];
}




template<typename TDomain>
void VMDisc<TDomain>::approximation_space_changed()
{
	// only do this the first time the approx changes (when it is initially set)
	if (m_bLocked) return;


	SmartPtr<MultiGrid> grid = this->approx_space()->domain()->grid();

	// create time attachment and accessor
	if (grid->has_vertex_attachment(m_aTime))
		UG_THROW("Time attachment necessary for Vm disc "
				 "could not be created, since it already exists.");
	grid->attach_to_vertices_dv(m_aTime, m_init_time);

	m_aaTime = Grid::AttachmentAccessor<Vertex, ANumber>(*grid, m_aTime);

	// create old solution attachment and accessor
	if (grid->has_vertex_attachment(m_aUold))
		UG_THROW("Old solution attachment necessary for Vm disc "
				 "could not be created, since it already exists.");
	grid->attach_to_vertices(m_aUold);

	m_aaUold = Grid::AttachmentAccessor<Vertex, AVector4>(*grid, m_aUold);

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



template <typename TDomain>
void VMDisc<TDomain>::
prep_timestep(number time, VectorProxyBase* upb)
{
	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;

	// TODO: implement me further!
}


// ///////////////////////////////////////////////////////////
// TODO														//
// It would be preferable to do this in one loop instead of	//
// element-wise. This would also enable us to call the		//
// update_presyn() method of the synapse_handler class from	//
// here instead of from the script.							//
// ///////////////////////////////////////////////////////////
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
	MGSubsetHandler& ssh = *this->approx_space()->domain()->subset_handler();

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
				d(m_v_, co) += -m_flux_value[i];
			}
		}

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		/// if a synapse provider is available
		if	(m_spSH.valid())
		{
			number current = 0;
			//UG_LOG_ALL_PROCS("In synapse Provider!!!" << "!"<<std::endl);
			// ... and assemble to defect if synapse present
			if (m_spSH->synapse_on_edge(pElem, co, time, current))
			{
				//UG_LOG_ALL_PROCS("Setting Current" << "!"<<std::endl);
				//UG_LOG_ALL_PROCS("Current: " << current << std::endl);
				d(m_v_, co) += 1e-12*current; // scaling from nA to C/ms
			}
		}
#endif

#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
		/// if a synapse distributor is available
		if	(m_spSD.valid())
		{
			number current = 0;

			// ... and assemble to defect if synapse present
			if (m_spSD->has_active_synapses(pElem, co, time, current))
			{
				//UG_LOG_ALL_PROCS("Setting Current" << "!"<<std::endl);
				//UG_LOG_ALL_PROCS("Current: " << current << std::endl);
				d(m_v_, co) += current;
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

			// getting subset of vertex
			size_t siElem = ssh.get_subset_index(pElem);
			std::string SName = ssh.get_subset_name(siElem);

			for (size_t j = 0; j<Subsets.size(); j++)
			{
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
		d(m_v_, co) += scv.volume()*PI*Diam * allOutCurrentValues[0];
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
			VecScaleAppend(grad_c, u(m_v_,sh), scvf.global_grad(sh));

		// scalar product with normal
		number grad_normal = VecDot(grad_c, scvf.normal());

		// scale by 1/resistance and by length of element
		number diff_flux = grad_normal * element_length / (m_spec_res*pre_resistance);

		// debug
		UG_ASSERT(std::fabs(diff_flux) < 1.0,
				  "m_spec_res: " << m_spec_res << "   pre_res: " << pre_resistance << std::endl);

		// add to local defect of VM
		d(m_v_, scvf.from()) -= diff_flux;
		d(m_v_, scvf.to()  ) += diff_flux;

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
		d(m_v_, co) += PI*diam*scv.volume()*u(m_v_, co)*spec_capacity;

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
			J(m_v_, scvf.from(), m_v_, sh) -= d_diff_flux;
			J(m_v_, scvf.to()  , m_v_, sh) += d_diff_flux;

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
		//UG_COND_THROW(fabs(Diam) <= 1e-12, "Diam zero!\n");


		//spec_capa has to be set later on in an varialbe

		// get spec capacity
		number spec_capacity = m_spec_cap;
		for (size_t k = 1; k < m_numb_funcs+1; k++)
		{
			J(k, co, k, co) += scv.volume()*0.25*PI*Diam*Diam;
		}
		// potential equation
		J(m_v_, co, m_v_, co) += PI*Diam*scv.volume()*spec_capacity;
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
	return this->approx_space()->fct_id_by_name(s.c_str());
}





////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
void VMDisc<TDomain>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();

	// register prepare_timestep functionality separately, only for CPU1
	size_t aid = bridge::AlgebraTypeIDProvider::instance().id<CPUAlgebra>();
	this->set_prep_timestep_fct(aid, &VMDisc<TDomain>::prep_timestep);
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


} // namespace cable
} // namespace ug
