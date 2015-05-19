/*
 * VM_Disc.cpp
 *
 *  Created on: 26.11.2014
 *      Author: Pgottmann
 */


#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__VM_DISC_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__VM_DISC_H__

#include "lib_grid/lg_base.h"
#include "lib_grid/grid/grid_base_objects.h"

#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/local_transfer_interface.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/interpolate.h"

#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "bridge/util_domain_dependent.h"

#include <vector>
#include <stdio.h>

#include "bindings/lua/lua_user_data.h"

#include "common/util/smart_pointer.h"
#include "common/util/vector_util.h"

#include "channel_interface.h"

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
	#include "../synapse_handler/synapse_handler.h"
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
	#include "../synapse_distributor/synapse_distributor.h"
#endif

namespace ug {
namespace synapse_handler {

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
	// forward declaration
	template <typename TDomain>
	class NETISynapseHandler;
#endif
}

// forward declaration
template <typename TDomain>
class IChannel;

template <typename TDomain>
class VMDisc
	: public IElemDisc<TDomain>
{
	public:

	/// test only!
#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		synapse_handler::NETISynapseHandler<TDomain> sh;
#endif
		/// type for channels
		typedef IChannel<TDomain> TIChannel;

		///	world dimension
		static const int dim = IElemDisc<TDomain>::dim;

	public:
		// indices for functions in this elemDisc
		static const size_t _v_ = 0;
		static const size_t _k_  = 1;
		static const size_t _na_ = 2;
		static const size_t _ca_ = 3;

		// outer concentrations
		const number k_out;		// mol/m^3 = mM
		const number na_out;	// mol/m^3 = mM
		const number ca_out;	// mol/m^3 = mM

		// temperature [in degrees C]
		double celsius;

		double m_v, m_na, m_k, m_ca;

		double m_ena, m_ek, m_eca;

		// dendritic params
		number m_spec_res;	// mV * ms * m / C
		number m_spec_cap;	// C / (mV * m^2)

		// diffusion coefficients
		std::vector<number> m_diff;

		// lists for all influxes
		std::vector<MathVector<dim> > m_coords;
		std::vector<number> m_flux_value;
		std::vector<number> m_beg_flux;
		std::vector<number> m_dur_flux;

		// influx ac
		number m_influx_ac;

		// channel machanisms
		std::vector<SmartPtr<TIChannel> > m_channel;

		// number of ions
		static const size_t m_numb_funcs = 3;

	public:

		///	constructor
		// TODO: rework this!
		// We generally only have the following functions: v, k, na, ca.
		// We should begin with generating the necessary functions in a hard-coded way
		// and might then check whether some of them are not in fact needed.
		VMDisc
		(
			//const char* functions,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
			const number init_time = 0.0
		)
		: IElemDisc<TDomain>("v, k, na, ca", subsets),
		  k_out(2.5), na_out(140.0), ca_out(1.5), //celsius(37),
		  m_v(0), m_na(0), m_k(0), m_ca(0),
		  m_ena(0), m_ek(0), m_eca(0),
		  m_spec_res(1.0e6), m_spec_cap(1.0e-5), m_influx_ac(1e-9),
		  m_bDiamNotSet(true),
#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		  m_spSP(SPNULL),
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
		  m_spSD(SPNULL),
#endif
		  m_spApproxSpace(approx), m_spDD(m_spApproxSpace->dof_distribution(GridLevel::TOP)),
		  m_bNonRegularGrid(false),
		  m_init_time(init_time)
		{
			// set diff constants
			m_diff.resize(3);
			m_diff[0] = 1.0e-12;
			m_diff[1] = 1.0e-12;
			m_diff[2] = 2.2e-13;

			// create time attachment and accessor
			if (m_spApproxSpace->domain()->grid()->has_vertex_attachment(m_aTime))
				UG_THROW("Time attachment necessary for Vm disc "
						 "could not be created, since it already exists.");
			m_spApproxSpace->domain()->grid()->attach_to_vertices_dv(m_aTime, init_time);

			m_aaTime = Grid::AttachmentAccessor<Vertex, ANumber>(*m_spApproxSpace->domain()->grid(), m_aTime);

			// create old solution attachment and accessor
			if (m_spApproxSpace->domain()->grid()->has_vertex_attachment(m_aUold))
				UG_THROW("Old solution attachment necessary for Vm disc "
						 "could not be created, since it already exists.");
			m_spApproxSpace->domain()->grid()->attach_to_vertices(m_aUold);

			m_aaUold = Grid::AttachmentAccessor<Vertex, AVector4>(*m_spApproxSpace->domain()->grid(), m_aUold);

			// handle diameter attachment
			m_aDiameter = GlobalAttachments::attachment<ANumber>(std::string("diameter"));
			SmartPtr<MultiGrid> grid = m_spApproxSpace->domain()->grid();
			if (!grid->has_attachment<Vertex>(m_aDiameter))
				grid->attach_to_vertices_dv(m_aDiameter, 1e-6);
			else
				m_bDiamNotSet = false;

			m_aaDiameter = Grid::AttachmentAccessor<Vertex, ANumber>(*grid, m_aDiameter);

#if 0
			//SmartPtr<GridFunction<TDomain, TAlgebra> > spGridPtr;
			std::vector<std::string> allFct = m_spApproxSpace->names();

			// loop channel types
			for (size_t k = 0; k < m_channel.size(); k++)
			{
				size_t new_fct_pos = 0;
				// get vector of functions of the channel
				std::vector<std::string> ch_fcts = TokenizeString(m_channel[k]->m_funcs);

				// loop functions to test whether they already exist here
				for (size_t j = 0; j < ch_fcts.size(); ++j)
				{
					bool exists = false;
					for (size_t i = 0; i < allFct.size(); i++)
					{
						if (allFct[i] == ch_fcts[j])
						{
							exists = true;
							break;
						}
					}

					// add new function to approx space
					if (!exists)
					{
						// add function to space
						add_func(ch_fcts[j]);
						m_spApproxSpace->add(ch_fcts[j].c_str(), "Lagrange", 1);
						allFct.push_back(ch_fcts[j]);

						// update local functions
						std::vector<std::string> vFct = this->symb_fcts();
						vFct.push_back(ch_fcts[j]);
						this->set_functions(vFct);

						// for every new species we also need the diffusion coefficients
						std::vector<number> diff_coeffs = m_channel[k]->get_diff();

						if (new_fct_pos < diff_coeffs.size())
						{
							m_diff.push_back(diff_coeffs[new_fct_pos]);
							++new_fct_pos;
						}
						else
						{
							UG_THROW("Need to set diffusion coefficient for ion '" + ch_fcts[j] + "' with set_diff(<LuaTable>)" +
								" if you do not want diffusion set to 0.");
						}
					}
				}

				/*
				size_t position = 0;
				std::string test = m_channel[k]->m_funcs;
				size_t start = 0;
				std::string new_func;
				size_t end;
				bool exists;
				SmartPtr<TIChannel> Channel = m_channel[k];
				std::vector<std::string> allFct = m_spApproxSpace->names();
				std::cout << "in VM channel diff " << m_channel[k]->m_funcs << std::endl;
				//std::cout << "in VM channel diff" << m_channel[k]->m_diff << std::endl;
				std::vector<number> neuer;
				neuer.resize(m_channel[k]->get_diff().size());
				neuer = m_channel[k]->get_diff();
				//std::cout << neuer.size() << "geht doch!" << std::endl;

			/// Adding all needed functions if not existent
				while (test.find(",", start) != test.npos)
				{
					exists = false;
					end = test.find(",", start);
					new_func = test.substr(start, (end-start));
					//std::cout << "test: " << test << "start: "<< start << "end: " << end << std::endl;
					start = end+2;
					//std::cout << "neue funktion" << new_func << std::endl;

					// proves if function exists
					for (size_t i = 0; i < allFct.size(); i++)
						{if (allFct[i] == new_func) { exists = true;}}
					if (exists == false)
					{
							// add function to space
							//std::cout << "tryes to add" << std::endl;
							add_func(new_func.c_str());
							m_spApproxSpace->add(new_func.c_str(), "Lagrange", 1);
							//Interpolate(0.0 , spGridPtr, new_func.c_str(), 0.0);
							//std::cout << "new function added" << std::endl;


							//std::cout << "m_diff size: " << m_channel[k]->m_diff.size() << std::endl;
							// for every channel we need also the diff coeffizients but only if it not was added before
							std::cout << "pushback m_diff postion" << std::endl;
							if (position <= neuer.size())
							{
								m_diff.push_back(neuer[position]);
								position +=1;
							} else {UG_THROW("Need to set diffusion coeffizient for ion " + new_func + " with set_diff(<LuaTable>)" +
									" if you do not need diffusion set to 0.");}
					}
				}
				// adding last func
				new_func = test.substr(start);
				// proves if function exists
				for (size_t i = 0; i < allFct.size(); i++)
					{if (allFct[i] == new_func) { exists = true;}}

				if (exists == false)
				{
					// add function to space
					//std::cout << "tryes to add" << std::endl;
					add_func(new_func.c_str());
					m_spApproxSpace->add(new_func.c_str(), "Lagrange", 1);
					//Interpolate(0.0 , spGridPtr, new_func.c_str(), 0.0);
					//std::cout << "new function added" << std::endl;


					//std::cout << "pushback m_diff" << std::endl;
					if (position <= neuer.size())
					{
						m_diff.push_back(neuer[position]);
						position +=1;
					} else {UG_THROW("Need to set diffusion coeffizient for ion " + new_func + " with set_diff(<LuaTable>)" +
							" if you do not need diffusion set to 0.");}
				}
				//std::cout << "neue funktion" << new_func << std::endl;
				*/
			}
#endif
		}

		///	destructor
		virtual ~VMDisc() {};

		/// @copydoc IElemDisc::approximation_space_changed()
		virtual void approximation_space_changed()
		{
#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
			// call update function for synapse_handler
			m_spSP->update();
#endif
		}

		// functions for setting params
		/// set constant diameter for dendrites
		void set_diameter(const number d);

		// This is done in the constructor.
		//void set_diameterGeo();

		//void set_diameter_attachment(ANumber diameter);

		/// set spec_resistance
		void set_spec_res(number val);

		/// set specific capacity
		void set_spec_cap(number val);

		/// set diffusion coefficients
		void set_diff_coeffs(const std::vector<number>& diff_coeffs);

		/// set influx position accuracy
		void set_influx_ac(number influx_ac);

		/// setting write temperature
		void set_celsius(number cels);

		/// set influx params (flux value, coordinates, beginning, duration)
		void set_influx(number Flux, number x, number y, number z, number beg, number dur);

		VMDisc<TDomain>* get_VmDisc();

		/// functions to get different ion fluxes
		double get_flux_ca();
		double get_flux_v();
		double get_flux_k();
		double get_flux_na();


		/// functions for different reversal potentials
		double get_eca();
		double get_ena();
		double get_ek();
		void set_eca(double value);
		void set_ena(double value);
		void set_ek(double value);


#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		void set_synapse_handler(synapse_handler::NETISynapseHandler<TDomain>* sp);
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
		void set_synapse_distributor(SmartPtr<SynapseDistributor> sd);
#endif
		/// adding a channel
		void add_channel(SmartPtr<IChannel<TDomain> > Channel);
#if 0
		/// add func
		void add_func(std::string func);
#endif
		SmartPtr<ApproximationSpace<TDomain> > approx_space();

	private:
		/// determines the function index based on its name
		size_t get_index(std::string s);

		/// update time in time attachments of an edge
		void update_time(const number newTime, Edge* edge);

		/// save old solution to attachments
		void save_old_sol(const LocalVector& u, Edge* edge);

	public:
		/// get vm
		void get_vm(std::vector<number>& outValues, Edge* edge) const;
		number get_vm(Vertex* vrt) const;

		/// get approx space
		ConstSmartPtr<ApproximationSpace<TDomain> > get_approximation_space() const;


	// inherited from IElemDisc
	public:

		///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		/**
		 * @brief Prepares for time step assemblings
		 *
		 * This method will be called before each time step assembling process.
		 * It can be used if any time-dependent modifications have to be made to an element
		 * before the defects/Jacobians for the specific time step can be calculated.
		 * This is especially useful if grid attachments have to be updated
		 * before each time step.
		 *
		 * @param time	new point in time
		 * @param u		local vector of unknowns
		 * @param elem	the element modifications can be made for
		 */
		virtual void prep_timestep_elem
		(
			const number time,
			const LocalVector& u,
			GridObject* elem,
			const MathVector<dim> vCornerCoords[]
		);

		///	prepares the loop over all elements
		/**
		 * This method prepares the loop over all elements. It resizes the Position
		 * array for the corner coordinates and schedules the local ip positions
		 * at the data imports.
		 */
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	prepares the element for assembling
		/**
		 * This methods prepares an element for the assembling. The Positions of
		 * the Element Corners are read and the Finite Volume Geometry is updated.
		 * The global ip positions are scheduled at the data imports.
		 */
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[]);

		///	finishes the loop over all elements
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

		///	assembles stiffness part of local defekt
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles mass part of local defect
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		/// assembles jacobian of stiffnes part
		template<typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		/// assembles jacobian of mass part
		template<typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	public:
		/// attachment and accessor for current time in each vertex (needs to be accessible by IChannel)
		ANumber m_aTime;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaTime;

	protected:
		/// dendritic radius attachment and accessor
		ANumber m_aDiameter;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaDiameter;

		bool m_bDiamNotSet;

		/// attachment and accessor for old solution in each vertex
		AVector4 m_aUold;
		Grid::AttachmentAccessor<Vertex, AVector4> m_aaUold;

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		SmartPtr<synapse_handler::NETISynapseHandler<TDomain> > m_spSP;
#endif

#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
		SmartPtr<SynapseDistributor> m_spSD;
#endif

		/// approx space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;
		ConstSmartPtr<DoFDistribution> m_spDD;

		///	current regular grid flag
		bool m_bNonRegularGrid;

		/// init time
		number m_init_time;

	private:
		///	register utils
		///	\{
		void register_all_funcs(bool bHang);

		template <typename TElem, typename TFVGeom>
		void register_func();
		/// \}
};


} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__VM_DISC_H__
