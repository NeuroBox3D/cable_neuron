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






namespace ug
{

// forward declaration
template <typename TDomain, typename TAlgebra>
class IChannel;


template <typename TDomain, typename TAlgebra>
class VMDisc
	: public IElemDisc<TDomain>
{
	public:
		/// type for channels
		typedef IChannel<TDomain, TAlgebra> TIChannel;

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

		// dendritic params
		number m_spec_res;	// mV * ms * m / C
		number m_spec_cap;	// C / (mV * m^2)

		// celsius
		number celsius;

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
		// We generally only have the following functions: VM, K, Na, Ca.
		// We should begin with generating the necessary functions in a hard-coded way
		// and might then check whether some of them are not in fact needed.
		VMDisc
		(
			//const char* functions,
			const char* subsets,
			const std::vector<SmartPtr<TIChannel> >& channels,
			SmartPtr<ApproximationSpace<TDomain> > approx
		)
		: IElemDisc<TDomain>("v, k, na, ca", subsets),
		  k_out(3.3918292968), na_out(107.796654), ca_out(1.5),
		  m_spec_res(1.0e6), m_spec_cap(1.0e-5),
		  m_influx_ac(1e-9), m_channel(channels), m_aDiameter("diameter"),
		  m_spApproxSpace(approx), m_spDD(m_spApproxSpace->dof_distribution(GridLevel::TOP)),
		  m_bNonRegularGrid(false)
		{
			// set this vm disc to each of its channels
			for (size_t i = 0; i < m_channel.size(); i++)
				m_channel[i]->set_vm_disc(this);

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

		// functions for setting params
		/// set constant diameter for dendrites
		void set_diameter(const number d);

		/// set spec_resistance
		void set_spec_res(number val);

		/// set specific capacity
		void set_spec_cap(number val);

		/// set diffusion coefficients
		void set_diff_coeffs(const std::vector<number>& diff_coeffs);

		/// set influx position accuracy
		void set_influx_ac(number influx_ac);

		/// set influx params (flux value, coordinates, beginning, duration)
		void set_influx(number Flux, number x, number y, number z, number beg, number dur);

#if 0
		/// adding a channel
		void add_channel(SmartPtr<IChannel<TDomain, TAlgebra> > Channel);

		/// add func
		void add_func(std::string func);

#endif
		SmartPtr<ApproximationSpace<TDomain> > approx_space();


		/// determines the function index based on its name
		size_t get_index(std::string s);

		void update(number dt, ConstSmartPtr<typename TAlgebra::vector_type> uOld);


	// inherited from IElemDisc
	public:

		///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

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

	protected:
		/// dendritic radius attachment and accessor
		ADouble m_aDiameter;
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaDiameter;

		/// approx space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;
		ConstSmartPtr<DoFDistribution> m_spDD;

		/// pointer to old solution	for later use in explicit assemblings
		ConstSmartPtr<typename TAlgebra::vector_type> m_uOld;

		///	current regular grid flag
		bool m_bNonRegularGrid;

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
