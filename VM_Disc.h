
/*
 * VM_Disc.cpp
 *
 *  Created on: 26.11.2014
 *      Author: Pgottmann
 */


#ifndef VM_DISC_H_
#define VM_DISC_H_

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

template <typename TDomain, typename TAlgebra>
class VMDisc
	: public IElemDisc<TDomain>
{
	public:
		VMDisc();

	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	Own type
		typedef VMDisc<TDomain, TAlgebra> this_type;

	/// Type for Gridfunction
		typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	/// Type for Channels
		typedef IChannel<TDomain, TAlgebra> TIChannel;

	public:
	///	World dimension
		static const int dim = base_type::dim;

		// Params dendrit
		number m_spec_res;
		number m_spec_cap;

		// Params for HH-Fluxes
		number m_g_K;
		number m_g_Na;
		number m_g_I;

		// reversal Pot of Sodium and Potassium
		number m_sodium;
		number m_potassium;

		// params gatting
		double m_accuracy;

		//lists for all influxes
		std::vector<MathVector<dim> > m_coords;
		std::vector<number> m_flux_value;
		std::vector<number> m_beg_flux;
		std::vector<number> m_dur_flux;

		// influx ac
		number m_influx_ac;

		//list with all channels
		std::vector<SmartPtr<TIChannel> > m_channel;

		//funcs params
		std::vector<const char*> m_funcs;
		number m_numb_funcs;

	public:
	///	Constructor
		VMDisc(SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct, const char* functions, const char* subsets)
		 : IElemDisc<TDomain>(functions, subsets), m_spGridFct(spGridFct), m_numb_funcs(0), _VM_(0), _K_(1), _Na_(2)
		   {
			m_bNonRegularGrid = false;
			register_all_funcs(m_bNonRegularGrid);
		   }

	///	Destructor
		virtual ~VMDisc() {};

	/// Functions for Standard VM and dendrit dependet settings
		// set diameter for dendrit
		void set_diameter(const number d);

		// set spec_resistance
		void set_spec_res(number val);

		// set spec capa
		void set_spec_cap(number val);

		// set influx ac
		void set_influx_ac(number influx_ac);

		// set influx params (Flux-value, koordinates, duration, beginning)
		void set_influx(number Flux, number x, number y, number z, number beg, number dur);

		//Adding function for channels
		void add_channel(SmartPtr<IChannel<TDomain, TAlgebra> > Channel, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct);

		//Add func
		void add_func(const char* func);


	// inherited from IElemDisc
	public:
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
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns if hanging nodes are needed
		//virtual bool use_hanging() const;

	protected:

		///	current regular grid flag
		bool m_bNonRegularGrid;

		SmartPtr<TDomain> m_dom;					//!< underlying domain
		SmartPtr<MultiGrid> m_mg;					//!< underlying multigrid
		SmartPtr<DoFDistribution> m_dd;				//!< underlying surface dof distribution
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;
		SmartPtr<GridFunction<TDomain, TAlgebra> > m_spGridFct;

	/// dendritic radius attachment and accessor
		ANumber m_aDiameter;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaDiameter;





	private:
		// VM is needed by every Channel
		size_t _VM_;
		size_t _Na_;
		size_t _K_;









	///	register utils
	///	\{
		void register_all_funcs(bool bHang);
		template <typename TElem, typename TFVGeom>
		void register_func();
	/// \}

};


} // namespace ug

#endif // VM_DISC_H_
