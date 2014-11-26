/*
 * channel_interface.h
 *
 *  Created on: 29.10.2014
 *      Author: mbreit
 */

#ifndef CHANNEL_INTERFACE_H_
#define CHANNEL_INTERFACE_H_

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




namespace ug
{

template <typename TDomain, typename TAlgebra>
class IChannel
	: public IElemDisc<TDomain>
{
	public:
		IChannel();

	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	Own type
		typedef IChannel<TDomain, TAlgebra> this_type;

	/// Type for Gridfunction
		typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	public:
	///	World dimension
		static const int dim = base_type::dim;



	public:
	///	Constructor
	// TODO: - define which variables will be influenced
	// TODO: - define which variables are needed for flux computation
		IChannel(SmartPtr<ApproximationSpace<TDomain> > approx, const char* functions, const char* subsets)
		 : IElemDisc<TDomain>(functions, subsets), m_spApproxSpace(approx)
		   {
			m_bNonRegularGrid = false;
			register_all_funcs(m_bNonRegularGrid);
		   }

	///	Destructor
		virtual ~IChannel() {};



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

	///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// functions for setting some HH params
		set_accuracy(double ac);

		set_consts(number Na, number K, number L);

		set_rev_pot(number R_Na, number R_K);



	/// initializes the defined channel type
	/** During the initialization, the necessary attachments are attached to the vertices
	 *	and their values calculated by the equilibrium state for the start membrane potential.
	**/

		virtual void init(number time, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) = 0;

	/// updates the gating parameters
		virtual void update_gating(number newTime, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) = 0;

	/// provides the ionic current (mol*s^-1) at a given vertex
		virtual void ionic_current(Vertex* v, std::vector<number>& outCurrentValues) = 0;

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
		//AdaptionSurfaceGridFunction<TDomain> m_AdaptSGF;
		//using IElemDisc<TDomain>::dim;




	private:
		// VM is needed by every Channel
		static const size_t _VM_ = 0;







	///	register utils
	///	\{
		void register_all_funcs(bool bHang);
		template <typename TElem, typename TFVGeom>
		void register_func();
	/// \}

};


template <typename TDomain, typename TAlgebra>
class ChannelHH
	: public IChannel<TDomain, TAlgebra>
{
	public:
	/// constructor
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

		IFunction<number>* m_Injection;

		ChannelHH	  (	SmartPtr<ApproximationSpace<TDomain> > approx,
						const char* functions,
						const char* subsets
					  )
		: IChannel<TDomain, TAlgebra>(approx, functions, subsets),
		  			  m_dom(approx->domain()), m_mg(approx->domain()->grid()), m_dd(approx->dof_distribution(GridLevel::TOP)),
		  			  m_bNonRegularGrid(false)
		  			  {
						register_all_funcs(m_bNonRegularGrid);
		  			  };

		/// destructor
		virtual ~ChannelHH() {};

		static const int dim = IChannel<TDomain, TAlgebra>::dim;

		// inherited from IChannel
		virtual void init(number time, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct);
		virtual void update_gating(number newTime, SmartPtr<ApproximationSpace<TDomain> > approx, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct);
		virtual void ionic_current(Vertex* v, std::vector<number>& outCurrentValues);

	///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	protected:
		SmartPtr<TDomain> m_dom;					//!< underlying domain
		SmartPtr<Grid> m_mg;						//!< underlying multigrid
		SmartPtr<DoFDistribution> m_dd;				//!< underlying surface dof distribution
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;
		bool m_bNonRegularGrid;


		///	register utils
		///	\{
			void register_all_funcs(bool bHang);
			template <typename TElem, typename TFVGeom>
			void register_func();
		/// \}

	private:
		// one attachment per state variable
		ADouble m_MGate;							//!< activating gating "particle"
		ADouble m_HGate;							//!< inactivating gating "particle"
		ADouble m_NGate;
		ADouble m_Vm;								//!< membrane voltage (in mili Volt)
		number m_rate_h;
		number m_rate_m;
		number m_rate_n;

		Grid::AttachmentAccessor<Vertex, ADouble> m_aaMGate;	//!< accessor for activating gate
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaHGate;	//!< accessor for inactivating gate
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaNGate;	//!< accessor for inactivating gate
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaVm;		//!< accessor for membrane potential

	/// Base type
		typedef IChannel<TDomain, TAlgebra> base_type;

	///	Own type
		typedef ChannelHH<TDomain, TAlgebra> this_type;

	/// GridFunction type
		typedef GridFunction<TDomain, TAlgebra> TGridFunction;
};

} // namespace ug

#endif // CHANNEL_INTERFACE_H_
