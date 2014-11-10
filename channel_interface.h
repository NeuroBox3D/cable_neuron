/*
 * channel_interface.h
 *
 *  Created on: 29.10.2014
 *      Author: mbreit
 */

#ifndef CHANNEL_INTERFACE_H_
#define CHANNEL_INTERFACE_H_

#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include <vector>
#include <stdio.h>
#include "bindings/lua/lua_user_data.h"
#include "common/util/smart_pointer.h"

namespace ug
{

template <typename TDomain>
class IChannel
	: public IElemDisc<TDomain>
{
	public:
		IChannel();
		virtual ~IChannel();

	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	Own type
		typedef IChannel<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor
	// TODO: - define which variables will be influenced
	// TODO: - define which variables are needed for flux computation
		IChannel(const char* functions, const char* subsets, SmartPtr<ApproximationSpace<TDomain> > approx)
		 : IElemDisc<TDomain>(functions,subsets)
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
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// initializes the defined channel type
	/** During the initialization, the necessary attachments are attached to the vertices
	 *	and their values calculated by the equilibrium state for the start membrane potential.
	**/
		// TODO: somehow pass arbitrary other unknowns
		virtual void init(number time) = 0;

	/// updates the gating parameters
	// TODO: somehow pass arbitrary other unknowns
		virtual void update_gating(number newTime) = 0;

	/// provides the ionic current (mol*s^-1) at a given vertex
		// TODO: somehow pass arbitrary other unknowns
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
		//using IElemDisc<TDomain>::dim;

	///	register utils
	///	\{
		void register_all_funcs(bool bHang);
		template <typename TElem, typename TFVGeom>
		void register_func();
	/// \}

};


template <typename TDomain>
class ChannelHH
	: public IChannel<TDomain>
{
	public:
	/// constructor
		///	World dimension

		ChannelHH	  (	const char* functions,
						const char* subsets,
						SmartPtr<ApproximationSpace<TDomain> > approx
					  )
		: IChannel<TDomain>(functions, subsets, approx),
		  			  m_dom(approx->domain()), m_mg(approx->domain()->grid()), m_dd(approx->dof_distribution(GridLevel::TOP)),
		  			  m_bNonRegularGrid(false), m_spApproxSpace(approx)
		  			  {
						register_all_funcs(m_bNonRegularGrid);
		  			  };

		/// destructor
		virtual ~ChannelHH() {};

		static const int dim = IChannel<TDomain>::dim;

		// inherited from IChannel
		virtual void init(number time);
		virtual void update_gating(number newTime);
		virtual void ionic_current(Vertex* v, std::vector<number>& outCurrentValues);

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

		Grid::AttachmentAccessor<Vertex, ADouble> m_aaMGate;	//!< accessor for activating gate
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaHGate;	//!< accessor for inactivating gate
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaNGate;	//!< accessor for inactivating gate
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaVm;		//!< accessor for membrane potential


		typedef IChannel<TDomain> base_type;

	///	Own type
		typedef ChannelHH<TDomain> this_type;
};

} // namespace ug

#endif // CHANNEL_INTERFACE_H_
