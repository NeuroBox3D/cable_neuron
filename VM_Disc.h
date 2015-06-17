/*
 * VM_Disc.cpp
 *
 *  Created on: 26.11.2014
 *      Author: Pgottmann
 */


#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__VM_DISC_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__VM_DISC_H__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "bridge/util_algebra_dependent.h"

#include "diam_attachment_handler.h"	// attachment handling for diameter attachment
#include "channel_interface.h"

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
	#include "../synapse_handler/synapse_handler.h"
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
	#include "../synapse_distributor/synapse_distributor.h"
#endif



namespace ug {

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
namespace synapse_handler {
	// forward declaration
	template <typename TDomain>
	class NETISynapseHandler;
}
#endif


namespace cable {


// forward declaration
template <typename TDomain>
class IChannel;

template <typename TDomain>
class VMDisc
	: public IElemDisc<TDomain>
{
	public:
		/// type for channels
		typedef IChannel<TDomain> TIChannel;

		///	world dimension
		static const int dim = IElemDisc<TDomain>::dim;

	
	public:

		///	constructor
		// TODO: rework this?
		// We generally only have the following functions: v, k, na, ca.
		// We should begin with generating the necessary functions in a hard-coded way
		// and might then check whether some of them are not in fact needed.
		VMDisc(const char* subsets, const number init_time = 0.0);

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

		/// setting write temperature
		void set_celsius(number cels);

		number celsius();

		/// set influx params (flux value, coordinates, beginning, duration)
		void set_influx(number Flux, number x, number y, number z, number beg, number dur);

		VMDisc<TDomain>* get_VmDisc();

		/// functions to get different ion fluxes
		number flux_ca();
		number flux_v();
		number flux_k();
		number flux_na();


		/// functions for reversal potentials
		number eca();
		number ena();
		number ek();
		number eleak();

		void set_ek(number value);
		void set_ena(number value);
		void set_eca(number value);
		void set_eleak(number value);

		/// functions for outer concentrations
		number k_out();
		number na_out();
		number ca_out();

		/// functions for Function size_t indexes
		size_t _v_();
		size_t _k_();
		size_t _na_();
		size_t _ca_();


#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		void set_synapse_handler(SmartPtr<synapse_handler::NETISynapseHandler<TDomain> > sh);
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

	private:

	    // outer concentrations
	    const number m_k_out;             // mol/m^3 = mM
	    const number m_na_out;    // mol/m^3 = mM
	    const number m_ca_out;    // mol/m^3 = mM

	    // temperature [in degrees C]
	    number m_celsius;

	    number m_v, m_na, m_k, m_ca;

	    // reversible potentials
	    number m_ena, m_ek, m_eca;

	    // leakeage Term
	    number m_eleak;

	    /// settings for dendrite
	    number m_spec_res, m_spec_cap;

	    number m_influx_ac;


	/// values for influxes
	std::vector<number> m_flux_value, m_beg_flux, m_dur_flux;

	/// vector for influx coordinates x, y, z
	std::vector<MathVector<dim> > m_coords;

	/// vector for diffusion consts of k, na and ca
	std::vector<number> m_diff;

    /// List of Channels
	std::vector<SmartPtr<TIChannel> > m_channel;

    ///     world dimension
    static const size_t m_v_ = 0;
    static const size_t m_k_ = 1;
    static const size_t m_na_ = 2;
    static const size_t m_ca_ = 3;

    //
    static const size_t m_numb_funcs = 3;

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

	// inherited from IElemDisc
	public:
		/// @copydoc IElemDisc::approximation_space_changed()
		virtual void approximation_space_changed();

		///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		/// prepare the timestep
		virtual void prep_timestep(number time, VectorProxyBase* up);

		/// prepares elements for time step assembling
		virtual void prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	prepares the loop over all elements
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	prepares the element for assembling
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[]);

		///	finishes the loop over all elements
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

		///	assembles stiffness part of local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles mass part of local defect
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		/// assembles jacobian of stiffness part
		template<typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		/// assembles jacobian of mass part
		template<typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	public:
		// TODO: Is this really necessary?
		// As far as I can tell, in every channel that we have, the time is only used to compute a time step dt
		// which we could easily get by passing dt directly to update_gating() instead of newTime.
		/// attachment and accessor for current time in each vertex (needs to be accessible by IChannel)
		ANumber m_aTime;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaTime;

	protected:
		/// dendritic radius attachment and accessor
		ANumber m_aDiameter;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaDiameter;

		/// handler for multigrid usage of attachment
		DiamAttachmentHandler m_dah;

		number m_constDiam;
		bool m_bConstDiamSet;

		/// attachment and accessor for old solution in each vertex
		AVector4 m_aUold;
		Grid::AttachmentAccessor<Vertex, AVector4> m_aaUold;

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		/// synapse handler
		SmartPtr<synapse_handler::NETISynapseHandler<TDomain> > m_spSH;
#endif

#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
		/// synapse distributor
		SmartPtr<SynapseDistributor> m_spSD;
#endif

		///	current regular grid flag
		bool m_bNonRegularGrid;

		/// init time
		number m_init_time;

		/// assembling time
		number m_ass_time;

		/// flag indicating whether approx space has been set
		bool m_bLocked;

	private:



		///	register utils
		///	\{
		void register_all_funcs(bool bHang);

		template <typename TElem, typename TFVGeom>
		void register_func();
		/// \}
};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__VM_DISC_H__
