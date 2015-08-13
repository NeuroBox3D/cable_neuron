/*
 * VM_Disc.cpp
 *
 *  Created on: 26.11.2014
 *      Author: Pgottmann, mbreit
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

#include "../../diam_attachment_handler.h"	// attachment handling for diameter attachment
#include "channel_interface.h"

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
	#include "../../../synapse_handler/synapse_handler.h"
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
	#include "../../../synapse_distributor/synapse_distributor.h"
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

		// indices for unknowns
	    static const size_t _v_;
	    static const size_t _k_;
	    static const size_t _na_;
	    static const size_t _ca_;

	    static const size_t m_numb_ion_funcs;

	    // constants
	    const number R;	///< universal gas constant
	    const number F; ///< Faraday constant
	
	public:
		///	constructor
		VMDisc(const char* subsets, const number init_time = 0.0);

		///	destructor
		virtual ~VMDisc() {};

		// /////////////////////
		// setting parameters //
		// /////////////////////
		/// set constant diameter in units of m
		void set_diameter(const number d);

		/// set specific resistance in units of mV/(C/ms)*m = 1e-6 Ohm*m
		void set_spec_res(number val);

		/// set specific capacity in units of C/mV/m^2 = 1e-3 F/m^2
		void set_spec_cap(number val);

		/// set outer ion concentrations (in units of mol/m^3 = mM)
		///	\{
		void set_k_out(number value);
		void set_na_out(number value);
		void set_ca_out(number value);
		/// \}

		/// set diffusion coefficients in units of m^2/ms
		void set_diff_coeffs(const std::vector<number>& diff_coeffs);

		/// set Nernst potentials for ion species in units of mV
		///	\{
		void set_ek(number value);
		void set_ena(number value);
		void set_eca(number value);
		/// \}

		/// set temperature in units of K
		void set_temperature(number kelvin);

		/// set temperature in units of degrees C
		void set_temperature_celsius(number cels);

		// /////////////////////
		// getting parameters //
		// /////////////////////
		/// get constant diameter in units of m
		number diameter();

		/// get specific resistance in units of mV/(C/ms)*m = 1e-6 Ohm*m
		number spec_res();

		/// get specific capacity in units of C/mV/m^2 = 1e-3 F/m^2
		number spec_cap();

		/// get outer ion concentrations (in units of mol/m^3 = mM)
		///	\{
		number k_out();
		number na_out();
		number ca_out();
		/// \}

		/// get diffusion coefficients in units of m^2/ms
		const std::vector<number>& diff_coeffs();

		/// get Nernst potentials for ion species in units of mV
		///	\{
		number ek();
		number ena();
		number eca();
		/// \}

		/// get temperature in units of K
		number temperature();

		/// get temperature in units of degrees C
		number temperature_celsius();

		// ////////////////////////////
		// setters for functionality //
		// ////////////////////////////
		/// set influx position accuracy
		void set_influx_ac(number influx_ac);

		/// set influx params (flux value, coordinates, beginning, duration)
		void set_influx(number Flux, number x, number y, number z, number beg, number dur);

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		void set_synapse_handler(SmartPtr<synapse_handler::NETISynapseHandler<TDomain> > sh);
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
		void set_synapse_distributor(SmartPtr<SynapseDistributor> sd);
#endif

		/// adding a channel
		void add_channel(SmartPtr<IChannel<TDomain> > Channel);

		// ////////////////////////////////
		// getters for functional values //
		// ////////////////////////////////
		/// TODO: Writting function with position dependence
		/// functions to get different ion fluxes
		///	\{
		/*number flux_k();
		number flux_na();
		number flux_ca();
		number flux_v();*/
		/// \}

		/// get current time
		number time();

		/// get membrane potential at a vertex
		number get_vm(Vertex* vrt) const;

		/// write all gating values for a position to file
		void write_gatings_for_position(number x, number y, number z, std::string pfad);

		/// sets Vars for writing output
		void set_output(bool output, number gating_x, number gating_y, number gating_z, std::string gating_pfad);

		/// estimate time step size for next step
		template <typename TVector>
		number estimate_cfl_cond(ConstSmartPtr<TVector> u);

	public:
		// ///////////////////////////
		// inherited from IElemDisc //
		// ///////////////////////////
		/// @copydoc IElemDisc::approximation_space_changed()
		virtual void approximation_space_changed();

		///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		/// prepare the time step
		void prep_timestep(number time, VectorProxyBase* up);

		///	prepares the loop over all elements
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	prepares the element for assembling
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID id, const MathVector<dim> vCornerCoords[]);

		///	assembles stiffness part of local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles mass part of local defect
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		/// assembles jacobian of stiffness part
		template<typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		/// assembles jacobian of mass part
		template<typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	finishes the loop over all elements
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	protected:
		///	register utils
		///	\{
		void register_all_funcs(bool bHang);

		template <typename TElem, typename TFVGeom>
		void register_func();
		/// \}

	protected:
		ANumber m_aDiameter;										///< dendritic diameter attachment
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaDiameter;		///< dendritic diameter attachment accessor
		DiamAttachmentHandler m_dah;								///< handler for multigrid usage of diameter attachment
		number m_constDiam;											///< constant diameter (if set)
		bool m_bConstDiamSet;										///< whether const diameter is set

	private:
		number m_spec_res;					///< specific resistance in units of mV/(C/ms)*m = 1e-6 Ohm*m
		number m_spec_cap;					///< specific capacity in units of C/mV/m^2 = 1e-3 F/m^2

		number m_k_out;     				///< outer [K] in units of  mol/m^3 = mM
		number m_na_out;     				///< outer [Na] in units of  mol/m^3 = mM
		number m_ca_out;     				///< outer [Ca] in units of  mol/m^3 = mM

		std::vector<number> m_diff;			///< vector for diffusion constants of K, Na and Ca in units of m^2/ms

		number m_ek;						///< reversal potential K in units of mV
		number m_ena;						///< reversal potential Na in units of mV
		number m_eca;						///< reversal potential Ca in units of mV

		number m_temperature;				///< temperature in units of K


		number m_influx_ac;

		bool m_output;
		number m_gating_x, m_gating_y, m_gating_z;
		std::string m_gating_pfad;


	protected:
		std::vector<number> m_flux_value, m_beg_flux, m_dur_flux;		///< values describing influxes
		std::vector<MathVector<dim> > m_coords;							///< vector for influx coordinates x, y, z

#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
		SmartPtr<synapse_handler::NETISynapseHandler<TDomain> > m_spSH;	///< synapse handler
#endif

#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
		SmartPtr<SynapseDistributor> m_spSD;							///< synapse distributor
#endif

		std::vector<SmartPtr<TIChannel> > m_channel;					///< list of channels

		number m_v, m_na, m_k, m_ca;

		const number m_init_time;						///< time of initialization
		number m_time;									///< current time (the old solution is valid for)

		SmartPtr<CPUAlgebra::vector_type> m_spUOld;		///< old solution vector


		bool m_bNonRegularGrid;				///< current regular grid flag
		bool m_bLocked;						///< flag indicating whether approximation space has been set

		std::vector<SmartPtr<TIChannel> > m_channelsOnCurrSubset;
		std::vector<std::vector<size_t> > m_vvCurrChWFctInd;
		std::vector<number> m_currVrtValues[domain_traits<dim>::MaxNumVerticesOfElem];
};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__VM_DISC_H__
