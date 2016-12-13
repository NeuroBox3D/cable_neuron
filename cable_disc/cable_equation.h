/*
 * cable_equation.h
 *
 *  Created on: 26.11.2014
 *  Author: Pgottmann, mbreit
 */


#ifndef __UG__PLUGINS__CABLE_NEURON__CABLE_DISC__CABLE_EQUATION_H__
#define __UG__PLUGINS__CABLE_NEURON__CABLE_DISC__CABLE_EQUATION_H__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

// library intern headers
#include "../membrane_transport/cable_membrane_transport_interface.h"
#include "../util/diam_attachment_handler.h"	// attachment handling for diameter attachment
#include "../util/cable_ass_tuner.h"
#include "../split_synapse_handler/split_synapse_handler.h"
//#include "../synapse_handler/synapse_handler.h"


namespace ug {
namespace cable_neuron {

namespace synapse_handler
{
	template <typename TDomain>
	class SplitSynapseHandler;

	// forward declaration
//	template <typename TDomain>
//	class NETISynapseHandler;
}


// forward declaration
template <typename TDomain>
class ICableMembraneTransport;

template <typename TDomain>
class CableEquation
	: public IElemDisc<TDomain>
{
	public:
		/// type for channels
		typedef ICableMembraneTransport<TDomain> TIChannel;

		///	world dimension
		static const int dim = IElemDisc<TDomain>::dim;

		// indices for unknowns
	    static const size_t _v_;
	    static const size_t _k_;
	    static const size_t _na_;
	    static const size_t _ca_;

	    const size_t m_numb_ion_funcs;

	    // constants
	    const number R;	///< universal gas constant
	    const number F; ///< Faraday constant
	
	public:
		///	constructor
		CableEquation(const char* subsets, bool withConcs = true, number init_time = 0.0);

		///	destructor
		virtual ~CableEquation() {};

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
		// TODO: make that set_rev_pot(string ion, number value)
		void set_rev_pot_k(number value);
		void set_rev_pot_na(number value);
		void set_rev_pot_ca(number value);
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
		number conc_out(size_t ion_spec);
		/// \}

		/// get all flux-values
		number flux_ca();
		number flux_na();
		number flux_k();

		/// get diffusion coefficients in units of m^2/ms
		const std::vector<number>& diff_coeffs();

		/// get Nernst potentials for ion species in units of mV
		///	\{
		number rev_pot_k();
		number rev_pot_na();
		number rev_pot_ca();
		/// \}

		/// get temperature in units of K
		number temperature();

		/// get temperature in units of degrees C
		number temperature_celsius();

		// ////////////////////////////
		// setters for functionality //
		// ////////////////////////////

		/// setter for influx via subset
		void set_influx_subset(int influx_subset, number input, number dur, number start);

		/// set influx position accuracy
		void set_influx_ac(number influx_ac);

		/// set influx params (flux value, coordinates, beginning, duration)
		void set_influx(number Flux, number x, number y, number z, number beg, number dur);

		/// set synapse handler
		void set_synapse_handler(SmartPtr<synapse_handler::SplitSynapseHandler<TDomain> > sh);

		/// adding a channel
		void add(SmartPtr<ICableMembraneTransport<TDomain> > transportMechanism);

		// ////////////////////////////////
		// getters for functional values //
		// ////////////////////////////////
		/// TODO: Writing function with position dependence
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
		number vm(Vertex* vrt) const;

		/// get current assembling subset
		int current_subset_index() const;

		/// write all internal state values for a position to file
		void write_states_for_position(number x, number y, number z, std::string pfad);

		/// sets Vars for writing output
		void set_output_point_and_path(bool output, number x, number y, number z, std::string outPath);

		/// estimate time step size for next step
		template <typename TVector>
		number estimate_cfl_cond(ConstSmartPtr<TVector> u);

		/// get a vector of all surface vertices
		const std::vector<Vertex*>& surface_vertices() const;

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
		MessageHub::SPCallbackId m_spGridDistributionCallbackID;
		void grid_distribution_callback(const GridMessage_Distribution& gmd);

		void update_surface_verts();

	protected:
		ANumber m_aDiameter;										///< dendritic diameter attachment
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaDiameter;		///< dendritic diameter attachment accessor
		DiamAttachmentHandler m_dah;								///< handler for multigrid usage of diameter attachment
		number m_constDiam;											///< constant diameter (if set)
		bool m_bConstDiamSet;										///< whether const diameter is set

	private:
		number m_spec_res;					///< specific resistance in units of Ohm*m
		number m_spec_cap;					///< specific capacitance in units of F/m^2

		number m_k_out;     				///< outer [K] in units of  mol/m^3 = mM
		number m_na_out;     				///< outer [Na] in units of  mol/m^3 = mM
		number m_ca_out;     				///< outer [Ca] in units of  mol/m^3 = mM

		std::vector<number> m_diff;			///< vector for diffusion constants of K, Na and Ca in units of m^2/s

		number m_ek;						///< reversal potential K in units of mV
		number m_ena;						///< reversal potential Na in units of mV
		number m_eca;						///< reversal potential Ca in units of mV

		number m_eqConc_ca;					///< equilibrium concentration Ca
		number m_reactionRate_ca;			///< reaction rate in Ca source/sink term (units of 1/s)

		number m_temperature;				///< temperature in units of K

		number m_influx_ac;

		bool m_bOutput;
		number m_output_x, m_output_y, m_output_z;
		std::string m_outputPath;

		// vars for influx via subset
		std::vector<int> m_influx_subset;
		std::vector<number> m_vSubsetInflux, m_vSubsetInfluxStart, m_vSubsetInfluxDur;

	protected:
		std::vector<number> m_vCurrent, m_vCurrentStart, m_vCurrentDur;		///< values describing influxes
		std::vector<MathVector<dim> > m_vCurrentCoords;					///< vector for influx coordinates x, y, z

		//synapse handler
		SmartPtr<synapse_handler::SplitSynapseHandler<TDomain> > m_spSH;

		std::vector<SmartPtr<TIChannel> > m_channel;					///< list of channels

		number m_v, m_na, m_k, m_ca;

		const number m_init_time;						///< time of initialization
		number m_time;									///< current time (the old solution is valid for)

		SmartPtr<CPUAlgebra::vector_type> m_spUOld;		///< old solution vector


		bool m_bNonRegularGrid;				///< current regular grid flag
		bool m_bLocked;						///< flag indicating whether approximation space has been set

		/// vector of subset indices this disc is declared on
		std::vector<int> m_vSI;

		int m_si;
		std::vector<SmartPtr<TIChannel> > m_channelsOnCurrSubset;
		std::vector<std::vector<size_t> > m_vvCurrChWFctInd;
		std::vector<number> m_currVrtValues[domain_traits<dim>::MaxNumVerticesOfElem];

		std::vector<Vertex*> m_vSurfVrt;
};


} // namespace cable_neuron
} // namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__CABLE_DISC__CABLE_EQUATION_H__
