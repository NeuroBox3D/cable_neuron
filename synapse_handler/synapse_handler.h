/*!
 * \file synapse_handler.h
 * \brief synapse handler interface and implementations
 *
 *	@authors stephanmg, mbreit, mstepnie
 *  \date last major update, 22/07/2015 (mbreit, mstepnie)
 */

// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SYNAPSE_HANDLER_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SYNAPSE_HANDLER_H__

// includes
#include <map>
#include <string>
#include <algorithm>

// boost includes for random numbers
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

// ug
#include "common/log.h"
#include "lib_grid/lib_grid.h"
#include "lib_grid/global_attachments.h"
#include "lib_grid/tools/copy_attachment_handler.h"

// synapses
#include "synapse.h"
//#include "synapse_factory.h"
//#include "selector/synapse_selector.h"
//#include "selector/synapse_mapping.h"
//#include "function/synapse_factory_producer.h"
#include "function/types.h"
#include "grid/synapse_info_attachment_handler.h"
#include "grid/synapse_info_io_traits.h"

// synapse distributor and cable equation
#include "../synapse_distributor/synapse_distributor.h"
#include "../cable_disc/cable_equation.h"


namespace ug {
namespace cable_neuron {


template <typename TDomain>
class CableEquation;
template <typename TDomain>
class ICableMembraneTransport;


namespace synapse_handler {

/*! \defgroup sh_plugin Synapse Handler plugin
 * \ingroup plugins_experimental
 * \{
 */

////////////////////////////////////////////////////////////////////////
/// SYNAPSE HANLDLER INTERFACE
////////////////////////////////////////////////////////////////////////
/*!
 * \brief synapse handler interface
 */
template <typename TDomain>
class ISynapseHandler {
	public:
		/*!
		 * \brief identify the synapse handler
		 * \remark this should be a unique identifier
		 *
		 * \return \c std::string
		 */
		virtual std::string name() const = 0;

		/*!
		 * \brief indicate if synapse is present at the given coordinates and time.
		 * \remark In addition the outward current is returned.
		 *
		 * \param[in] vec 		the position coordinates
		 * \param[in] time 		the current time
		 * \param[out] current 	outward current
		 *
		 * \return \c bool
		 */
		virtual bool synapse_at_location(const MathVector<TDomain::dim>& vec, number time, number& current) const = 0;

		/*!
		 * \brief indicate if synapse is present on the given sub-control
		 * volume of the given edge at the given time.
		 * \remark in addition the outward current is to be returned.
		 *
		 * In addition it should be checked on which site the
		 * synapse locates, since in FV we assemble the SCVs,
		 * where these edges are only "half-edges" then.
		 *
		 * \param[in] edge		the edge to test
		 * \param[in] scv 		the
		 * \param[in] time 		the current time
		 * \param[out] current	outward current
		 *
		 * \return \c bool
		 */
		virtual bool synapse_on_edge(const Edge* const edge, size_t scv, number time, number& current) = 0;

		/*!
		 * \brief mandatory vdtor
		 */
		virtual ~ISynapseHandler() {}
};



/*!
 * \brief synapse handler implementation for the synapse distributor plugin
 * \author stephanmg
 */
template <typename TDomain>
class SynapseDistributorSynapseHandler : public ISynapseHandler<TDomain> {
private:
	ConstSmartPtr<SynapseDistributor> m_spSD;

public:
	SynapseDistributorSynapseHandler(ConstSmartPtr<SynapseDistributor> spSD) : m_spSD(spSD) { }
	SynapseDistributorSynapseHandler() { }

	~SynapseDistributorSynapseHandler() {
	}

	///@copydoc ISynapseHandler<TDomain>::name()
	std::string name() const;

	///@copydoc ISynapseHandler<TDomain>::synapse_on_edge()
	bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current);

	///@copydoc ISynapseHandler<TDomain>::synapse_at_location()
	bool synapse_at_location(const ug::MathVector<TDomain::dim>& vec, number time, number& current) const;

	/// set underlying synapse distributor
	void set_sd(ConstSmartPtr<SynapseDistributor> sd);
};



/*!
 * \brief synapse handler for geometries created using the neuronal topology importer
 * \authors stephanmg, mbreit
 */

template <typename TDomain>
class NETISynapseHandler : public ISynapseHandler<TDomain>
{
	public:

		typedef Attachment<std::vector<SynapseInfo> > AVSynapse;

		/// constructor
		NETISynapseHandler();

		/// destructor
		virtual ~NETISynapseHandler() {};

		/// set presnyaptic subset name
		void set_presyn_subset(const char* presynSubset);

		/// set the CableEquation object for which this synapse handler handles the synapses
		void set_ce_object(SmartPtr<CableEquation<TDomain> > disc);

		/**
		 * @brief sets alpha synapses activity pattern (randomly)
		 *
		 * The settings are only written to the base level. A call to propagate_synapses_to_levels()
		 * is required afterwards to ensure that all grid levels are up to date.
		 *
		 * @param start_time		average start time of activity
		 * @param duration			average duration of activity
		 * @param start_time_dev	deviation of start time
		 * @param duration_dev		deviation of duration
		 * @param peak_cond			maximal conductivity
		 * @param constSeed			if true: take 0 as seed; if false: take time-dependent seed (default: true)
		 */
		void set_activation_timing
		(
			number start_time,
			number duration,
			number start_time_dev,
			number duration_dev,
			number peak_cond,
			bool constSeed
		);
		void set_activation_timing(number start_time, number duration, number start_time_dev, number duration_dev, number peak_cond)
			{set_activation_timing(start_time, duration, start_time_dev, duration_dev, peak_cond, true);}
		void set_activation_timing(number start_time, number duration, number start_time_dev, number duration_dev)
			{set_activation_timing(start_time, duration, start_time_dev, duration_dev, 6e-4, true);}

		/**
		 * @brief sets alpha synapses activity pattern (randomly)
		 *
		 * The settings are only written to the base level. A call to propagate_synapses_to_levels()
		 * is required afterwards to ensure that all grid levels are up to date.
		 *
		 * @param onset_mean		average onset of activity
		 * @param tau1_mean			average tau1 time constant
		 * @param tau2_mean			average tau2 time constant
		 * @param onset_dev			deviation of onset
		 * @param onset_tau1		deviation of tau1
		 * @param onset_tau2		deviation of tau2
		 * @param peak_cond			maximal conductivity
		 * @param constSeed			if true: take 0 as seed; if false: take time-dependent seed (default: true)
		 */
		void set_activation_timing_biexp
		(
			number onset_mean,
			number tau1_mean,
			number tau2_mean,
			number onset_dev,
			number tau1_dev,
			number tau2_dev,
			number peak_cond,
			bool constSeed
		);
		void set_activation_timing_biexp(number onset_mean, number tau1_mean, number tau2_mean,
			number onset_dev, number tau1_dev, number tau2_dev, number peak_cond)
		{
			set_activation_timing_biexp(onset_mean, tau1_mean, tau2_mean,
										onset_dev, tau1_dev, tau2_dev, peak_cond, true);
		}
		void set_activation_timing_biexp(number onset_mean, number tau1_mean, number tau2_mean,
			number onset_dev, number tau1_dev, number tau2_dev)
		{
			set_activation_timing_biexp(onset_mean, tau1_mean, tau2_mean,
										onset_dev, tau1_dev, tau2_dev, 6e-4, true);
		}


		/**
		 * @brief functionality executed when the grid is first known to synapse handler
		 *
		 * This method will be called be the assigned CableEquation object when its approximation space is valid,
		 * i.e. exactly when the CableEquation is added to the domain discretization and its method
		 * CableEquation<TDomain>::approximation_space_changed() is called.
		 *
		 */
		void grid_first_available();


		/**
		 *	@brief update synaptic information
		 *	This method is to be called (by the associated vmDisc object) whenever the Vm values
		 *	or the approximation space have changed.
		 *	It will then make sure that every proc has all and current necessary information.
		 */
		virtual void update_presyn();


		/// @copydoc ISynapseHandler::synapse_on_edge()
		bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current);

		/// @copydoc ISynapseHandler::synapse_at_location()
		bool synapse_at_location(const MathVector<TDomain::dim>& vec, number time, number& current) const;

		/// @copydoc ISynapseHandler::name()
		std::string name() const;

		/**
		 * @brief prints out synapse statistics
		 * @param soma_si subset index for soma subset
		 *
		 * @note this only works correctly if no neuron is cut by parallel distribution
		 * @todo dirty implementation -- should be re-implemented properly somewhere else
		 */
		void print_synapse_statistics(size_t soma_si);

		/**
		 * @brief writes synapse and somatic activity to file
		 *
		 * @note this only works correctly if no neuron is cut by parallel distribution
		 * @todo dirty implementation -- should be re-implemented properly somewhere else
		 */
		void write_activity_to_file(const std::string& fileName, number time);

	protected:
		struct MyIndexSort
		{
			MyIndexSort(const std::vector<size_t>& ind) : m_ind(ind) {}
			bool operator() (size_t a, size_t b) {
				UG_COND_THROW(a >= m_ind.size(), "Requested index "<< a
						<< " from vector containing only "<< m_ind.size()<<"."<<std::endl);
				UG_COND_THROW(b >= m_ind.size(), "Requested index "<< b
						<< " from vector containing only "<< m_ind.size()<<"."<<std::endl);
				return m_ind.at(a) < m_ind.at(b);}

			private:
				const std::vector<size_t>& m_ind;
		};

		struct LocalESynapseInfo
		{
			LocalESynapseInfo(Edge* e, size_t i, size_t f, size_t t)
				: edge(e), ind(i), from(f), to(t) {}

			Edge* edge;
			size_t ind;
			size_t from;
			size_t to;
		};
		struct LocalASynapseInfo
		{
			LocalASynapseInfo(Edge* e, size_t i, size_t t)
				: edge(e), ind(i), to(t) {}

			Edge* edge;
			size_t ind;
			size_t to;
		};

		bool is_active(LocalASynapseInfo locInfo, number t)
		{
			const std::vector<SynapseInfo>& vSI = m_aaSynapseInfo[locInfo.edge];
			const SynapseInfo& info = vSI[locInfo.ind];

			typedef synapse_traits<AlphaSynapse> STA;
			typedef synapse_traits<Exp2Syn> STB;
			typedef synapse_traits<void> STV;

			UG_COND_THROW(STV::type(info) != ALPHA_SYNAPSE, "Wrong synapse type! Expected ALPHA_SYNAPSE.")

			return t >= STA::onset(info) && t <= STA::onset(info) + 5*STA::tau(info);
		}
		bool is_active(LocalESynapseInfo locInfo, number t)
		{
			std::vector<SynapseInfo>& vSI = m_aaSynapseInfo[locInfo.edge];
			const SynapseInfo& info = vSI[locInfo.ind];

			typedef synapse_traits<AlphaSynapse> STA;
			typedef synapse_traits<Exp2Syn> STB;
			typedef synapse_traits<void> STV;

			UG_COND_THROW(STV::type(info) != EXP2_SYNAPSE, "Wrong synapse type! Expected EXP2_SYNAPSE.")

			return STB::activated(info);
		}

		/**
		 * @brief explicitly deactivates all bi-exponential synapses
		 *
		 * This can be used to ascertain that no synapse is active without being activated,
		 * which is necessary since SynapseInfo::m_onset value might be incorrectly set
		 * by topology importer.
		 * The deactivation is only performed on the base level. A call to propagate_synapses_to_levels()
		 * is required afterwards to ensure that all grid levels are up to date.
		 */
		void deactivate_all_biexp();

		/// performs the actual setting of activation timing when the grid is available
		void set_activation_timing_with_grid();

		void print_biexp_info(size_t syn_index);

		/**
		 * @brief finds out how many presynaptic indices are present in the geometry
		 *
		 * The local vector m_vPresynVmValues for storage of presynaptic potential
		 * values is resized accordingly.
		 *
		 * propagate_synapses_to_levels() has to be called beforehand!
		 */
		void resize_presyn_vector();

		/**
		 * @brief checks, if bi-exp. synapses are present in the geometry
		 */
		bool has_EXP2_SYNAPSE();

	private:
		/// CableEquation
		SmartPtr<CableEquation<TDomain> > m_spCEDisc;

		/// approx space
		SmartPtr<ApproximationSpace<TDomain> > m_spApprox;

		/// grid
		SmartPtr<MultiGrid> m_spGrid;

		AVSynapse m_aSynInfo;			///< synapse connection info attachment
		AUInt m_aPSI; 			///< presynaptic index attachment

		Grid::EdgeAttachmentAccessor<AVSynapse> m_aaSynapseInfo;
		Grid::VertexAttachmentAccessor<AUInt> m_aaPSI;

		/// presynaptic subset name
		std::string m_presynSubset;

		/// presynaptic subset index
		int m_presynSI;

		/// handling of presynaptic vm values (in units of V)
		std::vector<number> m_vPresynVmValues;

		/// alpha synapse timing params
		///	@{
		number m_start_time;
		number m_start_time_dev;
		number m_duration;
		number m_duration_dev;
		number m_peak_cond;
		bool m_constSeed;
		/// @}

		/// bi-exp primary synapse timing params
		///	@{
		number m_prim_biexp_onset_mean;
		number m_prim_biexp_tau1_mean;
		number m_prim_biexp_tau2_mean;
		number m_prim_biexp_onset_dev;
		number m_prim_biexp_tau1_dev;
		number m_prim_biexp_tau2_dev;
		number m_prim_biexp_peak_cond;
		bool m_prim_biexp_constSeed;
		/// @}

		/// attachment handler for presynapse index attachment
		CopyAttachmentHandler<Vertex, AUInt> m_cah;

		/// attachment handler for vertex SynapseInfo attachment
		SynapseInfoAttachmentHandler m_siah;

		bool m_bInited;
		bool m_bEXP2_SYNAPSE;

		std::vector<LocalASynapseInfo> m_locASynInfo;
		std::vector<LocalESynapseInfo> m_locESynInfo;
		std::vector<Vertex*> m_vSomaVertices;
		std::vector<size_t> m_vNeuronType;
};


///<! \}

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SYNAPSE_HANDLER_H__
