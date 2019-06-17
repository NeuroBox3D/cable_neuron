/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Lukas Reinhardt
 * Creation date: 2016-03-23
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_HANDLER_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_HANDLER_H

// includes
#include <map>
#include <vector>

// ug
#include "../cable_disc/cable_equation.h"
#include "lib_grid/lib_grid.h"
#include "lib_grid/global_attachments.h"
#include "synapse_attachment_handler.h"
#include "lib_grid/algorithms/serialization.h" // GeomObjDataSerializer
#include "synapses/base_synapse.h"
#include "synapses/pre_synapse.h"
#include "synapses/post_synapse.h"

namespace ug {
namespace cable_neuron {

template <typename TDomain>
class CableEquation;
template <typename TDomain>
class ICableMembraneTransport;

namespace synapse_handler {

/**
 * Provides iterator-like access to synapses of type TSyn.
 * Use SynapseHandler::begin<TSyn> and SynapseHandler::end<TSyn>
 * to get iterators
 */
template <typename TSyn>
class SynapseIter
{
	std::vector<IBaseSynapse*>::iterator m_it;
public:
    SynapseIter(std::vector<IBaseSynapse*>::iterator it) : m_it(it) {}
    SynapseIter(const SynapseIter& sit) : m_it(sit.m_it) {}

	TSyn* operator*() {
		return dynamic_cast<TSyn*>(*m_it);
	}

	SynapseIter<TSyn>& operator++() {
		++m_it;
		return *this;
	}

	bool operator!=(const SynapseIter& rhs) {
		return m_it != rhs.m_it;
	}

	/// wrappers for lua binding
	void next() {operator++();}
	bool inequal(ConstSmartPtr<SynapseIter> rhs) {return operator!=(*rhs);}
};

/**
 * SynapseHandler is the central management class for everything related to Synapses.
 * Provides iterators for each registered synapse type.
 */
template <typename TDomain>
class SynapseHandler
{
private:
	typedef Attachment<std::vector<IBaseSynapse*> > AVSynapse;
	Grid::VertexAttachmentAccessor<APosition> m_aaPosition;


	AVSynapse m_aSyn;
	Grid::EdgeAttachmentAccessor<AVSynapse> m_aaSyn;

	bool m_bInited;
	SmartPtr<MultiGrid> m_spGrid;
	SmartPtr<CableEquation<TDomain> > m_spCEDisc;
	SynapseAttachmentHandler m_ssah;
	SmartPtr<ApproximationSpace<TDomain> > m_spApprox;

	std::vector<IBaseSynapse*> m_vAllSynapses;

	std::map<synapse_id, IPreSynapse*> m_mPreSynapses; // " ... "
	std::map<synapse_id, IPostSynapse*> m_mPostSynapses; // " ... "

	std::map<synapse_id, IPreSynapse*> m_mActivePreSynapses; // internal memory of currently active pre-synapses
	std::map<synapse_id, IPostSynapse*> m_mActivePostSynapses; // internal memory of currently active post-synapses
	std::map<synapse_id, Edge*> m_mPreSynapseIdToEdge; // maps presynapses to their edge
	std::map<synapse_id, Edge*> m_mPostSynapseIdToEdge; // maps presynapses to their edge

	SmartPtr<GeomObjDataSerializer<Edge> > m_spSAS;
	MessageHub::SPCallbackId m_spGridDistributionCallbackID;

	/// alpha synapse timing params
    /// @{
    number m_primary_alpha_onset_mean;
    number m_primary_alpha_onset_dev;
    number m_primary_alpha_tau_mean;
    number m_primary_alpha_tau_dev;
    number m_primary_alpha_peak_cond_mean;
    number m_primary_alpha_peak_cond_dev;
    bool m_primary_alpha_constSeed;
    /// @}

    /// bi-exp primary synapse timing params
    /// @{
    number m_prim_biexp_onset_mean;
    number m_prim_biexp_onset_dev;
    number m_prim_biexp_tau1_mean;
    number m_prim_biexp_tau1_dev;
    number m_prim_biexp_tau2_mean;
    number m_prim_biexp_tau2_dev;
    number m_prim_biexp_peak_cond_mean;
    number m_prim_biexp_peak_cond_dev;
    bool m_prim_biexp_constSeed;
    /// @}

    ///// FIXME: This should be refactored into something accordingly I suggest:
    /**
     * struct ballTimings {
     *	 template <typename TSyn>
     *	 void handle_ball_timing(TSyn begin, TSyn end) {
     *	 	/// handle timing for synapse depending on synapse type
     *	 	/// e.g. iterate only over stimulation balls which correspond to synapse type
     *	 }
     *	 private:
     *	 typedef std::vector<std::pair<std::vector<number>, std::vector<number> > >  TIMINGS;
     *	 std::pair<SYNAPSE_TYPE, std::vector<TIMINGS> > m_pTimings;
     *	 void add_ball_timing(const TIMINGS& timing, const SYNAPSE_TYPE& type) {
     *	 	 m_pTimings[type].push_back(timing);
     *	 }
     *
     *  Then: synapse_handler.cpp can call handle_ball_timings<SynapseType>(it_begin, it_end);
     *  and all responsibilities for the ball timings are hidden properly...
     *
     */
    /// alpha synapse ball timings
    std::vector<std::pair<std::vector<number>, std::vector<number> > > m_vTimingAlphaBalls;

    /// biexp synapse ball timings
    std::vector<std::pair<std::vector<number>, std::vector<number> > > m_vTimingBiExpBalls;

private:
	// do not use copy ctor
	SynapseHandler(const SynapseHandler& sh);

	/**
	 * Fills the member m_vAllSynapses with every synapse currently on the grid.
	 * Called by grid_first_available.
	 */
	void collect_synapses_from_grid();

	/// used for sorting
	struct __comp{
		bool operator() (IBaseSynapse* a, IBaseSynapse* b) {return a->type() < b->type();}
	};


public:
	/// constructor
	SynapseHandler();

	/// destructor
	virtual ~SynapseHandler() {}


    /// set the CableEquation object for which this synapse handler handles the synapses
	void set_ce_object(SmartPtr<CableEquation<TDomain> > disc) {
		UG_COND_THROW(m_bInited, "The CableEquation object associated to this synapse handler "
					  "must not be changed\nafter addition of the original CableEquation object "
					  "to the domain discretization.");
		m_spCEDisc = disc;
	}


    /**
     * @brief sets alpha post synapses activity pattern for primary synapses (randomly)
     *
     * The settings are only written to the base level. A call to propagate_synapses_to_levels()
     * is required afterwards to ensure that all grid levels are up to date.
     *
     * @param onset             average start time of activity
     * @param tau               average duration of activity
     * @param onset_dev         deviation of start time
     * @param tau_dev           deviation of duration
     * @param peak_cond         maximal conductivity
     * @param constSeed         if true: take 0 as seed; if false: take time-dependent seed (default: true)
     */
    void set_activation_timing_alpha
    (
        number onset,
        number tau,
        number peak_cond,
        number onset_dev,
        number tau_dev,
        number peak_cond_dev,
        bool constSeed
    );
    void set_activation_timing_alpha(number onset, number tau, number onset_dev, number tau_dev, number peak_cond, bool constSeed)
        {set_activation_timing_alpha(onset, tau, peak_cond, onset_dev, tau_dev, 0.0, constSeed);}
    void set_activation_timing_alpha(number onset, number tau, number onset_dev, number tau_dev, number peak_cond)
        {set_activation_timing_alpha(onset, tau, peak_cond, onset_dev, tau_dev, 0.0, true);}
    void set_activation_timing_alpha(number onset, number tau, number onset_dev, number tau_dev)
        {set_activation_timing_alpha(onset, tau, 6e-4, onset_dev, tau_dev, 0.0, true);}

    /**
     * @brief sets exp2 post synapses activity pattern for primary synapses (randomly)
     *
     * The settings are only written to the base level. A call to propagate_synapses_to_levels()
     * is required afterwards to ensure that all grid levels are up to date.
     *
     * @param onset_mean        average onset of activity
     * @param tau1_mean         average tau1 time constant
     * @param tau2_mean         average tau2 time constant
     * @param onset_dev         deviation of onset
     * @param tau1_dev          deviation of tau1
     * @param tau2_dev          deviation of tau2
     * @param peak_cond_mean    average of maximal conductivity
     * @param peak_cond_dev     deviation of maximal conductivity
     * @param constSeed         if true: take 0 as seed; if false: take time-dependent seed (default: true)
     */
    void set_activation_timing_biexp
    (
        number onset_mean,
        number tau1_mean,
        number tau2_mean,
        number peak_cond_mean,
        number onset_dev,
        number tau1_dev,
        number tau2_dev,
        number peak_cond_dev,
        bool constSeed
    );

    void set_activation_timing_biexp(number onset_mean, number tau1_mean, number tau2_mean,
        number onset_dev, number tau1_dev, number tau2_dev, number peak_cond, bool constSeed)
    {
        set_activation_timing_biexp(onset_mean, tau1_mean, tau2_mean, peak_cond,
                                    onset_dev, tau1_dev, tau2_dev, 0.0, constSeed);
    }
    void set_activation_timing_biexp(number onset_mean, number tau1_mean, number tau2_mean,
        number onset_dev, number tau1_dev, number tau2_dev, number peak_cond)
    {
        set_activation_timing_biexp(onset_mean, tau1_mean, tau2_mean, peak_cond,
                                    onset_dev, tau1_dev, tau2_dev, 0.0, true);
    }
    void set_activation_timing_biexp(number onset_mean, number tau1_mean, number tau2_mean,
        number onset_dev, number tau1_dev, number tau2_dev)
    {
        set_activation_timing_biexp(onset_mean, tau1_mean, tau2_mean, 6e-4,
                                    onset_dev, tau1_dev, tau2_dev, 0.0, true);
    }


    /**
     * @brief Adds an activation timing for a ball region.
     * Specify alpha synapses timings by providing a vector with start time,
     * duration, start time deviation, duration deviation & peak conductance
     *
     * @param alpha_timings
     * @param ball
     */
    void add_activation_timing_alpha_ball
    (
        const std::vector<number>& alpha_timings,
        const std::vector<number>& ball
    );


   /**
     * @brief Adds an activation timing for a ball region.
     * Specify exp2 synapses timings by providing a vector with onset time,
     * tau1 mean, tau2 mean, peak conductance, onset deviation tau1 deviation
     * tau2 deviation and deviation of peak conductance
     *
     * @param biexp_timings
     * @param ball
     */
    void add_activation_timing_exp2_ball
    (
        const std::vector<number>& biexp_timings,
        const std::vector<number>& ball
    );



    /// all synapses located on a specified edge
    const std::vector<IBaseSynapse*>& get_synapses_on_edge(Edge* e) const {return m_aaSyn[e];}


    /// Dummy interface method for compatibility reasons, wraps current_on_edge.
    bool synapse_on_edge(const Edge* edge, size_t scv, number time, number& current) const;


    /// total synaptic current for sub-control volume scv on edge e at time t
    number current_on_edge(const Edge* e, size_t scv, number t) const;

    //
    number current(synapse_id) const;

#if 0
    /**
     * Returns vector of currents, and a vector of their corresponding synapse id's
     * at time t on neuron nid.
     */
    void get_currents(number t, int nid, std::vector<number>& vCurrOut, std::vector<synapse_id>& vSidOut);
#endif

    /**
     * @brief Functionality executed when the grid is first known to synapse handler.
     *
     * This method will be called by the assigned CableEquation object when its approximation space is valid,
     * i.e., exactly when the CableEquation is added to the domain discretization and its method
     * CableEquation<TDomain>::approximation_space_changed() is called.
     *
     * It performs several grid-dependent init operations,
     * notably, it gathers every synapse in m_vAllSynapses for quick access/iterators.
     */
    void grid_first_available();


    /**
     * @brief update synaptic information
     * This method is to be called (by the associated vmDisc object) whenever the Vm values
     * or the approximation space have changed.
     * It will then make sure that every proc has all and current necessary information.
     */
    void update_presyn(number time);

    // unused
    //const std::vector<synapse_id> active_presynapses() const;

    /// getter for active synapse IDs and their currents
    void active_postsynapses_and_currents
	(
		std::vector<synapse_id>& vActSynOut,
		std::vector<number>& vSynCurrOut,
		const std::vector<uint>& vNid,
		const Grid::VertexAttachmentAccessor<Attachment<uint> >& aaNID,
		number time
	) const;


    Edge* postsyn_edge(synapse_id postSynID) const;
    Edge* presyn_edge(synapse_id preSynID) const;

    IPreSynapse* pre_synapse(synapse_id id) const;
    IPostSynapse* post_synapse(synapse_id id) const;

    /**
	 * Returns a begin iterator to the desired synapse type.
	 * Templates have to be specialized for use in LUA.
	 * TODO: Improve at least to log(N) complexity (m_vAllSynapses is type-ordered!)
	 *       Better yet, save offsets for types for constant access.
	 */
	template <typename TSyn>
	SynapseIter<TSyn> begin() {
		std::vector<IBaseSynapse*>::iterator it = m_vAllSynapses.begin();
		for(; it != m_vAllSynapses.end(); ++it ) {
			if(dynamic_cast<TSyn*>(*it) ) {//found the first TSyn*
				 break;
			}
		}
		return SynapseIter<TSyn>(it);
	}

	/// wrapper for lua binding
	template <typename TSyn>
	SmartPtr<SynapseIter<TSyn> > begin_wrapper() {
        return SmartPtr<SynapseIter<TSyn> >(new SynapseIter<TSyn>(begin<TSyn>()));
    }

	/**
	 * Returns an end iterator to the desired synapse type.
	 * Templates have to be specialized for use in LUA.
	 * TODO: Improve at least to log(N) complexity (m_vAllSynapses is type-ordered!)
	 *       Better yet, save offsets for types for constant access.
	 */
	template <typename TSyn>
	SynapseIter<TSyn> end() {
		std::vector<IBaseSynapse*>::iterator it = m_vAllSynapses.begin();
		for(; it != m_vAllSynapses.end(); ++it ) {
			if(dynamic_cast<TSyn*>(*it) ) {//found the first TSyn*
				//std::cout << *it << std::endl;
				break;
			}
		}

		while(it != m_vAllSynapses.end()) {
			if(!dynamic_cast<TSyn*>(*it) ) {//found last TSyn*
				break;
			}
			++it;
		}
		return SynapseIter<TSyn>(it);
	}

    /// wrapper for lua binding
    template <typename TSyn>
    SmartPtr<SynapseIter<TSyn> > end_wrapper() {
        return SmartPtr<SynapseIter<TSyn> >(new SynapseIter<TSyn>(end<TSyn>()));
    }

	/// Prints list with all synapses and their parameters/id's etc.
	void show_status(number t);


protected:
    /// performs the actual setting of activation timing when the grid is available
    void set_activation_timing_with_grid();

    /**
     * @brief explicitly deactivates all bi-exponential synapses
     *
     * This can be used to ascertain that no synapse is active without being activated,
     * which is necessary since SynapseInfo::m_onset value might be incorrectly set
     * by topology importer.
     * The deactivation is only performed on the base level. A call to propagate_synapses_to_levels()
     * is required afterwards to ensure that all grid levels are up to date.
     */
    void deactivate_all_postSyns();

    /// callback to register at message hub
    void grid_distribution_callback(const GridMessage_Distribution& gmd);
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_HANDLER_H
