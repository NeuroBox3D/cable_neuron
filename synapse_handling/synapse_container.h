/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Lukas Reinhardt
 * Creation date: 2016-04-17
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

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_CONTAINER_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_CONTAINER_H

#include "common/types.h"
#include "synapses/base_synapse.h"   // IBaseSynapse, synapse_id
#include "synapses/alpha_post_synapse.h"
#include "synapses/exp2_post_synapse.h"
#include "synapses/onset_pre_synapse.h"
#include "synapses/threshold_pre_synapse.h"


namespace ug {
namespace cable_neuron {
namespace synapse_handler {


class ISynapseContainer {
public:
	virtual ~ISynapseContainer() {}
	virtual std::vector<IBaseSynapse*> get_synapses() = 0;
	virtual size_t size() = 0;
};


/**
 * Class for parameterization of a set of AlphaSynapses, for usage in lua-scripts.
 * Creates <num_synapses> pairs of a post- and presynapse, in that order.
 * Therefore the first ID will be <start_id>
 * and the last id will be <start_id> + 2*<num_synapses> - 1.
 */
class AlphaSynapses : ISynapseContainer
{
private:
	size_t m_start_id;
	number m_mean_onset;
	number m_dev_onset;
	number m_mean_gMax;
	number m_dev_gMax;
	number m_mean_tau;
	number m_dev_tau;
	number m_mean_rev;
	number m_dev_rev;

	std::vector<IBaseSynapse*> m_vSynapses;

public:
	/// constructor
	AlphaSynapses(const size_t start_id, const size_t num_synapses);

    /// constructor
	AlphaSynapses
	(
		const size_t start_id,
		const size_t num_synapses,
		number mean_onset,
		number dev_onset,
        number mean_gMax,
        number dev_gMax,
        number mean_tau,
        number dev_tau,
        number mean_rev,
        number dev_rev
    );

    virtual ~AlphaSynapses();

	void set_mean_onset(const number val);
	void set_dev_onset(const number val);
	void set_mean_gMax(const number val);
	void set_dev_gMax(const number val);
	void set_mean_tau(const number val);
	void set_dev_tau(const number val);
	void set_mean_rev(const number val);
	void set_dev_rev(const number val);

	size_t size();
	size_t start_id();

	/// assemble container
	std::vector<IBaseSynapse*> get_synapses();
};



/**
 * Class for parameterization of a set of Exp2Synapses, for usage in lua-scripts.
 * Creates <num_synapses> pairs of a post- and presynapse, in that order.
 * Therefore the first ID will be <start_id>
 * and the last id will be <start_id> + 2*<num_synapses> - 1.
 */
class Exp2Synapses : ISynapseContainer
{
private:
	size_t m_start_id;
	number m_mean_onset;
	number m_dev_onset;
	number m_mean_gMax;
	number m_dev_gMax;
	number m_mean_tau1;
	number m_dev_tau1;
	number m_mean_tau2;
	number m_dev_tau2;
	number m_mean_rev;
	number m_dev_rev;
	number m_mean_threshold;
	number m_dev_threshold;

	std::vector<IBaseSynapse*> m_vSynapses;

public:
	/// constructor
	Exp2Synapses(const size_t start_id, const size_t num_synapses);

    /// constructor
	Exp2Synapses
	(
        const size_t start_id,
        const size_t num_synapses,
        number mean_onset,
        number dev_onset,
        number mean_gMax,
        number dev_gMax,
        number mean_tau1,
        number dev_tau1,
        number mean_tau2,
        number dev_tau2,
        number mean_rev,
        number dev_rev,
        number mean_threshold,
        number dev_threshold
    );

	/// destructor
    virtual ~Exp2Synapses();

	void set_mean_onset(const number val);
	void set_dev_onset(const number val);
	void set_mean_gMax(const number val);
	void set_dev_gMax(const number val);
	void set_mean_tau1(const number val);
	void set_dev_tau1(const number val);
	void set_mean_tau2(const number val);
	void set_dev_tau2(const number val);
	void set_mean_rev(const number val);
	void set_dev_rev(const number val);
	void set_mean_threshold(const number val);
	void set_dev_threshold(const number val);

	size_t size();
	size_t start_id();


	/// assembles container
	std::vector<IBaseSynapse*> get_synapses();
};



class AlphaSynapsePair
{
    public:
        /// constructor
        AlphaSynapsePair();

        /// destructor
        ~AlphaSynapsePair();

        /// set id
        void set_id(size_t id);

        /// set onset
        void set_onset(number onset);

        /// set tau (and duration, which is 6*tau)
        void set_tau(number tau);

        /// set maximal conductance
        void set_gMax(number gMax);

        /// set reversal potential
        void set_reversal_potential(number ve);

        /// get pre-synapse
        OnsetPreSynapse* pre_synapse();

        /// get post-synapse
        AlphaPostSynapse* post_synapse();

    private:
        OnsetPreSynapse* m_pre;
        AlphaPostSynapse* m_post;
};



class Exp2SynapsePair
{
    public:
        /// constructor
        Exp2SynapsePair();

        /// destructor
        ~Exp2SynapsePair();

        /// set id
        void set_id(size_t id);

        /// set id
        void set_threshold(number t);

        /// set maximal conductance
        void set_gMax(number gMax);

        /// set reversal potential
        void set_reversal_potential(number ve);

        /// set tau1 and tau2 (and duration, which is determined by the taus)
        void set_taus(number tau1, number tau2);

        /// get pre-synapse
        ThresholdPreSynapse* pre_synapse();

        /// get post-synapse
        Exp2PostSynapse* post_synapse();

    private:
        ThresholdPreSynapse* m_pre;
        Exp2PostSynapse* m_post;
};



} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_CONTAINER_H
