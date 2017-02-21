/*
 * synapse_container.h
 *
 *  Created on: Apr 17, 2016
 *      Author: lreinhardt, mbreit
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
