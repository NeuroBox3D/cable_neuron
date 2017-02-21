/*
 * synapse_container.cpp
 *
 *  Created on: Apr 17, 2016
 *      Author: lreinhardt, mbreit
 */

#include "synapse_container.h"

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <cmath>


namespace ug {
namespace cable_neuron {
namespace synapse_handler {


AlphaSynapses::AlphaSynapses(const size_t start_id, const size_t num_synapses)
: m_start_id(start_id),
  m_mean_onset(0), m_dev_onset(0),
  m_mean_gMax(0), m_dev_gMax(0),
  m_mean_tau(0), m_dev_tau(0),
  m_mean_rev(0), m_dev_rev(0),
  m_vSynapses(std::vector<IBaseSynapse*>(2 * num_synapses))
{}

AlphaSynapses::AlphaSynapses
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
)
: m_start_id(start_id),
  m_mean_onset(mean_onset), m_dev_onset(dev_onset),
  m_mean_gMax(mean_gMax), m_dev_gMax(dev_gMax),
  m_mean_tau(mean_tau), m_dev_tau(dev_tau),
  m_mean_rev(mean_rev), m_dev_rev(dev_rev),
  m_vSynapses(std::vector<IBaseSynapse*>(2 * num_synapses))
{}


AlphaSynapses::~AlphaSynapses()
{
    // delete synapses in m_vSynapses which were created here
    for (size_t i = 0; i < m_vSynapses.size(); ++i)
        delete m_vSynapses[i];
}


void AlphaSynapses::set_mean_onset(const number val) {m_mean_onset=val;}
void AlphaSynapses::set_dev_onset(const number val) {m_dev_onset=val;}
void AlphaSynapses::set_mean_gMax(const number val) {m_mean_gMax=val;}
void AlphaSynapses::set_dev_gMax(const number val) {m_dev_gMax=val;}
void AlphaSynapses::set_mean_tau(const number val) {m_mean_tau=val;}
void AlphaSynapses::set_dev_tau(const number val) {m_dev_tau=val;}
void AlphaSynapses::set_mean_rev(const number val) {m_mean_rev=val;}
void AlphaSynapses::set_dev_rev(const number val) {m_dev_rev=val;}

size_t AlphaSynapses::size() {return m_vSynapses.size();}
size_t AlphaSynapses::start_id() {return m_start_id;}


std::vector<IBaseSynapse*> AlphaSynapses::get_synapses()
{
    boost::mt19937 onset_rng;
    boost::normal_distribution<number> onset_dist(m_mean_onset, m_dev_onset);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > onset_var(onset_rng, onset_dist);

    boost::mt19937 gmax_rng;
    boost::normal_distribution<number> gmax_dist(m_mean_gMax, m_dev_gMax);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > gmax_var(gmax_rng, gmax_dist);

    boost::mt19937 tau_rng;
    boost::normal_distribution<number> tau_dist(m_mean_tau, m_dev_tau);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau_var(tau_rng, tau_dist);

    boost::mt19937 rev_rng;
    boost::normal_distribution<number> rev_dist(m_mean_rev, m_dev_rev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > rev_var(rev_rng, rev_dist);


    // create a couple of synapses, pre and post synapse together.
    for(size_t i = 0; i < m_vSynapses.size(); i=i+2) {
        number onset = onset_var();
        number tau = tau_var();
        size_t syn_id = m_start_id + i/2;
        IBaseSynapse *s1 = new AlphaPostSynapse(syn_id, 0.0, gmax_var(), tau, rev_var());
        IBaseSynapse *s2 = new OnsetPreSynapse(syn_id, 0.0, onset, 6 * tau);
        m_vSynapses[i] = s1;
        m_vSynapses[i+1] = s2;
    }

    return m_vSynapses;
}




Exp2Synapses::Exp2Synapses(const size_t start_id, const size_t num_synapses)
: m_start_id(start_id),
  m_mean_onset(0), m_dev_onset(0),
  m_mean_gMax(0), m_dev_gMax(0),
  m_mean_tau1(0), m_dev_tau1(0),
  m_mean_tau2(0), m_dev_tau2(0),
  m_mean_rev(0), m_dev_rev(0),
  m_mean_threshold(0), m_dev_threshold(0),
  m_vSynapses(std::vector<IBaseSynapse*>(2 * num_synapses))
{}

Exp2Synapses::Exp2Synapses
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
)
: m_start_id(start_id),
  m_mean_onset(mean_onset), m_dev_onset(dev_onset),
  m_mean_gMax(mean_gMax), m_dev_gMax(dev_gMax),
  m_mean_tau1(mean_tau1), m_dev_tau1(dev_tau1),
  m_mean_tau2(mean_tau2), m_dev_tau2(dev_tau2),
  m_mean_rev(mean_rev), m_dev_rev(dev_rev),
  m_mean_threshold(mean_threshold), m_dev_threshold(dev_threshold),
  m_vSynapses(std::vector<IBaseSynapse*>(2 * num_synapses))
{}


Exp2Synapses::~Exp2Synapses()
{
    // delete synapses in m_vSynapses which were created here
    for (size_t i = 0; i < m_vSynapses.size(); ++i)
        delete m_vSynapses[i];
}


void Exp2Synapses::set_mean_onset(const number val) {m_mean_onset=val;}
void Exp2Synapses::set_dev_onset(const number val) {m_dev_onset=val;}
void Exp2Synapses::set_mean_gMax(const number val) {m_mean_gMax=val;}
void Exp2Synapses::set_dev_gMax(const number val) {m_dev_gMax=val;}
void Exp2Synapses::set_mean_tau1(const number val) {m_mean_tau1=val;}
void Exp2Synapses::set_dev_tau1(const number val) {m_dev_tau1=val;}
void Exp2Synapses::set_mean_tau2(const number val) {m_mean_tau2=val;}
void Exp2Synapses::set_dev_tau2(const number val) {m_dev_tau2=val;}
void Exp2Synapses::set_mean_rev(const number val) {m_mean_rev=val;}
void Exp2Synapses::set_dev_rev(const number val) {m_dev_rev=val;}
void Exp2Synapses::set_mean_threshold(const number val) {m_mean_threshold=val;}
void Exp2Synapses::set_dev_threshold(const number val) {m_dev_threshold=val;}

size_t Exp2Synapses::size() {return m_vSynapses.size();}
size_t Exp2Synapses::start_id() {return m_start_id;}


std::vector<IBaseSynapse*> Exp2Synapses::get_synapses()
{
    boost::mt19937 onset_rng;
    boost::normal_distribution<number> onset_dist(m_mean_onset, m_dev_onset);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > onset_var(onset_rng, onset_dist);

    boost::mt19937 gMax_rng;
    boost::normal_distribution<number> gMax_dist(m_mean_gMax, m_dev_gMax);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > gMax_var(gMax_rng, gMax_dist);

    boost::mt19937 tau1_rng;
    boost::normal_distribution<number> tau1_dist(m_mean_tau1, m_dev_tau1);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau1_var(tau1_rng, tau1_dist);

    boost::mt19937 tau2_rng;
    boost::normal_distribution<number> tau2_dist(m_mean_tau2, m_dev_tau2);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau2_var(tau2_rng, tau2_dist);

    boost::mt19937 rev_rng;
    boost::normal_distribution<number> rev_dist(m_mean_rev, m_dev_rev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > rev_var(rev_rng, rev_dist);

    boost::mt19937 threshold_rng;
    boost::normal_distribution<number> threshold_dist(m_mean_threshold, m_dev_threshold);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > threshold_var(threshold_rng, threshold_dist);

    for(size_t i = 0; i < m_vSynapses.size(); i=i+2) {
        number tau1 = tau1_var();
        number tau2 = tau2_var();

        number tp = (tau1*tau2)/(tau2-tau1) * std::log(tau2/tau1);
        tp = tp + 3 * tau2;
        size_t syn_id = m_start_id + i/2;
        IBaseSynapse *s1 = new Exp2PostSynapse(syn_id, 0.0, gMax_var(), tau1, tau2, rev_var());
        IBaseSynapse *s2 = new ThresholdPreSynapse(syn_id, 0.0, tp, threshold_var());

        m_vSynapses[i] = s1;
        m_vSynapses[i+1] = s2;
    }

    return m_vSynapses;
}




AlphaSynapsePair::AlphaSynapsePair()
: m_pre(new OnsetPreSynapse(0, 0.0, 0.0, 2.4e-3)),
  m_post(new AlphaPostSynapse(0, 0.0, 1.2e-9, 4e-4, 0.0))
{}

AlphaSynapsePair::~AlphaSynapsePair()
{
    delete m_pre;
    delete m_post;
}

void AlphaSynapsePair::set_id(size_t id) {m_pre->set_id(id); m_post->set_id(id);}
void AlphaSynapsePair::set_onset(number onset) {m_pre->set_onset(onset);}
void AlphaSynapsePair::set_tau(number tau) {m_post->set_tau(tau); m_pre->set_duration(6*tau);}
void AlphaSynapsePair::set_gMax(number gMax) {m_post->set_gMax(gMax);}
void AlphaSynapsePair::set_reversal_potential(number ve) {m_post->set_rev(ve);}

OnsetPreSynapse* AlphaSynapsePair::pre_synapse() {return m_pre;}
AlphaPostSynapse* AlphaSynapsePair::post_synapse() {return m_post;}






Exp2SynapsePair::Exp2SynapsePair()
: m_pre(new ThresholdPreSynapse(0, 0.0, 5.585e-3, 0.0)),
  m_post(new Exp2PostSynapse(0, 0.0, 1.2e-9, 2e-4, 1.7e-3, 0.0))
{}

Exp2SynapsePair::~Exp2SynapsePair()
{
    delete m_pre;
    delete m_post;
}

void Exp2SynapsePair::set_id(size_t id) {m_pre->set_id(id); m_post->set_id(id);}
void Exp2SynapsePair::set_threshold(number t) {m_pre->set_threshold(t);}
void Exp2SynapsePair::set_gMax(number gMax) {m_post->set_gMax(gMax);}
void Exp2SynapsePair::set_reversal_potential(number ve) {m_post->set_rev(ve);}
void Exp2SynapsePair::set_taus(number tau1, number tau2)
{
    number tp = (tau1*tau2)/(tau2-tau1) * log(tau2/tau1);
    tp += 3 * tau2;

    m_post->set_tau1(tau1);
    m_post->set_tau2(tau2);
    m_pre->set_duration(tp);
}

ThresholdPreSynapse* Exp2SynapsePair::pre_synapse() {return m_pre;}
Exp2PostSynapse* Exp2SynapsePair::post_synapse() {return m_post;}



} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

