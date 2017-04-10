/*
 * SynapseHandler.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: lreinhardt, mbreit
 */

#include "synapse_handler.h"

#include "common/error.h"
#ifdef UG_PARALLEL
#include "pcl/pcl_process_communicator.h"
#endif
#include "../util/functions.h"   // neuron_identification
#include "synapse_info_io_traits.h" // needed in GlobalAttachments::attachment
#include "synapse_attachment_serializer.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "synapses/onset_pre_synapse.h"
#include "synapses/threshold_pre_synapse.h"
#include "synapses/alpha_post_synapse.h"
#include "synapses/exp2_post_synapse.h"

#include <algorithm>    // std::sort


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

template <typename TDomain>
SynapseHandler<TDomain>::SynapseHandler()
: m_aSyn(GlobalAttachments::attachment<AVSynapse>("synapses")),
  m_bInited(false),
  m_spGrid(SPNULL),
  m_spCEDisc(SPNULL),
  m_spApprox(SPNULL),
  m_mPreSynapses(std::map<synapse_id,IPreSynapse*>()),
  m_mPostSynapses(std::map<synapse_id,IPostSynapse*>()),
  m_mActivePreSynapses(std::map<synapse_id, IPreSynapse*>()),
  m_mPreSynapseIdToEdge(std::map<synapse_id,Edge*>()),
  m_mPostSynapseIdToEdge(std::map<synapse_id,Edge*>()),
  m_spSAS(SPNULL),
  m_primary_alpha_onset_mean(std::numeric_limits<number>::max()), m_primary_alpha_onset_dev(0.0),
  m_primary_alpha_tau_mean(0.0), m_primary_alpha_tau_dev(0.0),
  m_primary_alpha_peak_cond_mean(6e-4), m_primary_alpha_peak_cond_dev(0.0),
  m_primary_alpha_constSeed(true),
  m_prim_biexp_onset_mean(std::numeric_limits<number>::max()), m_prim_biexp_onset_dev(0.0),
  m_prim_biexp_tau1_mean(0.0), m_prim_biexp_tau1_dev(0.0), m_prim_biexp_tau2_mean(0.0), m_prim_biexp_tau2_dev(0.0),
  m_prim_biexp_peak_cond_mean(0.0), m_prim_biexp_peak_cond_dev(0.0),
  m_prim_biexp_constSeed(true)
{
	m_ssah.set_attachment(m_aSyn);
}


template <typename TDomain>
number SynapseHandler<TDomain>::current_on_edge(const Edge* e, size_t scv, number t) const
{
	const std::vector<IBaseSynapse*>& vSyns = m_aaSyn[e];
	number vm_postsyn = m_spCEDisc->vm(e->vertex(scv));
	number curr = 0.0;

	for(size_t i=0; i<vSyns.size(); ++i) {
		IBaseSynapse* s = vSyns[i];

		// only post synapses can invoke a current
		if(s->is_postsynapse()) {
			IPostSynapse* ps = static_cast<IPostSynapse*>(s);

			if(ps->is_active(t)) {
				curr += ps->current(t, vm_postsyn);
			}
		}
	}

	return curr;
}

template <typename TDomain>
number SynapseHandler<TDomain>::current(synapse_id id) const
{
	//find edge

	ConstEdgeIterator eIter = m_spGrid->begin<Edge>();
	ConstEdgeIterator eIterEnd = m_spGrid->end<Edge>();
	while(eIter != eIterEnd) {
		Edge* e = *eIter;
		std::vector<IBaseSynapse*> v = m_aaSyn[e];
		for(size_t i=0; i<v.size(); ++i) { //find synapse
			if( (v[i]->id()==id) &&  (v[i]->is_postsynapse())) {


				number rel_loc=v[i]->location();
				number t = m_spCEDisc->time();
				Vertex* v0 = e->vertex(0);
				Vertex* v1 = e->vertex(1);
				number vm = m_spCEDisc->vm( v0)*rel_loc + m_spCEDisc->vm(v1)*(1-rel_loc);
				number curr = ((IPostSynapse*)v[i])->current(t, vm);

				return curr;
			}
		}
		eIter++;
	}
	return 0;
}


template <typename TDomain>
bool SynapseHandler<TDomain>::synapse_on_edge(const Edge* edge, size_t scv, number time, number& current) const
{
	current = current_on_edge(edge, scv,time);

	return true;
}

template <typename TDomain>
void SynapseHandler<TDomain>::grid_first_available()
{
	// throw if this is called more than once
	UG_COND_THROW(m_bInited, "Second initialization call is not allowed.");

	// check availability of approxSpace and grid; set members
	UG_COND_THROW(!m_spCEDisc.valid(), "Given CableEquation SmartPtr is not valid.");

	m_spApprox = m_spCEDisc->approx_space();
	UG_COND_THROW(!m_spApprox.valid(), "No valid approximation space available in synapse handler.\n"
				  "Did you forget to set a CableEquation via set_ce_object()?");

	UG_COND_THROW(!m_spApprox->domain().valid(), "The approximation space of the given CableEquation object"
				  " does not contain a valid domain.\n");
	m_spGrid = m_spApprox->domain()->grid();
	UG_COND_THROW(!m_spGrid.valid(), "There is no grid associated to the CableEquation object passed.\n"
				  "Make sure you load the domain before setting the CableEquation.");

    // check that essential attachments exist in grid
    if (!m_spGrid->has_attachment<Edge>(m_aSyn))
    {
        // the synapse attachment is allowed to be missing, but only if the grid is empty
        // (which may happen in parallel, if the local proc has no grid)
        UG_COND_THROW(m_spGrid->num_vertices(), "No synapses attached to grid!");

        // formally attach synapse attachment
        m_spGrid->attach_to_edges(m_aSyn);
    }

	// position access
    m_aaPosition = Grid::VertexAttachmentAccessor<APosition>(*m_spGrid, aPosition);

    // synapse access
    m_aaSyn = Grid::EdgeAttachmentAccessor<AVSynapse>(*m_spGrid, m_aSyn);

    // propagate synapse attachments through levels
    m_ssah.set_grid(m_spGrid);

    // create a serialization object for distribution purposes
    m_spSAS = make_sp(new SynapseAttachmentSerializer(*m_spGrid, m_aSyn));

   	// gather all synapses from surface and build all maps
	collect_synapses_from_grid();

	// set activation timing for alpha synapses
    set_activation_timing_with_grid();

    // set all post-synapses to deactivated status (might not be the case)
    deactivate_all_postSyns();

    // set this class as event listener for distribution
    m_spGridDistributionCallbackID = m_spGrid->message_hub()->register_class_callback(this,
        &SynapseHandler<TDomain>::grid_distribution_callback);

    // TODO: We also need a refinement callback
    //       The pointers to synapses being held by this class need to be on the surface!

    // set init'ed flag
    m_bInited = true;
}

#if 0
template <typename TDomain>
void
SynapseHandler<TDomain>::get_currents
(
    number t,
    int nid,
    std::vector<number>& vCurrOut,
    std::vector<synapse_id>& vSidOut
)
{
	vCurrOut.clear();
	vSidOut.clear();

	// gather local synapse ids and their current currents
#ifdef UG_PARALLEL
    if (pcl::NumProcs() > 1)
    {
        std::vector<number> vlocalcurr;
        std::vector<synapse_id> vlocalsid;
        for (size_t i=0; i<m_vAllSynapses.size(); ++i) {
            if (m_vAllSynapses[i]->is_postsynapse() /*&& m_aaNID[m_vAllSynapses[i]] == nid*/) {
                IPostSynapse* s = (IPostSynapse*)m_vAllSynapses[i];

                number vm; // FIXME: get membrane potential at post_synapse from cable_equation!
                vlocalcurr.push_back(s->current(t, vm));
                vlocalsid.push_back(s->id());
            }
        }

        pcl::ProcessCommunicator com;
        com.allgatherv(vCurrOut, vlocalcurr);
        com.allgatherv(vSidOut, vlocalsid);
	}
    else
    {
#endif
        for (size_t i=0; i<m_vAllSynapses.size(); ++i) {
            if (m_vAllSynapses[i]->is_postsynapse() /*&& m_aaNID[m_vAllSynapses[i]] == nid*/) {
                IPostSynapse* s = (IPostSynapse*)m_vAllSynapses[i];

                number vm; // FIXME: get membrane potential at post_synapse from cable_equation!
                vCurrOut.push_back(s->current(t, vm));
                vSidOut.push_back(s->id());
            }
        }
#ifdef UG_PARALLEL
    }
#endif
}
#endif

template <typename TDomain>
void
SynapseHandler<TDomain>::collect_synapses_from_grid()
{
	// clear synapse lists
	m_vAllSynapses.clear();
	m_mPreSynapses.clear();
	m_mPostSynapses.clear();
	m_mPreSynapseIdToEdge.clear();
	m_mPostSynapseIdToEdge.clear();

	// (re)populate synapse lists from surface(!) edges
	ConstSmartPtr<DoFDistribution> dd = m_spApprox->dd(GridLevel());
	typename DoFDistribution::traits<Edge>::const_iterator eit = dd->begin<Edge>();
	typename DoFDistribution::traits<Edge>::const_iterator eit_end = dd->end<Edge>();
	for (; eit != eit_end; ++eit)
	{
		const std::vector<IBaseSynapse*>& syns = m_aaSyn[*eit];
		size_t sz = syns.size();
		for (size_t i = 0; i < sz; ++i)
		{
			IBaseSynapse* syn = syns[i];
			m_vAllSynapses.push_back(syn);
			if (syn->is_presynapse())
			{
				// presynapse map and save edge
				IPreSynapse* s = (IPreSynapse*) syn;
				m_mPreSynapses[s->id()] = s;
				m_mPreSynapseIdToEdge[s->id()] = *eit;
			}
			else
			{
				// postsynapse map
				IPostSynapse* s = (IPostSynapse*) syn;
				m_mPostSynapses[s->id()] = s;
				m_mPostSynapseIdToEdge[s->id()] = *eit;
			}
		}
	}

	std::sort(m_vAllSynapses.begin(), m_vAllSynapses.end(), __comp());

	// TODO: save start and end indices for particular types
	// this will allow const access via begin() and end()
}


/*template <typename TDomain>
void SynapseHandler<TDomain>::
set_ce_object(SmartPtr<CableEquation<TDomain> > disc)
{
	UG_COND_THROW(m_bInited, "The CableEquation object associated to this synapse handler "
				  "must not be changed\nafter addition of the original CableEquation object "
				  "to the domain discretization.");

	// set cable equation disc object
	m_spCEDisc = disc;
}*/



template <typename TDomain>
void SynapseHandler<TDomain>::
update_presyn(number time)
{
#ifdef UG_PARALLEL
    // collect postsynapse id's that became active/inactive
	std::vector<synapse_id> vNewActivePostSynapseIds_local; // becoming active on local process
	std::vector<synapse_id> vNewInactivePostSynapseIds_local; // becoming inactive on local process
#endif

	std::map<synapse_id, IPreSynapse*>::iterator it = m_mPreSynapses.begin();
	for (; it != m_mPreSynapses.end(); ++it) {
		IPreSynapse* s = it->second;
		Edge* e = m_mPreSynapseIdToEdge[it->first] ; // edge on which s is located

		Vertex* v1 = e->vertex(0);
		Vertex* v2 = e->vertex(1);

		// TODO: synapse might depend on other values from CE besides membrane potential
		std::vector<number> uAtSynapseLocation(1);
		uAtSynapseLocation[0] = s->location() * m_spCEDisc->vm(v2)
								+ (1.0 - s->location()) *  m_spCEDisc->vm(v1);

		s->update(time, uAtSynapseLocation);

		// pre-synapse becomes active if it is active and cannot be found in the active presynapses map
		if (s->is_active(time) && m_mActivePreSynapses.find(s->id()) == m_mActivePreSynapses.end()) {
			m_mActivePreSynapses[s->id()] = s;
#ifdef UG_PARALLEL
			if (pcl::NumProcs() > 1)
			    vNewActivePostSynapseIds_local.push_back(s->id());
			else
#endif
			    m_mPostSynapses[s->id()]->activate(time);

		// pre-synapse becomes inactive
		} else if (!s->is_active(time) && m_mActivePreSynapses.find(s->id()) != m_mActivePreSynapses.end()) {
		    m_mActivePreSynapses.erase(s->id());
#ifdef UG_PARALLEL
		    if (pcl::NumProcs() > 1)
		        vNewInactivePostSynapseIds_local.push_back(s->id());
		    else
#endif
		        m_mPostSynapses[s->id()]->deactivate();
		}
	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() == 1) return;

	// TODO: improve this; use pre-computed information on which process has the post synapse
	//       and send only to that process (distribute instead of allgatherv)
	// propagate to all processes
	pcl::ProcessCommunicator com;

    std::vector<synapse_id> vNewActivePostSynapseIds_global; // complete list of all newly active synapses
    std::vector<synapse_id> vNewInactivePostSynapseIds_global; // complete list of all newly inactive synapses

    com.allgatherv(vNewActivePostSynapseIds_global, vNewActivePostSynapseIds_local);
	com.allgatherv(vNewInactivePostSynapseIds_global, vNewInactivePostSynapseIds_local);

	// scan for synapse ids that are on the local process and have to be activated
	for (size_t i=0; i<vNewActivePostSynapseIds_global.size(); ++i) {
		synapse_id psid = vNewActivePostSynapseIds_global[i];

		std::map<synapse_id, IPostSynapse*>::iterator post_it = m_mPostSynapses.find(psid);
		if (post_it != m_mPostSynapses.end()) { // synapse id found on local process
			post_it->second->activate(time);
			m_mActivePostSynapses[psid] = post_it->second;
		}
	}

	// scan for synapse ids that are on the local process and have to be deactivated
	for (size_t i=0; i<vNewInactivePostSynapseIds_global.size(); ++i) {
		synapse_id psid = vNewInactivePostSynapseIds_global[i];

		std::map<synapse_id, IPostSynapse*>::iterator post_it = m_mPostSynapses.find(psid);
		if (post_it != m_mPostSynapses.end()) {
			post_it->second->deactivate();
			m_mActivePostSynapses.erase(psid);
		}
	}
#endif
}


template <typename TDomain>
void SynapseHandler<TDomain>::
deactivate_all_postSyns()
{
    std::map<synapse_id, IPostSynapse*>::iterator it = m_mPostSynapses.begin();
    std::map<synapse_id, IPostSynapse*>::iterator it_end = m_mPostSynapses.end();
    for (; it != it_end; ++it)
       it->second->deactivate();
}


template <typename TDomain>
void SynapseHandler<TDomain>::
set_activation_timing_alpha
(
    number onset,
    number tau,
    number peak_cond,
    number onset_dev,
    number tau_dev,
    number peak_cond_dev,
    bool constSeed
)
{
    UG_COND_THROW(m_bInited, "The activation timing cannot be changed after addition of the\n"
                  "original CableEquation object to the domain discretization.");

    m_primary_alpha_onset_mean = onset;
    m_primary_alpha_tau_mean = tau;
    m_primary_alpha_peak_cond_mean = peak_cond;
    m_primary_alpha_onset_dev = onset_dev;
    m_primary_alpha_tau_dev = tau_dev;
    m_primary_alpha_peak_cond_dev = peak_cond_dev;
    m_primary_alpha_constSeed = constSeed;
}


template <typename TDomain>
void SynapseHandler<TDomain>::
set_activation_timing_biexp
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
)
{
    UG_COND_THROW(m_bInited, "The activation timing cannot be changed after addition of the\n"
                  "original CableEquation object to the domain discretization.");

    m_prim_biexp_onset_mean = onset_mean;
    m_prim_biexp_tau1_mean = tau1_mean;
    m_prim_biexp_tau2_mean = tau2_mean;
    m_prim_biexp_peak_cond_mean = peak_cond_mean;
    m_prim_biexp_onset_dev = onset_dev;
    m_prim_biexp_tau1_dev = tau1_dev;
    m_prim_biexp_tau2_dev = tau2_dev;
    m_prim_biexp_peak_cond_dev = peak_cond_dev;
    m_prim_biexp_constSeed = constSeed;
}


template <typename TDomain>
void SynapseHandler<TDomain>::
add_activation_timing_alpha_ball
(
    const std::vector<number>& alpha_timings,
    const std::vector<number>& ball
)
{
    // consistency checks
    UG_COND_THROW(m_bInited, "The activation timing cannot be changed after addition of the\n"
                  "original CableEquation object to the domain discretization.");

    UG_COND_THROW(alpha_timings.size() != 6, "Expected 6 timing values for alpha synapses.");
    UG_COND_THROW(ball.size() != 4, "Expected 4 parameters (x, y, z, d) to describe the ball in 3d.");

    // add activation timing for a ball
    m_vTimingAlphaBalls.push_back(std::make_pair(alpha_timings, ball));
}

template <typename TDomain>
void SynapseHandler<TDomain>::
add_activation_timing_exp2_ball
(
    const std::vector<number>& biexp_timings,
    const std::vector<number>& ball
)
{
    // consistency checks
    UG_COND_THROW(m_bInited, "The activation timing cannot be changed after addition of the\n"
                  "original CableEquation object to the domain discretization.");

    UG_COND_THROW(biexp_timings.size() != 8, "Expected 8 timing values for exp2 synapses.");
    UG_COND_THROW(ball.size() != 4, "Expected 4 parameters (x, y, z, d) to describe the ball in 3d.");

    // add activation timing for a ball
    m_vTimingBiExpBalls.push_back(std::make_pair(biexp_timings, ball));
}

template <typename TDomain>
void SynapseHandler<TDomain>::
set_activation_timing_with_grid()
{
    // check availability of all needed structures
    UG_COND_THROW(!m_spGrid.valid(), "No valid grid. Make sure that the synapse handler has a grid to work on!");

    // activity timing setup: random normal distribution
    boost::mt19937 rng_alpha;
    boost::mt19937 rng_biexp;
#ifdef UG_PARALLEL
    if (m_primary_alpha_constSeed)
        // what to use here?
        // constant -- all procs would have the same pattern and outcome depends on distribution
        // procRank() -- outcome would still depend on distribution
        // neuronID -- seems to be a viable option
        rng_alpha.seed(pcl::ProcRank()); // use neuron id instead of 0 so that patterns are not identical for all procs
    else
        rng_alpha.seed(pcl::ProcRank()*time(NULL));

    if (m_prim_biexp_constSeed)
        rng_biexp.seed(pcl::ProcRank());
    else
        rng_biexp.seed(pcl::ProcRank()*time(NULL));

#else
    if (m_primary_alpha_constSeed)
        rng_alpha.seed(0);
    else
        rng_alpha.seed(time(NULL));

    if (m_prim_biexp_constSeed)
        rng_biexp.seed(0);
    else
        rng_biexp.seed(time(NULL));
#endif

    // alpha synapse distributions
    boost::normal_distribution<number> start_dist(m_primary_alpha_onset_mean, m_primary_alpha_onset_dev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_start(rng_alpha, start_dist);
    boost::normal_distribution<number> tau_dist(m_primary_alpha_tau_mean, m_primary_alpha_tau_dev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_tau(rng_alpha, tau_dist);
    boost::normal_distribution<number> cond_dist(m_primary_alpha_peak_cond_mean, m_primary_alpha_peak_cond_dev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_cond(rng_alpha, cond_dist);

    // biexp synapse distributions
    boost::normal_distribution<number> onset_dist(m_prim_biexp_onset_mean, m_prim_biexp_onset_dev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_onset(rng_biexp, onset_dist);
    boost::normal_distribution<number> tau1_dist(m_prim_biexp_tau1_mean, m_prim_biexp_tau1_dev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_tau1(rng_biexp, tau1_dist);
    boost::normal_distribution<number> tau2_dist(m_prim_biexp_tau2_mean, m_prim_biexp_tau2_dev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_tau2(rng_biexp, tau2_dist);
    boost::normal_distribution<number> peakCond_dist(m_prim_biexp_peak_cond_mean, m_prim_biexp_peak_cond_dev);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_peakCond(rng_biexp, peakCond_dist);

    // loop onset_pre_synapses
    SynapseIter<OnsetPreSynapse> it = begin<OnsetPreSynapse>();
    SynapseIter<OnsetPreSynapse> it_end = end<OnsetPreSynapse>();
    for (; it != it_end; ++it)
    {
        OnsetPreSynapse* onsetSyn = *it;

        // get corresponding post-synapse
        std::map<synapse_id, IPostSynapse*>::iterator postIt = m_mPostSynapses.find(onsetSyn->id());

        // Technically, it is possible that this post synapse is not located on the same proc.
        // This should not happen, as OnsetPreSynapses belong to primary post-synapses and those are located
        // on the same edge, but you never know what an evil user might do ...
        UG_COND_THROW(postIt == m_mPostSynapses.end(), "Post-synapse for OnsetPreSynapse not on the same proc.")
        IPostSynapse* post = postIt->second;

        // depending on type, set activation
        SynapseType type = post->type();
        switch (type)
        {
            case ALPHA_POST_SYNAPSE:
            {
                // if parameters have not been set, do nothing
                if (m_primary_alpha_onset_mean == std::numeric_limits<number>::max())
                    break;

                // cast to alpha post-synapse
                AlphaPostSynapse* alphaPost = dynamic_cast<AlphaPostSynapse*>(post);
                UG_COND_THROW(!alphaPost, "Synapse claiming to be alpha post synapse, but is not.");

                number t_onset = var_start();
                t_onset = t_onset < 0 ? 0 : t_onset;

                number tau = var_tau();
                tau = tau < 0 ? 0 : tau;

                number cond = var_cond();
                cond = cond < 0 ? 0 : cond;

                // This will parameterize the alpha_synapse in such a way that
                // it is active until the current decays to 5% of its maximal strength.
                onsetSyn->set_onset(t_onset);
                onsetSyn->set_duration(6.0*tau);
                alphaPost->set_tau(tau);
                alphaPost->set_gMax(cond);

                break;
            }

            case EXP2_POST_SYNAPSE:
            {
                // if parameters have not been set, do nothing
                if (m_prim_biexp_onset_mean == std::numeric_limits<number>::max())
                    break;

                // cast to exp2 post-synapse
                Exp2PostSynapse* exp2Post = dynamic_cast<Exp2PostSynapse*>(post);
                UG_COND_THROW(!exp2Post, "Synapse claiming to be exp2 post synapse, but is not.");

                number t_onset = var_onset();
                t_onset = t_onset < 0 ? 0 : t_onset;

                number tau1 = var_tau1();
                tau1 = tau1 < 0 ? 0 : tau1;

                number tau2 = var_tau2();
                tau2 = tau2 < 0 ? 0 : tau2;

                number peakCond = var_peakCond();
                peakCond = peakCond < 0 ? 0 : peakCond;

                number duration = tau1 * tau2 / (tau2 - tau1) * log(tau2 / tau1);
                duration += 3*tau2;


                onsetSyn->set_onset(t_onset);
                onsetSyn->set_duration(duration);
                exp2Post->set_tau1(tau1);
                exp2Post->set_tau2(tau2);
                exp2Post->set_gMax(peakCond);

                break;
            }

            default: break;
        }
    }


    // finally, evaluate alpha synapse activation balls
    // possibly overriding the previously set activation pattern
    typedef std::vector<std::pair<std::vector<number>, std::vector<number> > >::const_iterator ball_it_type;
    ball_it_type ball_it = m_vTimingAlphaBalls.begin();
    ball_it_type ball_it_end = m_vTimingAlphaBalls.end();
    for (; ball_it != ball_it_end; ++ball_it)
    {
        const std::vector<number>& alpha_timings = ball_it->first;
        const std::vector<number>& ball = ball_it->second;

        // timings
        number onset = alpha_timings[0];
        number onset_dev = alpha_timings[1];
        number tau = alpha_timings[2];
        number tau_dev = alpha_timings[3];
        number peak_cond = alpha_timings[4];
        number peak_cond_dev = alpha_timings[5];

        // ball region
        vector3 center(ball[0], ball[1], ball[2]);
        number radius = ball[3];

        // sample from distribution for this ball
        boost::normal_distribution<number> start_dist(onset, onset_dev);
        boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_start(rng_alpha, start_dist);
        boost::normal_distribution<number> tau_dist(tau, tau_dev);
        boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_tau(rng_alpha, tau_dist);
        boost::normal_distribution<number> cond_dist(peak_cond, peak_cond_dev);
        boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_cond(rng_alpha, cond_dist);

        SynapseIter<OnsetPreSynapse> it = begin<OnsetPreSynapse>();
        for (; it != it_end; ++it)
        {
            OnsetPreSynapse* onsetSyn = *it;
            // check whether synapse is located in within ball
            Edge* e = m_mPreSynapseIdToEdge[onsetSyn->id()];
            vector3 a = m_aaPosition[e->vertex(0)];
            vector3 b = m_aaPosition[e->vertex(1)];
            /// FIXME: unsafe comparison
            if (VecDistanceSq(a, center) >= radius*radius || VecDistanceSq(b, center) >= radius*radius)
                continue;

            // here, this is the case; now, get corresponding post-synapse
            std::map<synapse_id, IPostSynapse*>::iterator postIt = m_mPostSynapses.find(onsetSyn->id());

            // Technically, it is possible that this post synapse is not located on the same proc.
            // This should not happen, as OnsetPreSynapses belong to primary post-synapses and those are located
            // on the same edge, but you never know what an evil user might do ...
            UG_COND_THROW(postIt == m_mPostSynapses.end(), "Post-synapse for OnsetPreSynapse not on the same proc.")
            IPostSynapse* post = postIt->second;

            // depending on type, set activation
            SynapseType type = post->type();
            switch (type)
            {
                case ALPHA_POST_SYNAPSE:
                {
                    // cast to alpha post-synapse
                    AlphaPostSynapse* alphaPost = dynamic_cast<AlphaPostSynapse*>(post);
                    UG_COND_THROW(!alphaPost, "Synapse claiming to be alpha post synapse, but is not.");

                    number t_onset = var_start();
                    t_onset = t_onset < 0 ? 0 : t_onset;

                    number tau = var_tau();
                    tau = tau < 0 ? 0 : tau;

                    number cond = var_cond();
                    cond = cond < 0 ? 0 : cond;

                    // This will parameterize the alpha_synapse in such a way that
                    // it is active until the current decays to 5% of its maximal strength.
                    onsetSyn->set_onset(t_onset);
                    onsetSyn->set_duration(6.0*tau);
                    alphaPost->set_tau(tau);
                    alphaPost->set_gMax(cond);

                    break;
                }

                // exp2 balls do not exist (yet)

                default:
                    break;
            }
        }
    }

    // finally, evaluate biexp synapse activation balls
    // possibly overriding the previously set activation pattern
    typedef std::vector<std::pair<std::vector<number>, std::vector<number> > >::const_iterator ball_it_type;
    ball_it = m_vTimingBiExpBalls.begin();
    ball_it_end = m_vTimingBiExpBalls.end();
    for (; ball_it != ball_it_end; ++ball_it)
    {
        const std::vector<number>& biexp_timings = ball_it->first;
        const std::vector<number>& ball = ball_it->second;

        // timings
        number onset = biexp_timings[0];
        number tau1_mean = biexp_timings[1];
        number tau2_mean = biexp_timings[2];
        number peak_cond = biexp_timings[3];
        number onset_dev = biexp_timings[4];
        number tau1_dev = biexp_timings[5];
        number tau2_dev = biexp_timings[6];
        number peak_cond_dev = biexp_timings[7];

        // ball region
        vector3 center(ball[0], ball[1], ball[2]);
        number radius = ball[3];

        // sample from distribution for this ball
        boost::normal_distribution<number> start_dist(onset, onset_dev);
        boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_start(rng_biexp, start_dist);
        boost::normal_distribution<number> tau1_dist(tau1_mean, tau1_dev);
        boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_tau1(rng_biexp, tau1_dist);
        boost::normal_distribution<number> tau2_dist(tau2_mean, tau2_dev);
        boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_tau2(rng_biexp, tau2_dist);
        boost::normal_distribution<number> cond_dist(peak_cond, peak_cond_dev);
        boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > var_cond(rng_biexp, cond_dist);

        SynapseIter<OnsetPreSynapse> it = begin<OnsetPreSynapse>();
        for (; it != it_end; ++it)
        {
            OnsetPreSynapse* onsetSyn = *it;

            // check whether synapse is located in within ball
            Edge* e = m_mPreSynapseIdToEdge[onsetSyn->id()];
            vector3 a = m_aaPosition[e->vertex(0)];
            vector3 b = m_aaPosition[e->vertex(1)];
            /// FIXME: unsafe comparison
            if (VecDistanceSq(a, center) >= radius*radius || VecDistanceSq(b, center) >= radius*radius)
                continue;

            // here, this is the case; now, get corresponding post-synapse
            std::map<synapse_id, IPostSynapse*>::iterator postIt = m_mPostSynapses.find(onsetSyn->id());

            // Technically, it is possible that this post synapse is not located on the same proc.
            // This should not happen, as OnsetPreSynapses belong to primary post-synapses and those are located
            // on the same edge, but you never know what an evil user might do ...
            UG_COND_THROW(postIt == m_mPostSynapses.end(), "Post-synapse for OnsetPreSynapse not on the same proc.")
            IPostSynapse* post = postIt->second;

            // depending on type, set activation
            SynapseType type = post->type();
            switch (type)
            {
            case EXP2_POST_SYNAPSE:
                {
                    // cast to alpha post-synapse
                	Exp2PostSynapse* exp2Post = dynamic_cast<Exp2PostSynapse*>(post);
                    UG_COND_THROW(!exp2Post, "Synapse claiming to be exp2 post synapse, but is not.");

                    number t_onset = var_start();
                    t_onset = t_onset < 0 ? 0 : t_onset;

                    number tau1 = var_tau1();
                    tau1 = tau1 < 0 ? 0 : tau1;

                    number tau2 = var_tau2();
                    tau2 = tau2 < 0 ? 0 : tau2;

                    number cond = var_cond();
                    cond = cond < 0 ? 0 : cond;

                    number duration = tau1 * tau2 / (tau2 - tau1) * log(tau2 / tau1);
                    duration += 3*tau2;

                    onsetSyn->set_onset(t_onset);
                    onsetSyn->set_duration(duration);
                    exp2Post->set_tau1(tau1);
                    exp2Post->set_tau2(tau2);
                    exp2Post->set_gMax(cond);
                    break;
                }
                default:
                    break;
            	}
    	}
    }
}



template <typename TDomain>
void SynapseHandler<TDomain>::show_status(number t)
{
    for(size_t i=0; i<m_vAllSynapses.size(); ++i) {
    	if(m_vAllSynapses[i]->is_presynapse()) {
    		IPreSynapse* s = (IPreSynapse*)m_vAllSynapses[i];
    		std::cout << i << ":" << t << " " << m_vAllSynapses[i] << ( (s->is_active(t))?" active":" inactive")  << std::endl;
    	} else if(m_vAllSynapses[i]->is_postsynapse()) {
    		IPostSynapse* s = (IPostSynapse*)m_vAllSynapses[i];
    		std::cout << i << ":" << t << " " << m_vAllSynapses[i] << " " << current(s->id()) << std::endl;
    	}
    }
    std::cout << std::endl;
    //std::cout << "Synapse iterator begin(): " << &*(m_vAllSynapses.begin()) << std::endl;
    //std::cout << "Synapse iterator end(): " << &*(m_vAllSynapses.end()) << std::endl;
    //std::cout << "AlphaPreSynapse iterator end(): " << &*end<AlphaPreSynapse>() << std::endl;
}

#if 0
// unused
template <typename TDomain>
const std::vector<synapse_id> SynapseHandler<TDomain>::active_presynapses() const
{
	std::map<synapse_id, IPreSynapse*>::const_iterator it = m_mActivePreSynapses.begin();
	std::map<synapse_id, IPreSynapse*>::const_iterator itEnd = m_mActivePreSynapses.end();

	std::vector<synapse_id> syn_id_ret;

	while(it != itEnd) {
		syn_id_ret.push_back(it->first);
		++it;
	}

	return syn_id_ret;
}
#endif

template <typename TDomain>
void SynapseHandler<TDomain>::active_postsynapses_and_currents
(
	std::vector<synapse_id>& vActSynOut,
	std::vector<number>& vSynCurrOut,
	const std::vector<uint>& vNid,
	const Grid::VertexAttachmentAccessor<uint>& aaNID,
	number time
) const
{
	vActSynOut.clear();
	vSynCurrOut.clear();

	const size_t nNid = vNid.size();

	std::map<synapse_id, IPostSynapse*>::const_iterator it = m_mActivePostSynapses.begin();
	std::map<synapse_id, IPostSynapse*>::const_iterator itEnd = m_mActivePostSynapses.end();

	for (; it != itEnd; ++it)
	{
		synapse_id id = it->first;
		Edge* e = m_mPostSynapseIdToEdge.at(id);
		Vertex* v0 = e->vertex(0);

		// only treat neuron IDs of interest further
		uint nid = aaNID[v0];
		for (size_t n = 0; n < nNid; ++n)
			if (vNid[n] == nid)
				goto nid_of_interest;
		continue;

nid_of_interest:
		Vertex* v1 = e->vertex(1);

		number rel_loc = it->second->location();
		number vm = m_spCEDisc->vm(v0)*(1.0-rel_loc) + m_spCEDisc->vm(v1)*rel_loc;

		IPostSynapse* post = m_mPostSynapses.at(id);
		number curr = post->current(time, vm);

		vActSynOut.push_back(id);
		vSynCurrOut.push_back(curr);
	}
}


template <typename TDomain>
Edge* SynapseHandler<TDomain>::postsyn_edge(synapse_id postSynID) const
{
	std::map<synapse_id, Edge*>::const_iterator it = m_mPostSynapseIdToEdge.find(postSynID);
	UG_COND_THROW(it == m_mPostSynapseIdToEdge.end(), "Requested edge for post-synapse with ID "
		<< postSynID << " which is not present.");

	return it->second;
}

template <typename TDomain>
Edge* SynapseHandler<TDomain>::presyn_edge(synapse_id preSynID) const
{
	std::map<synapse_id, Edge*>::const_iterator it = m_mPreSynapseIdToEdge.find(preSynID);
	UG_COND_THROW(it == m_mPreSynapseIdToEdge.end(), "Requested edge for pre-synapse with ID "
		<< preSynID << " which is not present.");

	return it->second;
}


template <typename TDomain>
IPreSynapse* SynapseHandler<TDomain>::pre_synapse(synapse_id id) const
{
	return m_mPreSynapses.at(id);
}


template <typename TDomain>
IPostSynapse* SynapseHandler<TDomain>::post_synapse(synapse_id id) const
{
	return m_mPostSynapses.at(id);
}


template <typename TDomain>
void SynapseHandler<TDomain>::grid_distribution_callback(const GridMessage_Distribution& gmd)
{
	if (gmd.msg() == GMDT_DISTRIBUTION_STARTS)
	{
		GridDataSerializationHandler& gdsh = gmd.serialization_handler();
		gdsh.add(m_spSAS);
	}
	else if (gmd.msg() == GMDT_DISTRIBUTION_STOPS)
	{
		// gather all synapses from grid on the proc and build all tables
		collect_synapses_from_grid();
	}
}

// ////////////////////////////////////
//	explicit template instantiations //
// ////////////////////////////////////

#ifdef UG_DIM_1
	template class SynapseHandler<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class SynapseHandler<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class SynapseHandler<Domain3d>;
#endif

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
