/*
 * SplitSynapseHandler.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: lreinhardt
 */

#include "split_synapse_handler.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

SplitSynapseHandler::SplitSynapseHandler()
:m_aSSyn(GlobalAttachments::attachment<AVSSynapse>("SplitSynapses")),
 m_bInited(false),
 m_spGrid(SPNULL)
{
	m_ssah.set_attachment(m_aSSyn);
}

number SplitSynapseHandler::current_on_edge(Edge* e, number t)
{
	std::vector<IBaseSynapse*> vSyns(m_aaSSyn[e]);
	number curr = 0.0;

	for(int i=0; i<vSyns.size(); ++i) {
		//postsynapses
		if(!vSyns[i]->split_type() ) {
			curr += static_cast<IPostSynapse*>(vSyns[i])->current(t);
		}
	}
	return curr;
}

bool SplitSynapseHandler::synapse_on_edge(const Edge* edge, size_t scv, number time, number& current)
{
	// check that handler has been init'ed
	UG_COND_THROW(!m_bInited, "Cannot provide synapse information before being init'ed.");

	// get edge synapse info
	std::vector<IBaseSynapse*> vSyns = m_aaSSyn[edge];

	bool active_synapse = false;
	current = 0.0;

	for (size_t i = 0; i < vSyns.size(); ++i)
	{
		IBaseSynapse* s = vSyns[i];

		if ((s->location() < 0.5 && scv == 0) || (s->location() >= 0.5 && scv == 1))
		{
			// get vmDisc potential values for edge
			number vm_postsyn = m_spCEDisc->vm(edge->vertex(scv));

			switch (STV::type(info))
			{
				// treat alpha synapse
				case ALPHA_SYNAPSE:
				{
					typedef synapse_traits<AlphaSynapse> STA;
					//if (SynapseSelector::is_declared(ALPHA_SYNAPSE))
					//{
						// create a default alpha synapse and populate with parameters
						AlphaSynapse alpha(1e-6*STA::g_max(info), 1e-3*STA::onset(info),	// TODO units
										   1e-3*STA::tau(info), 1e-3*STA::v_rev(info), vm_postsyn);
						current += alpha.current(time);
						/*if (STA::onset(info) <= time)
						{
							UG_LOG_ALL_PROCS("current current (alphasyn): " << current << std::endl);
							UG_LOG_ALL_PROCS("synInfo: " << info << std::endl);
						}*/
						active_synapse = true;
					//}
					break;
				}
				// treat bi-exponential synapse
				case EXP2_SYNAPSE:
				{
					typedef synapse_traits<Exp2Syn> STB;
					//if (SynapseSelector::is_declared(EXP2_SYNAPSE))
					//{
						// create a default exp2syn and populate with parameters
						Exp2Syn exp(1e-3*STB::tau1(info), 1e-3*STB::tau2(info),				// TODO units
									1e-3*STB::v_rev(info), 1e-6*STB::g_max(info), vm_postsyn);

						// get potential for presynapse
						unsigned int presynInd = STB::presyn_ind(info);
						UG_ASSERT(presynInd < m_vPresynVmValues.size(),
									"Requested presynaptic potential value for unknown index "
									<< presynInd << " (only " << m_vPresynVmValues.size()
									<< " presynaptic potential values are available)!");
						number vm_presyn = m_vPresynVmValues[presynInd];

						// activate synapse if V_m > threshold
						if (vm_presyn > 1e-3*STB::threshold(info) && !STB::activated(info))	// TODO units
							STB::activate(info, time);

						// handle synapse time
						if (STB::activated(info))
						{
							// active synapse
							//UG_LOG_ALL_PROCS("Active synapse " << STB::name() << "!" << std::endl);
							number syn_time = time - STB::onset(info);

							// deactivate if time past activity
							if (syn_time >= 1e-3*STB::activity_time(info))					// TODO units
								STB::deactivate(info);

							current += exp.current(syn_time);
							active_synapse = true;
						}
						//UG_LOG_ALL_PROCS("current current (exp2syn): " << current << std::endl);
					//}
					break;
				}

				case JANA_SYNAPSE_FROM_MARKUS_WITH_LOVE:
				{
					// create a default Jana synapse and populate with parameters
					typedef synapse_traits<JanaSynapseFromMarkusWithLove> STJ;
					JanaSynapseFromMarkusWithLove jsyn(1e-3*STJ::tau1(info), 1e-3*STJ::tau2(info),	// TODO units
							1e-3*STJ::onset(info), 1e-3*STJ::v_rev(info), 1e-6*STJ::g_max(info), vm_postsyn);

					current += jsyn.current(time);
					active_synapse = true;

					break;
				}

				// treat default synapse
				default:
					break;
			}
		}
	}

	return active_synapse;
}

void SplitSynapseHandler::grid_first_available()
{
	// throw if this is called more than once
	UG_COND_THROW(m_bInited, "Second initialization call is not allowed.");

	// Global Attachment setup

	// Check existence
	if (!GlobalAttachments::is_declared("SplitSynapses")) {
		UG_THROW("GlobalAttachment 'Synapses' not available.");
	}


	// Attach to grid, if not already present
	if (!m_spGrid->has_edge_attachment(m_aSSyn)) {
		m_spGrid->attach_to_edges(m_aSSyn);
	}

	// check that essential attachments exist in grid, create accessors
	UG_COND_THROW(!m_spGrid->has_attachment<Edge>(m_aSSyn), "No synapse info attached to grid!");

	m_aaSSyn = Grid::EdgeAttachmentAccessor<AVSSynapse>(*m_spGrid, m_aSSyn);


	// propagate attachments through levels
	m_ssah.set_grid(m_spGrid);

	// set init'ed flag
	m_bInited = true;
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
