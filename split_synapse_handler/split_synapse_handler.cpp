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

template <typename TDomain>
SplitSynapseHandler<TDomain>::SplitSynapseHandler()
:m_aSSyn(GlobalAttachments::attachment<AVSSynapse>("SplitSynapses")),
 m_bInited(false),
 m_spGrid(SPNULL),
 m_spCEDisc(SPNULL),
 m_spApprox(SPNULL),
 m_vAllSynapses(std::vector<IBaseSynapse*>()),
 m_vPreSynapses(std::vector<IPreSynapse*>()),
 m_vPostSynapses(std::vector<IPostSynapse*>()),
 m_mPostSynapses(std::map<unsigned long long,IPostSynapse*>()),
 m_mActivePreSynapses(std::map<unsigned long long,IPreSynapse*>())
{
	m_ssah.set_attachment(m_aSSyn);
}

template <typename TDomain>
number SplitSynapseHandler<TDomain>::current_on_edge(const Edge* e, size_t scv, number t)
{
	std::vector<IBaseSynapse*> vSyns(m_aaSSyn[e]);
	number vm_postsyn = m_spCEDisc->vm(e->vertex(scv));
	number curr = 0.0;

	for(int i=0; i<vSyns.size(); ++i) {
		//postsynapses
		if(!vSyns[i]->is_presynapse()) {
			curr += static_cast<IPostSynapse*>(vSyns[i])->current(t, vm_postsyn);
		}
	}
	return curr;
}

template <typename TDomain>
bool SplitSynapseHandler<TDomain>::synapse_on_edge(const Edge* edge, size_t scv, number time, number& current)
{
	current = current_on_edge(edge, scv,time);

	return true;
}

template <typename TDomain>
void SplitSynapseHandler<TDomain>::grid_first_available()
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

	m_vAllSynapses = all_synapses();
	m_vPreSynapses = all_pre_synapses();
	m_vPostSynapses = all_post_synapses();

	//map postsyn id to pointer
	for(size_t i=0; i<m_vPostSynapses.size(); ++i) {
		m_mPostSynapses[m_vPostSynapses[i]->id()] = m_vPostSynapses[i];
	}

}

template <typename TDomain>
std::vector<IBaseSynapse*>
SplitSynapseHandler<TDomain>::all_synapses()
{
	std::vector<IBaseSynapse*> vSyn;

	for(Edge* e=m_spGrid->begin<Edge>(0); e != m_spGrid->end<Edge>(0); ++e) {
		for(size_t i = 0; i<m_aaSSyn[e].size(); ++i) {
			vSyn.push_back(m_aaSSyn[e][i]);
		}
	}
	return vSyn;
}

template <typename TDomain>
std::vector<IPreSynapse*>
SplitSynapseHandler<TDomain>::all_pre_synapses()
{
	std::vector<IPreSynapse*> vSyn;

	for(Edge* e=m_spGrid->begin<Edge>(0); e != m_spGrid->end<Edge>(0); ++e) {
		for(size_t i = 0; i<m_aaSSyn[e].size(); ++i) {
			if(m_aaSSyn[e][i]->is_presynapse()) {
				vSyn.push_back( static_cast<IPreSynapse*>(m_aaSSyn[e][i]) );
			}
		}
	}
	return vSyn;
}

template <typename TDomain>
std::vector<IPostSynapse*>
SplitSynapseHandler<TDomain>::all_post_synapses()
{
	std::vector<IPostSynapse*> vSyn;

	for(Edge* e=m_spGrid->begin<Edge>(0); e != m_spGrid->end<Edge>(0); ++e) {
		for(size_t i = 0; i<m_aaSSyn[e].size(); ++i) {
			if(!m_aaSSyn[e][i]->is_presynapse()) {
				vSyn.push_back( static_cast<IPostSynapse*>(m_aaSSyn[e][i]) );
			}
		}
	}
	return vSyn;
}

template <typename TDomain>
void SplitSynapseHandler<TDomain>::
set_ce_object(SmartPtr<CableEquation<TDomain> > disc)
{
	UG_COND_THROW(m_bInited, "The CableEquation object associated to this synapse handler "
				  "must not be changed\nafter addition of the original CableEquation object "
				  "to the domain discretization.");

	// set cable equation disc object
	m_spCEDisc = disc;
}

template <typename TDomain>
void SplitSynapseHandler<TDomain>::
update_presyn(number time)
{
#ifdef UG_PARALLEL
	// collect activated post-synapse ids
	std::vector<unsigned long> vPSID;
	for (size_t i = 0; i < m_vPreSynapses.size(); ++i)
	{
		if (m_vPreSynapses[i]->is_active(time)
			&& m_mActivePreSynapses.find(m_vPreSynapses[i]->id()) == m_mActivePresynapse.end())
		{
			vPSID.push_back(m_vPreSynapses[i]->postsynapse_id());
		}
		else
		{
			// TODO: fill deactivation list
		}
	}

	// communicate

	// TODO: More easily and efficiently using allgatherv?
	//       Then probably without communication of local sizes beforehand.
	unsigned long* locSizes, locSizesRcv;
	size_t sz = pcl::NumProcs();
	locSizes = new unsigned long[sz];
	locSizesRcv = new unsigned long[sz];
	locSizes[pcl::ProcRank()] = vPSID.size();
	com.allreduce(locSizes, locSizesRcv, sz, PCL_DT_UNSIGNED_LONG, PCL_RO_SUM);

	// sum up sizes up to ProcRank()-1
	size_t start_index = 0;
	for (size_t i = 0; i < pcl::ProcRank(); ++i)
		start_index += locSizesRcv[i];

	size_t nPS = start_index;
	for (size_t i = pcl::ProcRank(); i < pcl::NumProcs(); ++i)
		nPS += locSizesRcv[i];

	delete[] locSizes;
	delete[] locSizesRcv;

	unsigned long* locPSID = new unsigned long[nPS];
	unsigned long* globPSID = new unsigned long[nPS];
	for (size_t i = 0; i < nPS; ++i) locPSID[i] = 0;
	for (size_t i = 0; i < vPSID.size(); ++i)
		locPSID[start_index+i] = vPSID[i];

	com.allreduce(locPSID, globPSID, nPS, PCL_DT_UNSIGNED_LONG, PCL_RO_SUM);

	delete[] locPSID;

	for (size_t i = 0; i < nPS; ++i)
	{
		std::map<size_t, IPostSynapse*>::iterator it = m_mPostSynapses.find(globPSID[i]);
		if (it != m_mPostSynapses.end())
			it->second->activate(time);
	}

	delete[] globPSID;


	// TODO: deactivation accordingly

#else
	for(size_t i=0; i<m_vPreSynapses.size(); ++i) {

		//pre synapse becomes active
		if(m_vPreSynapses[i]->is_active(time)) {
			if(!m_mActivePreSynapses.find(m_vPreSynapses[i]->id())) {


				m_mActivePreSynapses[m_vPreSynapses[i]->id()] = m_vPreSynapses[i];
				m_mPostSynapses[m_vPreSynapses[i]->postsynapse_id()]->activate(time);
			}
		//pre synapse becomes inactive
		} else {
			if(m_mActivePreSynapses.find(m_vPreSynapses[i]->id())) {
				m_mActivePreSynapses.erase(m_vPreSynapses[i]->id());
				m_mPostSynapses[m_vPreSynapses[i]->postsynapse_id()]->deactivate(time);
			}
		}
	}
#endif
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
