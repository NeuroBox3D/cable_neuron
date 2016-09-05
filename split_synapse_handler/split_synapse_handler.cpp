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
 m_mPreSynapses(std::map<SYNAPSE_ID,IPreSynapse*>()),
 m_mPostSynapses(std::map<SYNAPSE_ID,IPostSynapse*>()),
 m_mActivePreSynapses(std::map<SYNAPSE_ID, IPreSynapse*>()),
 m_mPreSynapseIdToEdge(std::map<SYNAPSE_ID,Edge*>()),
 m_iterator_type("")
{
	m_ssah.set_attachment(m_aSSyn);
	//std::cout << "\nSSH instantiated\n";
}

template <typename TDomain>
number SplitSynapseHandler<TDomain>::current_on_edge(const Edge* e, size_t scv, number t)
{
	//TODO: update!! currents are nan if synapse inactive
	std::vector<IBaseSynapse*> vSyns(m_aaSSyn[e]);
	number vm_postsyn = m_spCEDisc->vm(e->vertex(scv));
	number curr = 0.0;

	for(size_t i=0; i<vSyns.size(); ++i) {
		IBaseSynapse* s = vSyns[i];
		//postsynapses
		if(s->is_postsynapse()) {
			IPostSynapse* ps = static_cast<IPostSynapse*>(s);

			if(ps->is_active(t)) {
				//std::cout << (ps->is_active(t)) << std::endl;
				//std::cout << s->name() << std::endl;
				curr += ps->current(t, vm_postsyn);
			}
		}
	}
	//std::cout << curr << std::endl << std::endl;
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
	//std::cout << "\nSSH grid_first_available\n";
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

	//gather all synapses from grid and build all tables
	all_synapses();

}

/**
 * returns vector containing all synapses on the grid
 */
template <typename TDomain>
void
SplitSynapseHandler<TDomain>::all_synapses()
{
	for(geometry_traits<Edge>::const_iterator e=m_spGrid->begin<Edge>(0); e != m_spGrid->end<Edge>(0); ++e) {
		for(size_t i = 0; i<m_aaSSyn[*e].size(); ++i) {
			m_vAllSynapses.push_back(m_aaSSyn[*e][i]);
			if(m_aaSSyn[*e][i]->is_presynapse()) { //presynapse map and save edge
				IPreSynapse* s = (IPreSynapse*)m_aaSSyn[*e][i];
				m_mPreSynapses[s->id()] = s;
				m_mPreSynapseIdToEdge[s->id()] = *e;
			} else {  //postsynapse map
				IPostSynapse* s = (IPostSynapse*)m_aaSSyn[*e][i];
				m_mPostSynapses[s->id()] = s;
			}
		}
	}
}


/*template <typename TDomain>
void SplitSynapseHandler<TDomain>::
set_ce_object(SmartPtr<CableEquation<TDomain> > disc)
{
	UG_COND_THROW(m_bInited, "The CableEquation object associated to this synapse handler "
				  "must not be changed\nafter addition of the original CableEquation object "
				  "to the domain discretization.");

	// set cable equation disc object
	m_spCEDisc = disc;
}*/



template <typename TDomain>
void SplitSynapseHandler<TDomain>::
update_presyn(number time)
{
#ifdef UG_PARALLEL

	//collect postsynapse id's that became active/inactive
	std::vector<SYNAPSE_ID> vNewActivePostSynapseIds_local; //becoming active on local process
	std::vector<SYNAPSE_ID> vNewInactivePostSynapseIds_local; //becomming inactive on local process

	std::vector<SYNAPSE_ID> vNewActivePostSynapseIds_global; //complete list of all newly active synapses
	std::vector<SYNAPSE_ID> vNewInactivePostSynapseIds_global; //complete list of all newly inactive synapses

	for(std::map<SYNAPSE_ID, IPreSynapse*>::iterator it = m_mPreSynapses.begin(); it != m_mPreSynapses.end(); ++it) {
		IPreSynapse* s = it->second;
		Edge* e = m_mPreSynapseIdToEdge[it->first] ; //To s corresponding Edge e

		Vertex* v1 = e->vertex(0);
		Vertex* v2 = e->vertex(1);

		// TODO: get other values too (like in commented-out code snippet below)
		std::vector<number> uAtSynapseLocation(1);
		uAtSynapseLocation[0] = s->location() * m_spCEDisc->vm(v2)
								+ (1.0 - s->location()) *  m_spCEDisc->vm(v1);
		/*
		for (size_t fct = 0; fct <= m_spCEDisc->m_numb_ion_funcs; ++fct)
		{
			std::vector<DoFIndex> dofIndex1;
			std::vector<DoFIndex> dofIndex2;
			dd->dof_indices(v1, fct, dofIndex1, false, false);
			dd->dof_indices(v2, fct, dofIndex2, false, false);

			UG_COND_THROW(dofIndex1.size() != 1, "Not exactly one dof index on vertex 0 for function " << fct << ".");
			UG_COND_THROW(dofIndex2.size() != 1, "Not exactly one dof index on vertex 1 for function " << fct << ".");

			uAtSynapseLocation[fct] = s->location() * DoFRef(u, dofIndex2[0])
									  + (1.0 - s->location()) * DoFRef(u, dofIndex1[0]);
		}
		*/

		s->update(time, uAtSynapseLocation);

		//pre synapse becomes just active if it is active and couldn't be found in the active presynapses map
		if( (s->is_active(time)) && (m_mActivePreSynapses.find(s->id()) == m_mActivePreSynapses.end()) ) {
			m_mActivePreSynapses[s->id()] = s;
			vNewActivePostSynapseIds_local.push_back(s->id());

		//pre synapse becomes inactive
		} else if (!s->is_active(time) && m_mActivePreSynapses.find(s->id()) != m_mActivePreSynapses.end()) {

			m_mActivePreSynapses.erase(s->id());
			vNewInactivePostSynapseIds_local.push_back(s->id());
		}
	}

	//propagate to all processes
	pcl::ProcessCommunicator com;
	com.allgatherv(vNewActivePostSynapseIds_global, vNewActivePostSynapseIds_local);
	com.allgatherv(vNewInactivePostSynapseIds_global, vNewInactivePostSynapseIds_local);

	//todo: maybe assume a structured distribution to increase performance of the scanning, tbd
	//scan for synapse ids that are on the local process and have to be activated
	for(size_t i=0; i<vNewActivePostSynapseIds_global.size(); ++i) {
		SYNAPSE_ID psid = vNewActivePostSynapseIds_global[i];

		if(m_mPostSynapses.find(psid) != m_mPostSynapses.end()) { //synapse id found on local process
			m_mPostSynapses[psid]->activate(time);
		}
	}

	//scan for synapse ids that are on the local process and have to be deactivated
	for(size_t i=0; i<vNewInactivePostSynapseIds_global.size(); ++i) {
		SYNAPSE_ID psid = vNewInactivePostSynapseIds_global[i];

		if(m_mPostSynapses.find(psid) != m_mPostSynapses.end()) {
			m_mPostSynapses[psid]->deactivate(time);
		}
	}




	/*// collect activated post-synapse ids
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


	// TODO: deactivation accordingly*/


#else

	for(std::map<SYNAPSE_ID, IPreSynapse*>::iterator it = m_mPreSynapses.begin(); it != m_mPreSynapses.end(); ++it) {
		IPreSynapse* s = it->second;
		Edge* e = m_mPreSynapseIdToEdge[it->first] ; //To s corresponding Edge e

		Vertex* v1 = e->vertex(0);
		Vertex* v2 = e->vertex(1);

		// TODO: get other values too (like in commented-out code snippet below)
		std::vector<number> uAtSynapseLocation(1);
		uAtSynapseLocation[0] = s->location() * m_spCEDisc->vm(v2)
								+ (1.0 - s->location()) *  m_spCEDisc->vm(v1);
		/*
		for (size_t fct = 0; fct <= m_spCEDisc->m_numb_ion_funcs; ++fct)
		{
			std::vector<DoFIndex> dofIndex1;
			std::vector<DoFIndex> dofIndex2;
			dd->dof_indices(v1, fct, dofIndex1, false, false);
			dd->dof_indices(v2, fct, dofIndex2, false, false);

			UG_COND_THROW(dofIndex1.size() != 1, "Not exactly one dof index on vertex 0 for function " << fct << ".");
			UG_COND_THROW(dofIndex2.size() != 1, "Not exactly one dof index on vertex 1 for function " << fct << ".");

			uAtSynapseLocation[fct] = s->location() * DoFRef(u, dofIndex2[0])
									  + (1.0 - s->location()) * DoFRef(u, dofIndex1[0]);
		}
		*/

		s->update(time, uAtSynapseLocation);

		//pre synapse becomes just active if it is active and couldn't be found in the active presynapses map
		if( (s->is_active(time)) && (m_mActivePreSynapses.find(s->id()) == m_mActivePreSynapses.end()) ) {
			m_mActivePreSynapses[s->id()] = s;
			m_mPostSynapses[s->id()]->activate(time);

		//pre synapse becomes inactive
		} else if (!s->is_active(time) && m_mActivePreSynapses.find(s->id()) != m_mActivePreSynapses.end()) {
			m_mActivePreSynapses.erase(s->id());
			m_mPostSynapses[s->id()]->deactivate(time);
		}
	}

#endif
}


// ////////////////////////////////////
//	explicit template instantiations //
// ////////////////////////////////////

#ifdef UG_DIM_1
	template class SplitSynapseHandler<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class SplitSynapseHandler<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class SplitSynapseHandler<Domain3d>;
#endif

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
