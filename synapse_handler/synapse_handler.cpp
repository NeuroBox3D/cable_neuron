/*!
 * \file synapse_handler.cpp
 */

#include "synapse_handler.h"


namespace ug {
namespace cable_neuron {
namespace synapse_handler {


template <typename TDomain>
bool SynapseDistributorSynapseHandler<TDomain>::
synapse_on_edge(const Edge* edge, size_t scv, number time, number& current)
{
	const SynapseDistributor* pSD = m_spSD.get();
	if (!pSD->has_active_synapses(edge, scv, time, current)) {
		current = 0;
		return false;
	}

	return true;
}

template <typename TDomain>
bool SynapseDistributorSynapseHandler<TDomain>::
synapse_at_location(const ug::MathVector<TDomain::dim>& vec, number time, number& current) const
{
	const SynapseDistributor* pSD = m_spSD.get();
	if (!pSD->has_active_synapses(vec, time, current)) {
		current = 0;
		return false;
	}

	return true;
}

template <typename TDomain>
void SynapseDistributorSynapseHandler<TDomain>::set_sd(ConstSmartPtr<SynapseDistributor> sd)
{
	this->m_spSD = sd;
}

template <typename TDomain>
std::string SynapseDistributorSynapseHandler<TDomain>::
name() const
{
	return "SynapseDistributor Synapse Handler";
}


// ////////////////////////////////////
//	explicit template instantiations //
// ////////////////////////////////////

#ifdef UG_DIM_1
	template class SynapseDistributorSynapseHandler<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class SynapseDistributorSynapseHandler<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class SynapseDistributorSynapseHandler<Domain3d>;
#endif






template <typename TDomain>
NETISynapseHandler<TDomain>::
NETISynapseHandler()
: m_spCEDisc(SPNULL), m_spApprox(SPNULL), m_spGrid(SPNULL),
  m_aSynInfo(GlobalAttachments::attachment<AVSynapse>("Synapses")),
  m_aPSI(GlobalAttachments::attachment<AUInt>("presyn_index")),
  m_presynSubset(""), m_presynSI(-1),
  m_start_time(std::numeric_limits<number>::max()), m_start_time_dev(0.0),
  m_duration(0.0), m_duration_dev(0.0),
  m_peak_cond(6e-4),
  m_constSeed(true),
  m_bInited(false),
  m_bEXP2_SYNAPSE(true)
{
	m_cah.set_attachment(m_aPSI);
	m_siah.set_attachment(m_aSynInfo);
}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
set_presyn_subset(const char* presynSubset)
{
	std::vector<std::string> vSss = TokenizeString(presynSubset, ',');
	UG_COND_THROW(vSss.size() != 1, "Exactly one subset containing the presynapses "
				  "must be indicated.");

	m_presynSubset = vSss[0];
}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
set_ce_object(SmartPtr<CableEquation<TDomain> > disc)
{
	UG_COND_THROW(m_bInited, "The CableEquation object associated to this synapse handler "
				  "must not be changed\nafter addition of the original CableEquation object "
				  "to the domain discretization.");

	// set cable equation disc object
	m_spCEDisc = disc;
}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
set_activation_timing
(
	number start_time,
	number duration,
	number start_time_dev,
	number duration_dev,
	number peak_cond,
	bool constSeed
)
{
	UG_COND_THROW(m_bInited, "The activation timing cannot be changed after addition of the\n"
				  "original CableEquation object to the domain discretization.");

	m_start_time = start_time;
	m_duration = duration;
	m_start_time_dev = start_time_dev;
	m_duration_dev = duration_dev;
	m_peak_cond = peak_cond;
	m_constSeed = constSeed;
}



template <typename TDomain>
void NETISynapseHandler<TDomain>::
grid_first_available()
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
	if (!GlobalAttachments::is_declared("Synapses")) {
		UG_THROW("GlobalAttachment 'Synapses' not available.");
	}

	if (!GlobalAttachments::is_declared("presyn_index")) {
		UG_THROW("GlobalAttachment 'presyn_index' not available.");
	}

	// Attach to grid, if not already present
	if (!m_spGrid->has_edge_attachment(m_aSynInfo)) {
		m_spGrid->attach_to_edges(m_aSynInfo);
	}

	if (!m_spGrid->has_vertex_attachment(m_aPSI)) {
		m_spGrid->attach_to_vertices(m_aPSI);
	}

	// check that essential attachments exist in grid, create accessors
	UG_COND_THROW(!m_spGrid->has_attachment<Edge>(m_aSynInfo), "No synapse info attached to grid!");
	UG_COND_THROW(!m_spGrid->has_attachment<Vertex>(m_aPSI), "No presynaptic index attached to grid!");

	m_aaSynapseInfo = Grid::EdgeAttachmentAccessor<AVSynapse>(*m_spGrid, m_aSynInfo);
	m_aaPSI = Grid::VertexAttachmentAccessor<AUInt>(*m_spGrid, m_aPSI);

	// set activation timing for AlphaSynapses
	set_activation_timing_with_grid();

	// set all Exp2Syn synapses to deactivated status (currently nor guaranteed by grid)
	deactivate_all_biexp();

	// propagate attachments through levels
	m_cah.set_grid(m_spGrid);
	m_siah.set_grid(m_spGrid);

	// presynapic vertex indexing in case of EXP2_SYNAPSE presence
	if(!has_EXP2_SYNAPSE())
	{
		m_bEXP2_SYNAPSE = false;
	}
	else
	{
		// set presynaptic subset index (from name given in constructor)
		UG_COND_THROW(m_presynSubset.empty(), "No presynaptic subset has been set.\n"
					  "Use set_presyn_subset() method to do so.");

		typedef typename TDomain::subset_handler_type sh_type;
		ConstSmartPtr<sh_type> sh = m_spApprox->domain()->subset_handler();
		UG_COND_THROW(!sh.valid(), "The subset handler assigned to the domain taken from the passed\n"
					  "CableEquation object is not valid.")

		bool found = false;
		for (int si = 0; si < sh->num_subsets(); ++si)
		{
			if (strcmp (m_presynSubset.c_str(), sh->get_subset_name(si)) == 0)
			{
				m_presynSI = si;
				found = true;
				break;
			}
		}
		UG_COND_THROW(!found, "The subset "<< m_presynSubset << "' is not defined in the domain taken from\n"
					  "the passed CableEquation object.");

		// find number of presynaptic indices
		resize_presyn_vector();
	}

	// set init'ed flag
	m_bInited = true;
}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
update_presyn()
{
//	update only in case of EXP2_SYNAPSE presence
	if(m_bEXP2_SYNAPSE == true)
	{
		// check that handler has been init'ed
		UG_COND_THROW(!m_bInited, "Cannot update before initialization.");

		// TODO: ensure consistency!?

		// reset presynaptic index values
		m_vPresynVmValues.assign(m_vPresynVmValues.size(), -std::numeric_limits<number>::max());

		// loop presynaptic subset and fill known values
		MGSubsetHandler& ssh = *m_spCEDisc->approx_space()->domain()->subset_handler();

		const std::vector<Vertex*>& sfv = m_spCEDisc->surface_vertices();
		size_t sfv_sz = sfv.size();
		for (size_t sv = 0; sv < sfv_sz; ++sv)
		{
			Vertex* vrt = sfv[sv];
			if (ssh.get_subset_index(vrt) != m_presynSI) continue;

			uint idx = m_aaPSI[vrt];
			m_vPresynVmValues[idx] = m_spCEDisc->vm(vrt);	// todo: adapt SH units to CE !
		}

		// communicate (all-to-all)
#ifdef UG_PARALLEL
		//UG_LOG_ALL_PROCS("m_vPresynVmValues.size() = " << m_vPresynVmValues.size() << std::endl);
		if (pcl::NumProcs() > 1 && m_vPresynVmValues.size())
		{
			pcl::ProcessCommunicator com;

			number* localData;
			size_t sz = m_vPresynVmValues.size();
			localData = new number[sz];
			memcpy(localData, &m_vPresynVmValues[0], sizeof(number)*sz);
			com.allreduce(localData, &m_vPresynVmValues[0],
						  sz, PCL_DT_DOUBLE, PCL_RO_MAX);
			delete[] localData;
		}

#ifndef NDEBUG
		// for debugging purposes: check that all values have been set
		for (size_t i = 0; i < m_vPresynVmValues.size(); ++i)
		{
			UG_ASSERT(m_vPresynVmValues[i] != -std::numeric_limits<number>::max(),
					  "Value for presynaptic index " << i << " is not set.");
		}
#endif
#endif
	}
}


template <typename TDomain>
bool NETISynapseHandler<TDomain>::
synapse_on_edge(const Edge* edge, size_t scv, number time, number& current)
{
	// check that handler has been init'ed
	UG_COND_THROW(!m_bInited, "Cannot provide synapse information before being init'ed.");

	// get edge synapse info
	std::vector<SynapseInfo>& vInfo = m_aaSynapseInfo[edge];

	bool active_synapse = false;
	current = 0.0;		// todo: adapt SH units to CE! here: current in units of A instead of nA

	typedef synapse_traits<> STV;
	typedef synapse_traits<AlphaSynapse> STA;
	typedef synapse_traits<Exp2Syn> STB;

	for (size_t i = 0; i < vInfo.size(); ++i)
	{
		SynapseInfo& info = vInfo[i];

		if ((STV::loc_coord(info) < 0.5 && scv == 0) || (STV::loc_coord(info) >= 0.5 && scv == 1))
		{
			// get vmDisc potential values for edge
			number vm_postsyn = m_spCEDisc->vm(edge->vertex(scv));

			switch (STV::type(info))
			{
				// treat alpha synapse
				case ALPHA_SYNAPSE:
				{
					//if (SynapseSelector::is_declared(ALPHA_SYNAPSE))
					//{
						// create a default alpha synapse and populate with parameters
						AlphaSynapse alpha(STA::g_max(info), STA::onset(info), STA::tau(info), STA::v_rev(info), vm_postsyn);
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
					//if (SynapseSelector::is_declared(EXP2_SYNAPSE))
					//{
						// create a default exp2syn and populate with parameters
						Exp2Syn exp(STB::tau1(info), STB::tau2(info), STB::v_rev(info), STB::g_max(info), vm_postsyn);

						// get potential for presynapse
						unsigned int presynInd = STB::presyn_ind(info);
						UG_ASSERT(presynInd < m_vPresynVmValues.size(),
									"Requested presynaptic potential value for unknown index "
									<< presynInd << " (only " << m_vPresynVmValues.size()
									<< " presynaptic potential values are available)!");
						number vm_presyn = m_vPresynVmValues[presynInd];

						// activate synapse if V_m > threshold
						if (vm_presyn > STB::threshold(info) && !STB::activated(info))
							STB::activate(info, time);

						// handle synapse time
						if (STB::activated(info))
						{
							// active synapse
							//UG_LOG_ALL_PROCS("Active synapse " << STB::name() << "!" << std::endl);
							number syn_time = time - STB::onset(info);


							// deactivate if time past activity
							if (syn_time >= STB::activity_time(info))
								STB::deactivate(info);

							current += exp.current(syn_time);
							active_synapse = true;
						}
						//UG_LOG_ALL_PROCS("current current (exp2syn): " << current << std::endl);
					//}
					break;
				}

				// treat default synapse
				default:
					break;
			}
			return true;
		}
	}

	return active_synapse;
}


template <typename TDomain>
bool NETISynapseHandler<TDomain>::
synapse_at_location(const MathVector<TDomain::dim>& vec, number time, number& current) const
{
	// check that handler has been init'ed
	UG_COND_THROW(!m_bInited, "Cannot provide synapse information before being init'ed.");

	/// TODO implement
	current = 0;
	UG_THROW("Not implemented!");
	return false;
}

template <typename TDomain>
std::string NETISynapseHandler<TDomain>::
name() const
{
	return "NETI Synapse Handler";
}



template <typename TDomain>
void NETISynapseHandler<TDomain>::
deactivate_all_biexp()
{
	typedef geometry_traits<Edge>::const_iterator iter_type;
	iter_type eIter = m_spGrid->begin<Edge>(0);
	iter_type eEnd = m_spGrid->end<Edge>(0);

	typedef synapse_traits<Exp2Syn> STB;
	typedef synapse_traits<void> STV;

	// loop base level edges and deactivate all biexp synapses
	for (; eIter != eEnd ; ++eIter)
	{
		Edge* e = *eIter;
		std::vector<SynapseInfo>& vSI = m_aaSynapseInfo[e];
		for (size_t i = 0; i < vSI.size(); ++i)
		{
			SynapseInfo& info = vSI[i];

			// only treat alpha synapses
			if (STV::type(info) != EXP2_SYNAPSE)
				continue;

			STB::deactivate(info);
		}
	}
}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
set_activation_timing_with_grid()
{
	// activity timing setup: random normal distribution
	boost::mt19937 rng;
#ifdef UG_PARALLEL
	if (m_constSeed)
		rng.seed(0);//pcl::ProcRank());
	else
		rng.seed(pcl::ProcRank()*time(NULL));
#else
	if (m_constSeed)
		rng.seed(0);
	else
		rng.seed(time(NULL));
#endif

	boost::normal_distribution<double> start_dist(m_start_time, m_start_time_dev);
	boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > var_start(rng, start_dist);

	boost::normal_distribution<double> duration_dist(m_duration, m_duration_dev);
	boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > var_duration(rng, duration_dist);

	// check availability of all needed structures
	UG_COND_THROW(!m_spGrid.valid(), "No valid grid. Make sure that the synapse handler has a grid to work on!");

	typedef geometry_traits<Edge>::const_iterator iter_type;
	iter_type eIter = m_spGrid->begin<Edge>(0);
	iter_type eEnd = m_spGrid->end<Edge>(0);

	typedef synapse_traits<AlphaSynapse> STA;
	typedef synapse_traits<void> STV;

	// loop edges and change alpha synapses
	for (; eIter != eEnd ; ++eIter)
	{
		Edge* e = *eIter;
		std::vector<SynapseInfo>& vSI = m_aaSynapseInfo[e];
		for (size_t i = 0; i < vSI.size(); ++i)
		{
			SynapseInfo& info = vSI[i];

			// only treat alpha synapses
			if (STV::type(info) != ALPHA_SYNAPSE)
				continue;

			number t_start = var_start();

			if (t_start < 0)
				t_start = 0;

			number dur = var_duration();
			dur = dur < 0 ? 0 : dur;

			// set start and duration time
			STA::onset(info) = t_start;
			// This will parameterize the alpha_synapse in such a way that
			// at time = 6*tau, the current will be < 0.05 times its maximal strength.
			STA::tau(info) = dur / 6.0;
			STA::nSpikes(info) = 1;
			STA::freq(info) = 1;
			STA::g_max(info) = m_peak_cond;
		}
	}
}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
print_biexp_info(size_t syn_index)
{
	// check availability of all needed structures
	UG_COND_THROW(!m_spGrid.valid(), "No valid grid. Make sure that the synapse handler has a grid to work on!");

	typedef geometry_traits<Edge>::const_iterator iter_type;
	iter_type eIter = m_spGrid->begin<Edge>(0);
	iter_type eEnd = m_spGrid->end<Edge>(0);

	typedef synapse_traits<Exp2Syn> STB;
	typedef synapse_traits<void> STV;

	// loop edges and change alpha synapses
	for (; eIter != eEnd ; ++eIter)
	{
		Edge* e = *eIter;
		std::vector<SynapseInfo>& vSI = m_aaSynapseInfo[e];
		for (size_t i = 0; i < vSI.size(); ++i)
		{
			SynapseInfo& info = vSI[i];

			// only treat exp2 synapses
			if (STV::type(info) != EXP2_SYNAPSE)
				continue;
			// find the one with correct index
			if (STB::presyn_ind(info) != syn_index)
				continue;

			// print info
			std::cout << std::endl << "Info for exp2-synapse " << syn_index << ":" << std::endl;
			info.print_to_log();
			return;
		}
	}
}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
print_synapse_statistics(size_t soma_si)
{
	// check availability of all needed structures
	UG_COND_THROW(!m_bInited, "print_synapse_statistics is not yet init'ed.");

	typename TDomain::position_accessor_type& aaPos
		= m_spApprox->domain()->position_accessor();
	UG_COND_THROW(TDomain::position_type::Size < 3, "Dimension must be 3d.");

	const size_t NT_L23 = 0;
	const size_t NT_L4 = 1;
	const size_t NT_L5A = 2;
	const size_t NT_L5B = 3;

	typedef synapse_traits<AlphaSynapse> STA;
	typedef synapse_traits<Exp2Syn> STB;
	typedef synapse_traits<void> STV;

	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	edge_list el;

	// vector holding neuron indices (to be accorded) for presynapse indices
	std::vector<size_t> vNeuronInd(m_vPresynVmValues.size(), -1);

	// vector holding number of postsynapses connected to a presynaptic index
	std::vector<unsigned long> vNumPostSyn(m_vPresynVmValues.size(), 0);

	// vector holding type of neuron
	std::vector<size_t> vNeuronType;

	// information regarding exp2 synapses
	// (first: neuron index postsynapse; second: synapse number) holds: presyn index
	std::vector<std::vector<size_t> > vSynInfo;

	// information regarding alpha synapses
	std::vector<size_t> vASynCount;

	// marker for treated vertices
	m_spGrid->begin_marking();

	// iterate over somatic subset
	typedef geometry_traits<Edge>::const_iterator iter_type;
	MGSubsetHandler& ssh = *m_spCEDisc->approx_space()->domain()->subset_handler();
	iter_type eIter = ssh.begin<Edge>(soma_si, 0);
	iter_type eEnd = ssh.end<Edge>(soma_si, 0);

	size_t nNeuron = 0;
	for (; eIter != eEnd ; ++eIter)
	{
		Edge* e = *eIter;

		// only treat untreated edges
		if (m_spGrid->is_marked(e)) continue;

		// we are on a new neuron
		size_t nInd = nNeuron++;

		// save vertex 0 as soma vertex for this neuron
		m_vSomaVertices.push_back(e->vertex(0));

		// get type by coordinate (empirically determined from NeuGen)
		size_t type;
		Vertex* v = e->vertex(0);
		number zcoord = aaPos[v][2];
		if (zcoord > 0.00095) type = NT_L23;
		else if (zcoord > 0.00071) type = NT_L4;
		else if (zcoord > 0.00037) type = NT_L5A;
		else type = NT_L5B;

		vNeuronType.push_back(type);
		vSynInfo.push_back(std::vector<size_t>());
		vASynCount.push_back(0);
		std::vector<size_t>& eSynVec = vSynInfo.back();
		size_t& aSynCount = vASynCount.back();

		// treat all connected edges and vertices of neuron
		std::queue<std::pair<Edge*, Vertex*> > q;
		q.push(std::make_pair(e, (Vertex*)NULL));

		// work off queue
		while (!q.empty())
		{
			std::pair<Edge*, Vertex*> ev = q.front();
			q.pop();

			e = ev.first;
			v = ev.second;

			m_spGrid->mark(e);

			// save info for possible postsynapses
			std::vector<SynapseInfo>& vSI = m_aaSynapseInfo[e];
			for (size_t i = 0; i < vSI.size(); ++i)
			{
				SynapseInfo& info = vSI[i];

				// only treat exp2 synapses
				switch (STV::type(info))
				{
					case ALPHA_SYNAPSE:
						++aSynCount;
						m_locASynInfo.push_back(LocalASynapseInfo(e,i,type));
						break;
					case EXP2_SYNAPSE:
						++(vNumPostSyn[STB::presyn_ind(info)]);
						eSynVec.push_back(STB::presyn_ind(info));
						m_locESynInfo.push_back(LocalESynapseInfo(e,i,-1,type));
						break;
					default:
						UG_LOG("Unknown synapse type");
				}
			}

			// treat vertices of edge
			Vertex* v0 = e->vertex(0);
			if (v0 != v)
			{
				// save neuron index for possible presynapse
				if (ssh.get_subset_index(v0) == m_presynSI)
					vNeuronInd[m_aaPSI[v0]] = nInd;

				// add other associated edges to queue
				m_spApprox->domain()->grid()->associated_elements(el, v0);
				for (size_t k = 0; k < el.size(); ++k)
					if (el[k] != e && !m_spGrid->is_marked(el[k])) q.push(std::make_pair(el[k], v0));
			}
			v0 = e->vertex(1);
			if (v0 != v)
			{
				// save neuron index for possible presynapse
				if (ssh.get_subset_index(v0) == m_presynSI)
					vNeuronInd[m_aaPSI[v0]] = nInd;

				// add other associated edges to queue
				m_spApprox->domain()->grid()->associated_elements(el, v0);
				for (size_t k = 0; k < el.size(); ++k)
					if (el[k] != e && !m_spGrid->is_marked(el[k])) q.push(std::make_pair(el[k], v0));
			}
		}
	}

	m_spGrid->end_marking();

	// copy neuron type vector to local variable
	m_vNeuronType = vNeuronType;

	// vector holding neuron types for presynapse indices
	std::vector<unsigned long> vPresynNType(m_vPresynVmValues.size(), -1);
	size_t sz_ps = vPresynNType.size();
	for (size_t i = 0; i < sz_ps; ++i)
		if (vNeuronInd[i] != (size_t) -1)
			vPresynNType[i] = vNeuronType[vNeuronInd[i]];

	// communicate
#ifdef UG_PARALLEL
	size_t nProcs = pcl::NumProcs();
	if (nProcs > 1 && vPresynNType.size())
	{
		pcl::ProcessCommunicator com;

		// vPresynNType
		unsigned long* localData;
		size_t sz = vPresynNType.size();
		localData = new unsigned long[sz];
		memcpy(localData, &vPresynNType[0], sizeof(unsigned long)*sz);
		com.allreduce(localData, &vPresynNType[0], sz, PCL_DT_UNSIGNED_LONG, PCL_RO_MIN);

		// vNeuronInd (indices are only valid locally, atm;
		// we can make them unique as we do not need their local values any longer)
		unsigned long locNumNeuron = vNeuronType.size();
		std::vector<unsigned long> vNeuronCounts(nProcs,0);
		com.allgather(&locNumNeuron, 1, PCL_DT_UNSIGNED_LONG,
					  &vNeuronCounts[0], 1, PCL_DT_UNSIGNED_LONG);

		size_t offset = 0;
		size_t rk = pcl::ProcRank();
		for (size_t i = 0; i < rk; ++i)
			offset += vNeuronCounts[i];

		for (size_t i = 0; i < sz; ++i)
			if (vNeuronInd[i] != (size_t) -1)
				vNeuronInd[i] += offset;

		memcpy(localData, &vNeuronInd[0], sizeof(unsigned long)*sz);
		com.allreduce(localData, &vNeuronInd[0], sz, PCL_DT_UNSIGNED_LONG, PCL_RO_MIN);
		delete[] localData;
	}
#endif

	// for debugging purposes: check that all values have been set
	for (size_t i = 0; i < vPresynNType.size(); ++i)
	{
		UG_COND_THROW(vNeuronInd[i] == (size_t) -1, "Neuron index for presynaptic index " << i << " is not set.");
		UG_COND_THROW(vPresynNType[i] == (size_t) -1, "Neuron type for presynaptic index " << i << " is not set.");
	}

	// now evaluate stored info
	size_t count[4];

	size_t minACount[4];
	size_t maxACount[4];
	size_t sumACount[4];
	number avgACount[4];

	size_t minECount[4][4];
	size_t maxECount[4][4];
	size_t sumECount[4][4];
	number avgECount[4][4];

	size_t minECountNeuronWise[4][4];
	size_t maxECountNeuronWise[4][4];
	size_t sumECountNeuronWise[4][4];
	number avgECountNeuronWise[4][4];

	for (size_t i = 0; i < 4; ++i)
	{
		count[i] = 0;
		minACount[i] = std::numeric_limits<size_t>::max();
		maxACount[i] = sumACount[i] = 0;
		for (size_t j = 0; j < 4; ++j)
		{
			minECount[i][j] = minECountNeuronWise[i][j] = std::numeric_limits<size_t>::max();
			maxECount[i][j] = sumECount[i][j] = 0;
			maxECountNeuronWise[i][j] = sumECountNeuronWise[i][j] = 0;
		}
	}

	size_t locESynInd = 0;
	size_t sz = vNeuronType.size();
	for (size_t n = 0; n < sz; ++n)
	{
		size_t type = vNeuronType[n];
		++(count[type]);

		size_t aCount = vASynCount[n];
		minACount[type] = std::min(aCount, minACount[type]);
		maxACount[type] = std::max(aCount, maxACount[type]);
		sumACount[type] += aCount;

		std::vector<size_t>& vESyn = vSynInfo[n];
		// sort vector of exp2syn for this neuron on presyn_neuron_index
		MyIndexSort myIndexSort(vNeuronInd);
		std::sort(vESyn.begin(), vESyn.end(), myIndexSort);
		size_t lastPSNeuronInd = vESyn.size() ? vESyn[0]+1 : 0;

		size_t sz2 = vESyn.size();
		size_t eCount[4];
		size_t eCountNeuronWise[4];
		for (size_t i = 0; i < 4; ++i) eCount[i] = eCountNeuronWise[i] = 0;
		for (size_t s = 0; s < sz2; ++s)
		{
			size_t type2 = vPresynNType[vESyn[s]];
			++(eCount[type2]);
			if (vNeuronInd[vESyn[s]] != lastPSNeuronInd)
			{
				++(eCountNeuronWise[type2]);
				lastPSNeuronInd = vNeuronInd[vESyn[s]];
			}
			m_locESynInfo[locESynInd].from = type2;
			++locESynInd;
		}

		for (size_t i = 0; i < 4; ++i)
		{
			minECount[type][i] = std::min(eCount[i], minECount[type][i]);
			maxECount[type][i] = std::max(eCount[i], maxECount[type][i]);
			sumECount[type][i] += eCount[i];

			minECountNeuronWise[type][i] = std::min(eCountNeuronWise[i], minECountNeuronWise[type][i]);
			maxECountNeuronWise[type][i] = std::max(eCountNeuronWise[i], maxECountNeuronWise[type][i]);
			sumECountNeuronWise[type][i] += eCountNeuronWise[i];
		}
	}

	// communicate results
#ifdef UG_PARALLEL
	if (nProcs > 1 && vPresynNType.size())
	{
		pcl::ProcessCommunicator com;

		unsigned long* locSumData = new unsigned long[40];
		unsigned long* locMinData = new unsigned long[36];
		unsigned long* locMaxData = new unsigned long[36];
		unsigned long* globSumData = new unsigned long[40];
		unsigned long* globMinData = new unsigned long[36];
		unsigned long* globMaxData = new unsigned long[36];

		// summable data
		memcpy(locSumData, &count[0], sizeof(unsigned long)*4);
		memcpy(locSumData+4, &sumACount[0], sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(locSumData+8+4*i, &sumECount[i][0], sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(locSumData+24+4*i, &sumECountNeuronWise[i][0], sizeof(unsigned long)*4);

		com.allreduce(locSumData, globSumData, 40, PCL_DT_UNSIGNED_LONG, PCL_RO_SUM);

		memcpy(&count[0], globSumData, sizeof(unsigned long)*4);
		memcpy(&sumACount[0], globSumData+4, sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(&sumECount[i][0], globSumData+8+4*i, sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(&sumECountNeuronWise[i][0], globSumData+24+4*i, sizeof(unsigned long)*4);

		// min-able data
		memcpy(locMinData, &minACount[0], sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(locMinData+4+4*i, &minECount[i][0], sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(locMinData+20+4*i, &minECountNeuronWise[i][0], sizeof(unsigned long)*4);

		com.allreduce(locMinData, globMinData, 36, PCL_DT_UNSIGNED_LONG, PCL_RO_MIN);

		memcpy(&minACount[0], globMinData, sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(&minECount[i][0], globMinData+4+4*i, sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(&minECountNeuronWise[i][0], globMinData+20+4*i, sizeof(unsigned long)*4);

		// max-able data
		memcpy(locMaxData, &maxACount[0], sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(locMaxData+4+4*i, &maxECount[i][0], sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(locMaxData+20+4*i, &maxECountNeuronWise[i][0], sizeof(unsigned long)*4);

		com.allreduce(locMaxData, globMaxData, 36, PCL_DT_UNSIGNED_LONG, PCL_RO_MAX);

		memcpy(&maxACount[0], globMaxData, sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(&maxECount[i][0], globMaxData+4+4*i, sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(&maxECountNeuronWise[i][0], globMaxData+20+4*i, sizeof(unsigned long)*4);

		delete[] locSumData;
		delete[] locMinData;
		delete[] locMaxData;
		delete[] globSumData;
		delete[] globMinData;
		delete[] globMaxData;

		unsigned long* localPresynCounts;
		size_t sz = vNumPostSyn.size();
		localPresynCounts = new unsigned long[sz];
		memcpy(localPresynCounts, &vNumPostSyn[0], sizeof(unsigned long)*sz);
		com.allreduce(localPresynCounts, &vNumPostSyn[0], sz, PCL_DT_UNSIGNED_LONG, PCL_RO_SUM);
		delete[] localPresynCounts;
	}
#endif

	for (size_t i = 0; i < 4; ++i)
	{
		avgACount[i] = sumACount[i] / (number)count[i];
		for (size_t j = 0; j < 4; ++j)
		{
			avgECount[i][j] = sumECount[i][j] / (number)count[i];
			avgECountNeuronWise[i][j] = sumECountNeuronWise[i][j] / (number)count[i];
		}
	}

	std::sort(vNumPostSyn.begin(), vNumPostSyn.end());

	// print out info
	UG_LOGN("-------------------------");
	UG_LOGN("-- synapse information --");
	UG_LOGN("-------------------------");
	UG_LOGN("");
	UG_LOGN("number of neurons:");
	UG_LOGN("  L2/3  "<<std::setw(5)<<std::setfill(' ')<<count[NT_L23]);
	UG_LOGN("  L4    "<<std::setw(5)<<std::setfill(' ')<<count[NT_L4]);
	UG_LOGN("  L5A   "<<std::setw(5)<<std::setfill(' ')<<count[NT_L5A]);
	UG_LOGN("  L5B   "<<std::setw(5)<<std::setfill(' ')<<count[NT_L5B]);
	UG_LOGN("");

	size_t sz_nps = vNumPostSyn.size();
	UG_LOGN("number of presynapse indices: " << sz_nps);
	UG_LOGN("maximal number of postsynapses connected to a presynapse: "
			<< (sz_nps ? vNumPostSyn[sz_nps-1] : -1));
	UG_LOGN("minimal number of postsynapses connected to a presynapse: "
			<< (sz_nps ? vNumPostSyn[0] : -1));
	UG_LOGN("");
	UG_LOGN("number of alpha synapses:");
	UG_LOGN("          min    max    avg");
	UG_LOGN("  L2/3  "<<std::setw(5)<<std::setfill(' ') <<minACount[NT_L23]<<"  "
					  <<std::setw(5)<<std::setfill(' ') <<maxACount[NT_L23]<<"  "
					  <<std::setw(5)<<std::setfill(' ') <<avgACount[NT_L23]);
	UG_LOGN("  L4    "<<std::setw(5)<<std::setfill(' ') <<minACount[NT_L4]<<"  "
					  <<std::setw(5)<<std::setfill(' ') <<maxACount[NT_L4]<<"  "
					  <<std::setw(5)<<std::setfill(' ') <<avgACount[NT_L4]);
	UG_LOGN("  L5A   "<<std::setw(5)<<std::setfill(' ') <<minACount[NT_L5A]<<"  "
					  <<std::setw(5)<<std::setfill(' ') <<maxACount[NT_L5A]<<"  "
					  <<std::setw(5)<<std::setfill(' ') <<avgACount[NT_L5A]);
	UG_LOGN("  L5B   "<<std::setw(5)<<std::setfill(' ') <<minACount[NT_L5B]<<"  "
					  <<std::setw(5)<<std::setfill(' ') <<maxACount[NT_L5B]<<"  "
					  <<std::setw(5)<<std::setfill(' ') <<avgACount[NT_L5B]);
	UG_LOGN("");
	UG_LOGN("number of exp2 synapses:");
	UG_LOGN("                  min    max    avg");
	UG_LOGN("  L2/3 <- L2/3  "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L23][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L23][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L23][NT_L23]);
	UG_LOGN("  L2/3 <- L4    "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L23][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L23][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L23][NT_L4]);
	UG_LOGN("  L2/3 <- L5A   "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L23][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L23][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L23][NT_L5A]);
	UG_LOGN("  L2/3 <- L5B   "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L23][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L23][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L23][NT_L5B]);
	UG_LOGN("  L4   <- L2/3  "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L4][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L4][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L4][NT_L23]);
	UG_LOGN("  L4   <- L4    "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L4][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L4][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L4][NT_L4]);
	UG_LOGN("  L4   <- L5A   "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L4][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L4][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L4][NT_L5A]);
	UG_LOGN("  L4   <- L5B   "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L4][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L4][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L4][NT_L5B]);
	UG_LOGN("  L5A  <- L2/3  "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L5A][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L5A][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L5A][NT_L23]);
	UG_LOGN("  L5A  <- L4    "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L5A][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L5A][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L5A][NT_L4]);
	UG_LOGN("  L5A  <- L5A   "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L5A][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L5A][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L5A][NT_L5A]);
	UG_LOGN("  L5A  <- L5B   "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L5A][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L5A][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L5A][NT_L5B]);
	UG_LOGN("  L5B  <- L2/3  "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L5B][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L5B][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L5B][NT_L23]);
	UG_LOGN("  L5B  <- L4    "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L5B][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L5B][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L5B][NT_L4]);
	UG_LOGN("  L5B  <- L5A   "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L5B][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L5B][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L5B][NT_L5A]);
	UG_LOGN("  L5B  <- L5B   "<<std::setw(5)<<std::setfill(' ')<<minECount[NT_L5B][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECount[NT_L5B][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECount[NT_L5B][NT_L5B]);
	UG_LOGN("");
	UG_LOGN("number of connected neurons:");
	UG_LOGN("                  min    max    avg");
	UG_LOGN("  L2/3 <- L2/3  "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L23][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L23][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L23][NT_L23]);
	UG_LOGN("  L2/3 <- L4    "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L23][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L23][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L23][NT_L4]);
	UG_LOGN("  L2/3 <- L5A   "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L23][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L23][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L23][NT_L5A]);
	UG_LOGN("  L2/3 <- L5B   "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L23][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L23][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L23][NT_L5B]);
	UG_LOGN("  L4   <- L2/3  "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L4][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L4][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L4][NT_L23]);
	UG_LOGN("  L4   <- L4    "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L4][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L4][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L4][NT_L4]);
	UG_LOGN("  L4   <- L5A   "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L4][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L4][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L4][NT_L5A]);
	UG_LOGN("  L4   <- L5B   "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L4][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L4][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L4][NT_L5B]);
	UG_LOGN("  L5A  <- L2/3  "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L5A][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L5A][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L5A][NT_L23]);
	UG_LOGN("  L5A  <- L4    "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L5A][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L5A][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L5A][NT_L4]);
	UG_LOGN("  L5A  <- L5A   "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L5A][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L5A][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L5A][NT_L5A]);
	UG_LOGN("  L5A  <- L5B   "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L5A][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L5A][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L5A][NT_L5B]);
	UG_LOGN("  L5B  <- L2/3  "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L5B][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L5B][NT_L23]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L5B][NT_L23]);
	UG_LOGN("  L5B  <- L4    "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L5B][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L5B][NT_L4]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L5B][NT_L4]);
	UG_LOGN("  L5B  <- L5A   "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L5B][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L5B][NT_L5A]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L5B][NT_L5A]);
	UG_LOGN("  L5B  <- L5B   "<<std::setw(5)<<std::setfill(' ')<<minECountNeuronWise[NT_L5B][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<maxECountNeuronWise[NT_L5B][NT_L5B]<<"  "
							  <<std::setw(5)<<std::setfill(' ')<<avgECountNeuronWise[NT_L5B][NT_L5B]);
}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
write_activity_to_file(const std::string& fileName, number time)
{
	unsigned long countS[4];
	unsigned long countA[4];
	unsigned long countE[4][4];
	for (size_t i = 0; i < 4; ++i)
	{
		countS[i] = 0;
		countA[i] = 0;
		for (size_t j = 0; j < 4; ++j) countE[i][j] = 0;
	}

	// somata
	size_t sz = m_vSomaVertices.size();
	UG_COND_THROW(sz != m_vNeuronType.size(), "size mismatch");
	for (size_t s = 0; s < sz; ++s)
	{
		number vm = m_spCEDisc->vm(m_vSomaVertices[s]);
		if (vm > -45.0)
		{
			UG_COND_THROW(m_vNeuronType[s]>=4,"invalid neuron type");
			++(countS[m_vNeuronType[s]]);
		}
	}

	// alpha synapses
	sz = m_locASynInfo.size();
	for (size_t a = 0; a < sz; ++a)
	{
		if (is_active(m_locASynInfo[a], time))
		{
			UG_COND_THROW(m_locASynInfo[a].to >= 4,"invalid neuron type");
			countA[m_locASynInfo[a].to] += 1;
		}
	}

	// exp2 synapses
	sz = m_locESynInfo.size();
	for (size_t e = 0; e < sz; ++e)
	{
		if (is_active(m_locESynInfo[e], time))
		{
			UG_COND_THROW(m_locESynInfo[e].to >= 4,"invalid neuron type");
			UG_COND_THROW(m_locESynInfo[e].from >= 4,"invalid neuron type");
			countE[m_locESynInfo[e].to][m_locESynInfo[e].from] += 1;
		}
	}

#ifdef UG_PARALLEL
	// communicate
	unsigned long* locSumData = new unsigned long[24];
	unsigned long* globSumData;

	bool isOutProc = GetLogAssistant().is_output_process();

	if (isOutProc) globSumData = new unsigned long[24];
	else globSumData = NULL;


	memcpy(locSumData, &countS[0], sizeof(unsigned long)*4);
	memcpy(locSumData+4, &countA[0], sizeof(unsigned long)*4);
	for (size_t i = 0; i < 4; ++i)
		memcpy(locSumData+8+4*i, &countE[i][0], sizeof(unsigned long)*4);

	pcl::ProcessCommunicator com;
	com.reduce(locSumData, globSumData, 24, PCL_DT_UNSIGNED_LONG,
			   PCL_RO_SUM, GetLogAssistant().get_output_process());

	delete[] locSumData;

	// check if this proc is output proc
	if (isOutProc)
	{
		// only out proc needs to receive
		memcpy(&countS[0], globSumData, sizeof(unsigned long)*4);
		memcpy(&countA[0], globSumData+4, sizeof(unsigned long)*4);
		for (size_t i = 0; i < 4; ++i)
			memcpy(&countE[i][0], globSumData+8+4*i, sizeof(unsigned long)*4);

		delete[] globSumData;
	// write to file
#endif
		std::string name[4];
		name[0] = std::string("L23");
		name[1] = std::string("L4");
		name[2] = std::string("L5A");
		name[3] = std::string("L5B");

		for (size_t i = 0; i < 4; ++i)
		{
			// somata
			std::string fName = fileName + std::string("_soma_") + name[i];
			std::ofstream outFileSoma;
			outFileSoma.open(fName.c_str(), std::ios_base::app);

			try {outFileSoma << time << "\t" << countS[i] << "\n";}
			UG_CATCH_THROW("Output file" << fName << "could not be written to.");
			outFileSoma.close();

			// alphas
			fName = fileName + std::string("_alpha_") + name[i];
			std::ofstream outFileAlpha;
			outFileAlpha.open(fName.c_str(), std::ios_base::app);

			try {outFileAlpha << time << "\t" << countA[i] << "\n";}
			UG_CATCH_THROW("Output file" << fName << "could not be written to.");
			outFileAlpha.close();

			// exp2s
			for (size_t j = 0; j < 4; ++j)
			{
				fName = fileName + std::string("_exp2_") + name[i] + std::string("_from_") + name[j];
				std::ofstream outFileExp2;
				outFileExp2.open(fName.c_str(), std::ios_base::app);

				try {outFileExp2 << time << "\t" << countE[i][j] << "\n";}
				UG_CATCH_THROW("Output file" << fName << "could not be written to.");
				outFileExp2.close();
			}
		}
#ifdef UG_PARALLEL
	}
#endif

}


template <typename TDomain>
void NETISynapseHandler<TDomain>::
resize_presyn_vector()
{
	// on each processor, get the maximal index, then communicate and resize.

	// loop presynaptic subset
	UG_COND_THROW(!m_spApprox.valid(), "No valid approximation space available.");

	// loop presynaptic subset
	MGSubsetHandler& ssh = *m_spApprox->domain()->subset_handler();

	const std::vector<Vertex*>& sfv = m_spCEDisc->surface_vertices();
	size_t sfv_sz = sfv.size();
	unsigned long max = 0;
	for (size_t sv = 0; sv < sfv_sz; ++sv)
	{
		Vertex* vrt = sfv[sv];
		if (ssh.get_subset_index(vrt) != m_presynSI) continue;

		uint idx = m_aaPSI[vrt];
		//UG_LOG_ALL_PROCS("idx: " << idx << std::endl);
		if (idx > max) max = idx;
	}

	// communicate
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		unsigned long localMax = max;
		com.allreduce(&localMax, &max, 1, PCL_DT_UNSIGNED_LONG, PCL_RO_MAX);
	}
#endif

	// resize
	m_vPresynVmValues.resize(max+1);
}



template <typename TDomain>
bool NETISynapseHandler<TDomain>::
has_EXP2_SYNAPSE()
{
	typedef synapse_traits<> STV;

	bool result = false;

	// Iterate over all base level edges
	for(EdgeIterator eIter = m_spGrid->begin<Edge>(0); eIter != m_spGrid->end<Edge>(0); ++eIter)
	{
		Edge* e = *eIter;
		std::vector<SynapseInfo>& vInfo = m_aaSynapseInfo[e];

		// Iterate over all edge SynapseInfos
		for(size_t i = 0; i < vInfo.size(); ++i)
		{
			SynapseInfo& info = vInfo[i];

			// return true, if SynapseInfo of type EXP2_SYNAPSE was found
			if(STV::type(info) == EXP2_SYNAPSE)
			{
				result = true;
				goto communicate;
			}
		}
	}

	communicate:
#ifdef UG_PARALLEL
	// in parallel, we need to communicate whether _any_ of the procs
	// has a biexp synapse; otherwise we will get problems with global communication
	// in resize_presyn_vector()!
	if (pcl::NumProcs() > 1)
	{
		int thisProcessHasSynapse = result;
		int anyProcessHasSynapse;
		pcl::ProcessCommunicator com;
		com.allreduce(&thisProcessHasSynapse, &anyProcessHasSynapse, 1, PCL_DT_INT, PCL_RO_SUM);
		result = anyProcessHasSynapse != 0;
	}
#endif

	return result;
}


// ////////////////////////////////////
//	explicit template instantiations //
// ////////////////////////////////////

#ifdef UG_DIM_1
	template class NETISynapseHandler<Domain1d>;
#endif

#ifdef UG_DIM_2
	template class NETISynapseHandler<Domain2d>;
#endif

#ifdef UG_DIM_3
	template class NETISynapseHandler<Domain3d>;
#endif



} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug


