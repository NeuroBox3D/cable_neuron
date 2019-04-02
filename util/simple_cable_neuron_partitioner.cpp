/*
 * simple_cable_neuron_partitioner.cpp
 *
 *  Created on: 2019-03-18
 *      Author: mbreit
 */

#include "simple_cable_neuron_partitioner.h"
#include "functions.h"  // for neuron_identification

#include "lib_grid/global_attachments.h"  // for global attachments


namespace ug {
namespace cable_neuron {


SimpleCableNeuronPartitioner::SimpleCableNeuronPartitioner()
: m_mg(NULL)
{
	m_spProcessHierarchy = SPProcessHierarchy(new ProcessHierarchy);
	m_spProcessHierarchy->add_hierarchy_level(0, 1);
}


SimpleCableNeuronPartitioner::~SimpleCableNeuronPartitioner()
{}


void SimpleCableNeuronPartitioner::set_next_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_spNextProcessHierarchy = procHierarchy;
}


ConstSPProcessHierarchy SimpleCableNeuronPartitioner::current_process_hierarchy() const
{
	return m_spProcessHierarchy;
}


ConstSPProcessHierarchy SimpleCableNeuronPartitioner::next_process_hierarchy() const
{
	return m_spNextProcessHierarchy;
}



bool SimpleCableNeuronPartitioner::supports_repartitioning() const
{
	// This implementation only handles the initial partitioning.
	return false;
}


void SimpleCableNeuronPartitioner::set_balance_weights(SPBalanceWeights balanceWeights)
{
	// ignore
}


bool SimpleCableNeuronPartitioner::partition(size_t baseLvl, size_t elementThreshold)
{
	// the partitioner only works on the base level and cannot handle refined grids
	UG_COND_THROW(m_mg->top_level() != 0,
		"The SimpleCableNeuronPartitioner only works on unrefined grids.");

	m_partitionSH.clear();

	const ProcessHierarchy* procH;
	if(m_spNextProcessHierarchy.valid())
		procH = m_spNextProcessHierarchy.get();
	else
		procH = m_spProcessHierarchy.get();

	// make sure every vertex has a neuron ID
	neuron_identification(*m_mg);

    typedef Attachment<uint> ANeuronID;
    if (!GlobalAttachments::is_declared("neuronID"))
        UG_THROW("GlobalAttachment 'NeuronID' not declared.");
	ANeuronID aNID = GlobalAttachments::attachment<ANeuronID>("neuronID");
	Grid::VertexAttachmentAccessor<ANeuronID> aaNID =
		Grid::VertexAttachmentAccessor<ANeuronID>(*m_mg, aNID);


	const size_t nHLvls = procH->num_hierarchy_levels();
	for (size_t hlvl = 0; hlvl < nHLvls; ++hlvl)
	{
		// do not distribute if hierarchy level starts on grid level > 0
		if (procH->grid_base_level(hlvl) > m_mg->top_level())
			break;

		// now count neurons and elements of each neuron on min level
		std::vector<size_t> edgeCounts;
		typedef geometry_traits<Edge>::const_iterator const_edge_iter;
		const_edge_iter it = m_mg->begin<Edge>(0);
		const_edge_iter itEnd = m_mg->end<Edge>(0);
		for (; it != itEnd; ++it)
		{
			const uint nid = aaNID[(*it)->vertex(0)];
			if ((size_t) nid >= edgeCounts.size())
				edgeCounts.resize(nid+1, 0);
			++edgeCounts[nid];
		}

		// if there are N neurons and P processors
		// every proc will either get N/P or N/P + 1 neurons
		const size_t nNeurons = edgeCounts.size();
		const size_t nProcs = procH->num_global_procs_involved(hlvl);

		const size_t minNeuronsPerProc = nNeurons / nProcs;
		const size_t nProcsWithMaxNeurons = nNeurons % nProcs;

		// sort neurons w.r.t. their edge counts
		std::vector<size_t> sorting(nNeurons);
		for (size_t i = 0; i < nNeurons; ++i)
			sorting[i] = i;
		std::sort(sorting.begin(), sorting.end(), NeuronSorting(edgeCounts));

		// first assign the smallest neurons to the procs that receive N/P+1 ones
		std::vector<size_t> assignedProcs(nNeurons);
		size_t cur = 0;
		for (size_t i = 0; i < nProcsWithMaxNeurons; ++i)
		{
			for (size_t j = 0; j < minNeuronsPerProc + 1; ++j)
			{
				assignedProcs[sorting[cur]] = i;
				++cur;
			}
		}

		// then assign the rest
		for (size_t i = nProcsWithMaxNeurons; i < nProcs; ++i)
		{
			for (size_t j = 0; j < minNeuronsPerProc; ++j)
			{
				assignedProcs[sorting[cur]] = i;
				++cur;
			}
		}

		// now assign partitions to each
		for (it = m_mg->begin<Edge>(0); it != itEnd; ++it)
		{
			Edge* e = *it;
			const uint nid = aaNID[e->vertex(0)];
			m_partitionSH.assign_subset(e, assignedProcs[nid]);
		}
	}

	return true;
}


SubsetHandler& SimpleCableNeuronPartitioner::get_partitions()
{
	return m_partitionSH;
}


const std::vector<int>* SimpleCableNeuronPartitioner::get_process_map() const
{
	return NULL;
}



void SimpleCableNeuronPartitioner::set_grid(MultiGrid* mg, Attachment<MathVector<3> > aPos)
{
	if (mg == m_mg)
		return;

	if (m_mg)
	{
		m_partitionSH.assign_grid(NULL);
		m_mg = NULL;
	}

	if (mg)
	{
		m_mg = mg;
		m_partitionSH.assign_grid(m_mg);
	}
}


} // namespace cable_neuron
} // namespace ug
