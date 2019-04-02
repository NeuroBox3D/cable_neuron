/*
 * simple_cable_neuron_partitioner.h
 *
 *  Created on: 2019-03-18
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__CABLE_NEURON__UTIL__SIMPLE_CABLE_NEURON_PARTITIONER_H
#define UG__PLUGINS__CABLE_NEURON__UTIL__SIMPLE_CABLE_NEURON_PARTITIONER_H


#include "common/math/math_vector_matrix/math_vector.h"  // for MathVector
#include "lib_grid/attachments/attachment_pipe.h"  // for Attachment
#include "lib_grid/multi_grid.h"  // for MultiGrid
#include "lib_grid/parallelization/partitioner.h"  // for IPartitioner
#include "lib_grid/parallelization/process_hierarchy.h"  // for ProcessHierarchy
#include "lib_grid/tools/subset_handler_grid.h"  // for SubsetHandler

#include <cstddef>  // for size_t
#include <vector>


namespace ug {
namespace cable_neuron {



class SimpleCableNeuronPartitioner
: public IPartitioner
{
	public:
		SimpleCableNeuronPartitioner();

		virtual ~SimpleCableNeuronPartitioner();

	// inherited from IPartitioner
	public:
		virtual void set_next_process_hierarchy(SPProcessHierarchy procHierarchy);
		virtual ConstSPProcessHierarchy current_process_hierarchy() const;
		virtual ConstSPProcessHierarchy next_process_hierarchy() const;

		virtual bool supports_repartitioning() const;

		virtual void set_balance_weights(SPBalanceWeights balanceWeights);

		virtual bool partition(size_t baseLvl, size_t elementThreshold);
		virtual SubsetHandler& get_partitions();
		virtual const std::vector<int>* get_process_map() const;

	// needed for DomainPartitioner<Domain3d, SimpleCableNeuronPartitioner>
	public:
		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<3> > aPos);

	private:
		struct NeuronSorting
		{
			NeuronSorting(const std::vector<size_t>& ec)
			: edgeCounts(ec) {}

			bool operator()(const size_t& a, const size_t& b)
			{
				return edgeCounts[a] < edgeCounts[b];
			}

			const std::vector<size_t>& edgeCounts;
		};

	protected:
		MultiGrid* m_mg;

		SPProcessHierarchy m_spProcessHierarchy;
		SPProcessHierarchy m_spNextProcessHierarchy;

		SubsetHandler m_partitionSH;
};



} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__UTIL__SIMPLE_CABLE_NEURON_PARTITIONER_H
