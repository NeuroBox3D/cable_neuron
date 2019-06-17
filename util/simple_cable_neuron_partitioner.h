/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-03-18
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
