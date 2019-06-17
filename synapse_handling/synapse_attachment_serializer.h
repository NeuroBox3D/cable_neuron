/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-12-02
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

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_SERIALIZER_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_SERIALIZER_H

#include "lib_grid/algorithms/serialization.h"   // GeomObjDataSerializer
#include "synapses/base_synapse.h"  // IBaseSynapse

#include <vector>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {


/**
 * As the synapse attachment is a vector only containing pointers
 * to the actual synapse objects, we need to manually take care of communicating
 * these objects during domain decomposition and redistribution.
 */
class SynapseAttachmentSerializer
: public GeomObjDataSerializer<Edge>
{
	public:
		typedef Attachment<std::vector<IBaseSynapse*> > SynAttach;

	public:
		SynapseAttachmentSerializer(Grid& g, SynAttach& a);

		virtual ~SynapseAttachmentSerializer();

		/**
		 * @brief Write synapses of a specific edge to binary buffer.
		 * When synapse is written to buffer it is guaranteed to be moved to another process.
		 * This is why we need to delete it here via its pointer on the current proc.
		 */
		virtual void write_data(BinaryBuffer& out, Edge* o) const;

		/**
		 * @brief Read synapses of a specific edge from binary buffer.
		 * The data received in the buffer will trigger the creation of a new synapse
		 * and storing of a pointer to it in the attachment.
		 */
		virtual void read_data(BinaryBuffer& in, Edge* o);

	protected:
		Grid::AttachmentAccessor<Edge, SynAttach> m_aa;

		bool m_blubb;
};


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_SERIALIZER_H
