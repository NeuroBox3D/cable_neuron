/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Lukas Reinhardt
 * Creation date: 2016-03-22
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

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_DEALER_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_DEALER_H

#include <map>
#include <vector>
#include <string>

#include "registry/registry.h"
#include "synapses/base_synapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class ISynapseGenerator
{
public:
	virtual IBaseSynapse* generate() = 0;
	virtual ~ISynapseGenerator() {}
};

template<typename TSyn>
class SynapseGenerator : public ISynapseGenerator
{
public:
	IBaseSynapse* generate() {return new TSyn();}
};

class SynapseDealer {
	private:
		std::map<std::string, size_t> m_synReg;
        std::vector<ISynapseGenerator*> m_vSynGen;
		std::vector<size_t> m_vSize;
		std::vector<std::vector<size_t> > m_vNoOverwriteBytes;
		static SynapseDealer *m_instance;
		~SynapseDealer();
		SynapseDealer();

	public:
        /// register pre-synpase at synapse type registry and UG registry
        template <typename TSyn>
        bridge::ExportedClass<TSyn>& register_pre_synapse_type
        (
            bridge::Registry* reg,
            const std::string& grp
        );

        /// register post-synpase at synapse type registry and UG registry
        template <typename TSyn>
        bridge::ExportedClass<TSyn>& register_post_synapse_type
        (
            bridge::Registry* reg,
            const std::string& grp
        );

        /// register only at synapse type registry, but not at UG registry
        template <typename TSyn>
        void register_synapse_type();

		size_t unique_id(const std::string& name);

		IBaseSynapse* deal(const std::string& t);
		IBaseSynapse* deal(size_t unique_id);

		size_t size_of(const std::string& name);
		size_t size_of(size_t unique_id);

		const std::vector<size_t>& non_data_bytes(size_t uid) const;

		static SynapseDealer* instance() {
			if(m_instance == 0) {
				m_instance = new SynapseDealer;
			}
			return m_instance;
		}

	protected:
		template <typename TSyn, typename SynSubtype>
        bridge::ExportedClass<TSyn>& register_synapse_type
        (
            const std::string& name,
            bridge::Registry* reg,
            const std::string& grp
        );

		size_t get_unique_id();

		void find_no_overwrite_bytes
		(
			const char* s,
			size_t uid
		);
};


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug


#include "synapse_dealer_impl.h"

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_DEALER_H
