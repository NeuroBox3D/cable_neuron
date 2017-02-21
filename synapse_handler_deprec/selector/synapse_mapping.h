/*!
 * \file synapse_mapping.h
 * \brief synapse mappings
 *
 *  Created on: Apr 22, 2015
 *      Author: stephan
 */

/// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SELECTOR__SYNAPSE_MAPPING_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SELECTOR__SYNAPSE_MAPPING_H__

/// includes
#include <common/error.h>
#include "../function/types.h"

/* \defgroup sh_plugin Synapse Handler plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace cable_neuron {
namespace synapse_handler {

/*!
 * \brief a s ynapse entry
 */
struct entry {
	const char* m_presence;
	const char* m_synapses;

	entry(const char* presence, const char* synapses)
	: m_presence(presence), m_synapses(synapses) { };

	entry() : m_presence(""), m_synapses("") { };
 };

/*!
 * \brief synapse mapping (synapse type to attachment types)
 */
class SynapseMapping {
public:
		/// get mapping for string
		entry get_attachment_names(SynapseType type) {
				std::map<SynapseType, entry>::
				const_iterator cit = m_mapping.find(type);

				UG_COND_THROW(cit == m_mapping.end(), "Non-registered synapse mapping to attachment queried.");
				return m_mapping[type];
		}


		/// get singleton instance
		static SynapseMapping& getInstance() {
				static SynapseMapping instance;
				return instance;
		}

private:
		/// stores mapping as follow: declared synapse selector ->
		/// declared attachment for synapse information on edge,
		/// declared attachment for synapse presence at vertex.
		std::map<SynapseType, entry> m_mapping;

		/// register all synapses types with mapping in correspondence
		/// to used attachments in serialized grid files (ugx)
		SynapseMapping() {
		   m_mapping.insert(m_mapping.begin(), std::make_pair(ALPHA_SYNAPSE, entry("synapses-presence-NETI", "synapses-NETI")));
		   m_mapping.insert(m_mapping.begin(), std::make_pair(EXP2_SYNAPSE, entry("synapses-presence-NETI", "synapses-NETI")));
		   m_mapping.insert(m_mapping.begin(), std::make_pair(EMPTY_SYNAPSE, entry("synapses-presence-SD", "synapses-SD")));
		};

		/// private (Scott Meyer's singleton)
		SynapseMapping(const SynapseMapping&);
		SynapseMapping& operator=(const SynapseMapping&);
};


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
//<! \}
#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SELECTOR__SYNAPSE_MAPPING_H__
