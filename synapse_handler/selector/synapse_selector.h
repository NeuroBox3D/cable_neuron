/*!
 * \file synapse_selector.h
 * \brief synapse selector
 *
 *  Created on: Apr 17, 2015
 *      Author: stephan
 */

/// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SELECTOR__SYNAPSE_SELECTOR_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SELECTOR__SYNAPSE_SELECTOR_H__

/// includes
#include <map>

#include <common/error.h>

#include "synapse_selector_traits.h"
#include "synapse_selector_info_traits.h"

#include "../function/types.h"

/* \defgroup sh_plugin Synapse Handler plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class SynapseSelector {
	/// public functions
	public:
	/*!
	 * \brief get the current of the associated synapse
	 */
	template <class TSynapse>
	static number current(TSynapse& synapse, const number& t) {
		return synapse_selector_trait<TSynapse>::current(synapse, t);
	}

	/*!
	 * \brief get the name of the type
	 */
	template <class TSynapse>
	static SynapseType name(TSynapse& synapse) {
		return synapse_selector_info_traits<TSynapse>::type_name();
	}


	/*!
	 * \brief synapse
	 */
	struct SynapseEntry {
		ISynapse*	synapse;
		SynapseType		type;

		SynapseEntry () : synapse(NULL), type(EMPTY_SYNAPSE) {
		}

		SynapseEntry (SynapseType& typ, ISynapse* syn) : synapse(syn), type(typ) {
		}

	};

	typedef std::map<SynapseType, SynapseEntry> SynapsesMap;
	std::vector<SynapseType>	m_synapses;
	SynapsesMap				m_synapsesMap;

	/*!
	 * \brief declare a synapse to be used
	 */
	template <class TSynapse>
	static void declare_synapse (const SynapseType& name) {
		SynapseType typeName = synapse_selector_info_traits<TSynapse>::type_name();
		if(synapses()[name].synapse != NULL) {
				UG_COND_THROW(
					dynamic_cast<TSynapse*> (synapses()[name].synapse) == NULL,
						  "Synapse with name '" << name
						  << "' was already declared in SynapseSelector with a different type. "
						  << "Old type: " << synapses()[name].type <<
						  ", new type: " << typeName);
				return;
		}

		synapses_names().push_back(name);
		synapses()[name] = SynapseEntry(typeName, new TSynapse());
	}

	/*!
	 * \brief return all names of declared synapses
	 */
	static const std::vector<SynapseType>& declared_synapses_names () {
		return synapses_names();
	}

	/*!
	 * \brief get the declared synapses as a map
	 */
	static SynapsesMap& synapses() {
		return inst().m_synapsesMap;
	}

	/*!
	 * \brief determine if the synapse with name was declared
	 */
	static bool is_declared(const SynapseType& name) {
		return synapses().find(name) != synapses().end();
	}


	/*
	 * \brief query synapse with given name and return the synapse
	 */
	template <class TSynapse>
	static TSynapse synapse (const SynapseType& name) {
		SynapseEntry& e = synapse_entry(name);
		TSynapse* a = dynamic_cast<TSynapse*>(e.synapse);
		UG_COND_THROW(!a, "Synapse with invalid type queried. Given type "
					  "is " << e.type << ", queried type is " <<
					  synapse_selector_info_traits<TSynapse>::type_name());
		return *a;
	}

	/*!
	 * \brief get the type name of given synapse name
	 */
	static SynapseType type_name (const SynapseType& name) {
		SynapseEntry& e = synapse_entry(name);
		return e.type;
	}

	static SynapseSelector& inst () {
		static SynapseSelector h;
		return h;
	}

	SynapseSelector() {
	}

	~SynapseSelector() {
	}

	static std::vector<SynapseType>& synapses_names() {
		return inst().m_synapses;
	}

	static SynapseEntry& synapse_entry(const SynapseType& name) {
		SynapseEntry& e = synapses()[name];
		UG_COND_THROW(!e.synapse, "Undeclared synapse queried: " << name);
		return e;
	}
};


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
//<! \}

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SELECTOR__SYNAPSE_SELECTOR_H__
