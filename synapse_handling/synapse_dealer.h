/*
 * SynapseDealer.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
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
