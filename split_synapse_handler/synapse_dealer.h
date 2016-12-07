/*
 * SynapseDealer.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SYNAPSE_DEALER_H_
#define SPLIT_SYNAPSE_HANDLER_SYNAPSE_DEALER_H_

#include <map>
#include <string>

#include "base_synapse.h"

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
		std::map<std::string, size_t> m_register;
		std::vector<ISynapseGenerator*> m_vSynGen;
		std::vector<size_t> m_vSize;
		static SynapseDealer *m_instance;
		~SynapseDealer();
		SynapseDealer();

	public:
		template <typename TSyn> void register_synapse_type(std::string name);

		size_t unique_id(std::string name);

		IBaseSynapse* deal(std::string t);
		IBaseSynapse* deal(size_t unique_id);

		size_t size_of(std::string name);
		size_t size_of(size_t unique_id);

		static SynapseDealer* instance() {
			if(m_instance == 0) {
				m_instance = new SynapseDealer;
			}
			return m_instance;
		}

	protected:
		size_t get_unique_id();
};


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */


#include "synapse_dealer_impl.h"

#endif /* SPLIT_SYNAPSE_HANDLER_SYNAPSE_DEALER_H_ */
