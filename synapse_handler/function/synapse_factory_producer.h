/*!
 * \file synapse_factory_producer.h
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__SYNAPSE_FACTORY_PRODUCER_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__SYNAPSE_FACTORY_PRODUCER_H__

/// includes
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/assign/list_of.hpp>

#include <common/error.h>

#include "../synapse.h"
#include "../synapse_factory.h"

/* \defgroup sh_plugin Synapse Provider plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace cable_neuron {
namespace synapse_handler {

/*!
 * \brief synapse factory producer
 */
class SynapseFactoryProducer {
private:
	/// synapse factory types
	static const std::string ALPHA_SYNAPSE_FACTORY;
	static const std::string EXP_SYNAPSE_FACTORY;
	static const std::string EMPTY_SYNAPSE_FACTORY;
	static const std::string JANA_SYNAPSE_FROM_MARKUS_WITH_LOVE_FACTORY;

public:
	/*!
	* \brief get factory
	*/
	SmartPtr<ISynapseFactory> getSynapseFactory(const std::string& description) {
		if (boost::iequals(ALPHA_SYNAPSE_FACTORY, description)) {
			return make_sp<AlphaSynapseFactory>(new AlphaSynapseFactory());
		} else if (boost::iequals(EXP_SYNAPSE_FACTORY, description)) {
			return make_sp<Exp2SynapseFactory>(new Exp2SynapseFactory());
		} else {
			UG_THROW("An unknown synapse factory was selected!");
		}
		return SPNULL;
	}

	/*
	 * \brief get default factory
	 */
	SmartPtr<ISynapseFactory> get_default_synapse_factory() {
		return make_sp(new AlphaSynapseFactory());
	}

	/*!
	 * \brief get alpha synapse factory
	 */
	SmartPtr<AlphaSynapseFactory> getAlphaSynapseFactory() {
		return make_sp<AlphaSynapseFactory>(new AlphaSynapseFactory());
	}

	/*!
	 * \brief get bi-exponential synapse factory
	 */
	SmartPtr<Exp2SynapseFactory> getExp2SynapseFactory() {
		return make_sp<Exp2SynapseFactory>(new Exp2SynapseFactory());
	}
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
//<! \}

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__SYNAPSE_FACTORY_PRODUCER_H__
