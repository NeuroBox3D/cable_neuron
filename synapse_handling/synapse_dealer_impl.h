/*
 * SynapseDealer_impl.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#include "common/error.h"
#include "common/log.h"
#include "synapse_handler.h"
#include <cstring>  // memset

namespace ug {
namespace cable_neuron {
namespace synapse_handler {


template <typename TSyn>
void SynapseDealer::register_synapse_type()
{
	// assign a unique id to each synapse type;
	// save corresponding generators and object sizes
	size_t uid = get_unique_id();
	ISynapseGenerator* sg = new SynapseGenerator<TSyn>();
	m_synReg[TSyn::name_string] = uid;
	m_vSynGen.resize(uid+1);
	m_vSynGen[uid] = sg;
	m_vSize.resize(uid+1);
	m_vSize[uid] = sizeof(TSyn);

	// compute and save byte positions not to be overwritten
	// when communicating this object type from one proc to another
	m_vNoOverwriteBytes.resize(uid+1);

	#ifdef UG_PARALLEL
	char* buf = new char[sizeof(TSyn)];
	std::memset(buf, 0, sizeof(TSyn));
	new (buf) TSyn();
	find_no_overwrite_bytes(buf, uid);
	delete[] buf;
	#endif

	/* DEBUGGING
	UG_LOGN("Current synapse registry");
	for (std::map<std::string, size_t>::const_iterator it = m_register.begin(); it != m_register.end(); ++it)
		UG_LOGN(it->first << ": " << it->second);
	UG_LOGN("");
	*/
}


template <typename TSyn, typename SynSubtype>
bridge::ExportedClass<TSyn>& SynapseDealer::register_synapse_type
(
    const std::string& name,
    bridge::Registry* reg,
    const std::string& grp
)
{
    // check if name already exists and throw if that is the case
    std::map<std::string, size_t>::const_iterator it = m_synReg.find(name);
    UG_COND_THROW(it != m_synReg.end(), "Tried to register synapse type named "
                  "'" << name << "', but this name already exists.");


    // register class at ug registry
    bridge::ExportedClass<TSyn>& expClass = reg->add_class_<TSyn, SynSubtype>(name, grp);

    // register a corresponding iterator class
    typedef SynapseIter<TSyn> TIter;
    reg->add_class_<TIter>(name + "_It", grp)
        .add_method("get", &TIter::operator*)
        .add_method("next", &TIter::next)
        .add_method("inequal", &TIter::inequal);

    // TODO: Technically, this is required for all domains (but we use only 3d).
    // register a begin() and end() method for SynapseHandler
    // returning a corresponding iterator
    typedef SynapseHandler<Domain3d> TSH;
    bridge::ExportedClass<TSH>* shClass3;
    shClass3 = dynamic_cast<bridge::ExportedClass<TSH>*>(reg->get_class("SynapseHandler3d"));
    UG_COND_THROW(!shClass3, "No class registered as SynapseHandler3d could be found in registry.");
    shClass3->add_method("begin_" + name, &TSH::template begin_wrapper<TSyn>);
    shClass3->add_method("end_" + name, &TSH::template end_wrapper<TSyn>);

    // register with this class' synapse "registry"
    register_synapse_type<TSyn>();

	return expClass;
}

template <typename TSyn>
bridge::ExportedClass<TSyn>& SynapseDealer::register_pre_synapse_type
(
   bridge::Registry* reg,
   const std::string& grp
)
{
    return register_synapse_type<TSyn, IPreSynapse>(TSyn::name_string, reg, grp);
}

template <typename TSyn>
bridge::ExportedClass<TSyn>& SynapseDealer::register_post_synapse_type
(
   bridge::Registry* reg,
   const std::string& grp
)
{
    return register_synapse_type<TSyn, IPostSynapse>(TSyn::name_string, reg, grp);
}


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
