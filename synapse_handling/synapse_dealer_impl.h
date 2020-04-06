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
    if (shClass3)  // This may be null, e.g., if we do not compile for 3d
    {
    	shClass3->add_method("begin_" + name, &TSH::template begin_wrapper<TSyn>);
    	shClass3->add_method("end_" + name, &TSH::template end_wrapper<TSyn>);
    }

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
