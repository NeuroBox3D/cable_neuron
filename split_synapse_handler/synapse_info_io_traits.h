/*!
 * \file synapse_io_traits.h
 * \brief read, write and info traits for struct synapse info
 *
 *  Created on: May 11, 2015
 *      Author: stephan
 */

/// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__GRID__SYNAPSE_INFO_IO_TRAITS_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__GRID__SYNAPSE_INFO_IO_TRAITS_H__

/// includes
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>

#include "lib_grid/attachments/attachment_info_traits.h"
#include "lib_grid/attachments/attachment_io_traits.h"
#include "lib_grid/global_attachments.h"

#include "../function/types.h"
#include "synapse_info.h"

namespace ug {
using namespace cable_neuron::synapse_handler;

/*!
 * \brief attachment_io_traits specialization for std::vector<SynapseInfo>
 */
template <>
struct attachment_io_traits<Attachment<std::vector<SynapseInfo> > > {
	/// attachment type
	typedef std::vector<ISuperSynapse*> value_type;

	/// serialize
	static void write_value (std::ostream& out, const value_type& v) {
		out << v.size() << " ";
		for (value_type::const_iterator cit = v.begin();
				cit != v.end(); ++cit) out << *cit << " ";
	}

	/// deserialize
	static void read_value(std::istream& in, value_type& v) {
		// determine number of synapses to be read and reserve
		size_t numSyn;
		in >> numSyn;
		v.resize(numSyn);

		for (size_t i = 0; i < numSyn; ++i)
		{
			// get synapse identifier (e.g. "SynapseA")

			// have SynapseDealer deal a SynapseA*
			ISuperSynapse iss = synapseDealer_singleton.deal("SynapseA");

			// have iss deal with the rest of the input
			in >> (*iss);

			v[i] = iss;
		}
	}
};

/*!
 * \brief info for the attachment used in the read and write traits
 */
template <>
struct attachment_info_traits<std::vector<SynapseInfo> > {
	static const char* type_name () {
		return "SynapseInfo";
	}
};

//!< declare the attachment
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<std::vector<SynapseInfo> >, "SynapseInfo");

} // namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__GRID__SYNAPSE_INFO_IO_TRAITS_H__
