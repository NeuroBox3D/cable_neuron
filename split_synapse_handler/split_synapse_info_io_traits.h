/*
 * split_synapse_info_io_traits.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_INFO_IO_TRAITS_H_
#define SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_INFO_IO_TRAITS_H_

#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include "lib_grid/attachments/attachment_info_traits.h"
#include "lib_grid/attachments/attachment_io_traits.h"
#include "lib_grid/global_attachments.h"

#include "../synapse_handler/function/types.h"
#include "base_synapse.h"
#include "synapse_dealer.h"

namespace ug {
using namespace cable_neuron::synapse_handler;


/*!
 * \brief attachment_io_traits specialization for std::vector<SynapseInfo>
 */
template <>
struct attachment_io_traits<Attachment<std::vector<IBaseSynapse*> > > {
	/// attachment type
	typedef std::vector<IBaseSynapse*> value_type;

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

		for (size_t i = 0; i < numSyn; ++i) {
			// create synapse
			std::string identifier;
			in >> identifier;
			for(size_t j = 0; j < identifier.size(); ++j) {
				in.unget();										//reset input stream
			}

			IBaseSynapse* s = SynapseDealer::instance()->deal(identifier); //todo
			in >> s;
			v[i] = s;

		}
	}
};

/*!
 * \brief info for the attachment used in the read and write traits
 */
template <>
struct attachment_info_traits<std::vector<IBaseSynapse*> > {
	static const char* type_name () {
		return "IBaseSynapse*";
	}
};

//!< declare the attachment
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<std::vector<IBaseSynapse*> >, "IBaseSynapse*");

} // namespace ug

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_INFO_IO_TRAITS_H_ */
