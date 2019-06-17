/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Lukas Reinhardt
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

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_INFO_IO_TRAITS_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_INFO_IO_TRAITS_H

#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include "lib_grid/attachments/attachment_info_traits.h"
#include "lib_grid/attachments/attachment_io_traits.h"
#include "lib_grid/global_attachments.h"

#include "types.h"
#include "synapses/base_synapse.h"
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
//			for(size_t j = 0; j < identifier.size(); ++j) {
//				in.unget();										//reset input stream
//			}

			IBaseSynapse* s;
			try {s = SynapseDealer::instance()->deal(identifier);}
			UG_CATCH_THROW("Synapse of type '" << identifier << "' could not be created.")

			in >> s;
			v[i] = s;

		}
	}
};

//!< declare the attachment
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<std::vector<IBaseSynapse*> >, "IBaseSynapse*");


} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_INFO_IO_TRAITS_H
