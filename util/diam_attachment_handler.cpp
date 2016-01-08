/*
 * diam_attachment_handler.cpp
 *
 *  Created on: 21.05.2015
 *      Author: mbreit
 */


#include "diam_attachment_handler.h"

namespace ug {
namespace cable_neuron {


void DiamAttachmentHandler::copy_from_other_elem_type(GridObject* parent, Vertex* child)
{
	// ensure that parent is an edge
	Edge* par = dynamic_cast<Edge*>(parent);
	if (!par) return;

	// average parent edge vertex diameters
	number diam = 0.0;
	size_t num_vrt = par->num_vertices();

	for (size_t j = 0; j < num_vrt; ++j)
		diam += m_aa[(*par)[j]];

	m_aa[child] = diam / num_vrt;
}

} // end namespace cable_neuron
} // end namespace ug

