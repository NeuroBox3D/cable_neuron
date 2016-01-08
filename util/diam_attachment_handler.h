/*
 * diam_attachment_handler.h
 *
 *  Created on: 21.05.2015
 *      Author: mbreit
 */
#ifndef __UG__PLUGINS__CABLE_NEURON__UTIL__DIAM_ATTACHMENT_HANDLER_H__
#define __UG__PLUGINS__CABLE_NEURON__UTIL__DIAM_ATTACHMENT_HANDLER_H__


#include "lib_grid/lib_grid.h"
#include "lib_grid/tools/copy_attachment_handler.h"

#include <vector>

namespace ug {
namespace cable_neuron {

/**
 * @brief handler for diameter attachment in the CableEquation multi-grid
 *
 * This class implements an attachment handler for the diameter attachment used
 * in the CableEquation class.
 *
 * The diameter vertex attachment needs to be propagated not only by copying from
 * base vertices to their child vertices and grand-child vertices (and so on);
 * it is also required on vertices that appear as children of refined edges.
 *
 * In these cases, the diameter needs to be calculated. The current implementation
 * uses the average value of the two vertices bounding the parent edge.
 *
 */
class DiamAttachmentHandler
: public CopyAttachmentHandler<Vertex, ANumber>
{
	public:

		/// constructor
		DiamAttachmentHandler() {};

		/// destructor
		virtual ~DiamAttachmentHandler() {};

	protected:
		virtual void copy_from_other_elem_type(GridObject* parent, Vertex* child);
};

} // end namespace cable_neuron
} // end namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__UTIL__DIAM_ATTACHMENT_HANDLER_H__
