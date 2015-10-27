/*
 * diam_attachment_handler.h
 *
 *  Created on: 21.05.2015
 *      Author: mbreit
 */
#ifndef PLUGINS__HHCABLE__DIAM_ATTACHMENT_HANDLER_H_
#define PLUGINS__HHCABLE__DIAM_ATTACHMENT_HANDLER_H_


#include "lib_grid/lib_grid.h"
#include "lib_grid/tools/copy_attachment_handler.h"

#include <vector>

namespace ug {
namespace cable {

/**
 * @brief handler for diameter attachment in the VMDisc multi-grid
 *
 * This class implements an attachment handler for the diameter attachment used
 * in the VMDisc class.
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

} // end namespace cable
} // end namespace ug

#endif // PLUGINS__HHCABLE__DIAM_ATTACHMENT_HANDLER_H_
