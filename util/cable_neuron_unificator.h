/*
 * cable_neuron_unificator.h
 *
 *  Created on: 2019-03-22
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__CABLE_NEURON__UTIL__CABLE_NEURON_UNIFICATOR_H
#define UG__PLUGINS__CABLE_NEURON__UTIL__CABLE_NEURON_UNIFICATOR_H


#include "../../Parmetis/src/unificator_interface.h"  // for IUnificator


namespace ug {
namespace cable_neuron {


class CableNeuronUnificator
: public parmetis::IUnificator<Edge>
{
	public:
		CableNeuronUnificator();
		virtual ~CableNeuronUnificator();

		typedef Edge::side side_t;
		typedef Attachment<int> AElemIndex;
		typedef Attachment<std::vector<int> > AElemIndices;

		// inherited from IUnificator
		virtual void unify
		(
			MultiGrid* mg,
			int lvl,
			int localOffset,
			const Grid::AttachmentAccessor<Edge, AElemIndex>& aaElemInd, // local indices!
			const Grid::AttachmentAccessor<side_t, AElemIndices>& aaSideElemInd, // global indices!
			std::vector<std::pair<int, int> >& unificationPairs // global indices!
		) const;
};


}  // namespace cable_neuron
}  // namespace ug

#endif  // UG__PLUGINS__CABLE_NEURON__UTIL__CABLE_NEURON_UNIFICATOR_H
