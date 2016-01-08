/*
 * util.h
 *
 *	Utilities for the cable project.
 *
 *  Created on: 08.08.2015
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__UTIL__CABLE_ASS_TUNER_H__
#define __UG__PLUGINS__CABLE_NEURON__UTIL__CABLE_ASS_TUNER_H__

#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_grid/tools/selector_grid.h"


namespace ug {
namespace cable_neuron {


template <typename TDomain, typename TAlgebra>
class CableAssTuner
{
	public:
		CableAssTuner
		(
			SmartPtr<DomainDiscretization<TDomain, TAlgebra> > domDisc,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		template <size_t dim>
		void remove_ghosts_from_assembling_iterator();

	private:
		SmartPtr<DomainDiscretization<TDomain, TAlgebra> > m_spDomDisc;
		SmartPtr<ApproximationSpace<TDomain> > m_spApprox;
		SmartPtr<Selector> m_sel;


#if 0

template <typename TDomain, typename TAlgebra>
class CableAssTuner
{
	public:
		CableAssTuner
		(
			SmartPtr<DomainDiscretization<TDomain, TAlgebra> > domDisc,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		template <size_t dim>
		void remove_ghosts_from_assembling_iterator();

		void update_elem_lists(uint elem_types);
		const std::vector<Vertex*>& vertex_vector() const;
		const std::vector<Edge*>& edge_vector() const;

	private:
		template <uint dim>
		void update_elem_list_for_elem();

		template <typename TElem, typename TDummy = void>
		struct StoreElem
		{
			StoreElem(TElem* elem, CableAssTuner* cas)
				{UG_THROW("Not implemented for this elem type.");}
		};
		template <typename TDummy>
		struct StoreElem<Vertex, TDummy>
		{
			StoreElem(Vertex* elem, CableAssTuner* cas)
				{cas->m_vVrt.push_back(elem);}
		};
		template <typename TDummy>
		struct StoreElem<Edge, TDummy>
		{
			StoreElem(Edge* elem, CableAssTuner* cas)
				{cas->m_vEdge.push_back(elem);}
		};
		template <typename TElem, typename TDummy>
		friend struct StoreElem;

	private:
		SmartPtr<DomainDiscretization<TDomain, TAlgebra> > m_spDomDisc;
		SmartPtr<ApproximationSpace<TDomain> > m_spApprox;

		bool m_bRemoveGhostsFromAssemblingIterator;
		SmartPtr<Selector> m_sel;

		std::vector<Vertex*> m_vVrt;
		std::vector<Edge*> m_vEdge;

#endif

};

} // namespace cable_neuron
} // namespace ug


#include "cable_ass_tuner_impl.h"


#endif // __UG__PLUGINS__CABLE_NEURON__UTIL__CABLE_ASS_TUNER_H__


