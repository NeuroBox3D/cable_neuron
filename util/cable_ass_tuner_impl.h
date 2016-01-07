/*
 * util_impl.h
 *
 *	Utilities for the cable project.
 *
 *  Created on: 08.08.2015
 *      Author: mbreit
 */

#include "cable_ass_tuner.h"
#include "lib_disc/domain_traits.h"
#include "lib_grid/tools/grid_level.h"
#include "lib_grid/grid/grid_base_objects.h"


namespace ug {
namespace cable {


template <typename TDomain, typename TAlgebra>
CableAssTuner<TDomain, TAlgebra>::
CableAssTuner
(
	SmartPtr<DomainDiscretization<TDomain, TAlgebra> > domDisc,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
:  m_spDomDisc(domDisc), m_spApprox(approx), m_sel(SPNULL)
{}

template <typename TDomain, typename TAlgebra>
template <size_t dim>
void CableAssTuner<TDomain, TAlgebra>::
remove_ghosts_from_assembling_iterator()
{
	typedef typename domain_traits<dim>::grid_base_object elem_type;
	typedef typename DoFDistribution::traits<elem_type>::const_iterator it_type;

	// create selector to store assemble elements
	SmartPtr<MultiGrid> mg = m_spApprox->domain()->grid();
	m_sel = make_sp(new Selector(1 << dim));
	m_sel->assign_grid(mg.get());
	m_sel->clear();

	// select elements
	GridLevel gl(GridLevel::TOP, GridLevel::SURFACE, false);
	ConstSmartPtr<DoFDistribution> dd = m_spApprox->dof_distribution(gl, false);
	it_type it = dd->template begin<elem_type>(SurfaceView::MG_ALL);
	it_type it_end = dd->template end<elem_type>(SurfaceView::MG_ALL);
	m_sel->select(it, it_end);

	SmartPtr<AssemblingTuner<TAlgebra> > assTuner = m_spDomDisc->ass_tuner();
	assTuner->set_selector(m_sel.get());
}


#if 0
template <typename TDomain, typename TAlgebra>
CableAssTuner<TDomain, TAlgebra>::
CableAssTuner
(
	SmartPtr<DomainDiscretization<TDomain, TAlgebra> > domDisc,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: m_spDomDisc(domDisc), m_spApprox(approx),
  m_bRemoveGhostsFromAssemblingIterator(false), m_sel(SPNULL)
{}


template <typename TDomain, typename TAlgebra>
template <size_t dim>
void CableAssTuner<TDomain, TAlgebra>::
remove_ghosts_from_assembling_iterator()
{
	m_bRemoveGhostsFromAssemblingIterator = true;
	update_elem_lists(1 << EDGE);
}


template <typename TDomain, typename TAlgebra>
const std::vector<Vertex*>& CableAssTuner<TDomain, TAlgebra>::
vertex_vector() const
{
	return m_vVrt;
}


template <typename TDomain, typename TAlgebra>
const std::vector<Edge*>& CableAssTuner<TDomain, TAlgebra>::
edge_vector() const
{
	return m_vEdge;
}


template <typename TDomain, typename TAlgebra>
void CableAssTuner<TDomain, TAlgebra>::
update_elem_lists(uint elem_types)
{
	if ((1 << VERTEX) & elem_types)
		update_elem_list_for_elem<VERTEX>();
	if ((1 << EDGE) & elem_types)
		update_elem_list_for_elem<EDGE>();
}


template <typename TDomain, typename TAlgebra>
template <uint dim>
void CableAssTuner<TDomain, TAlgebra>::
update_elem_list_for_elem()
{
	typedef typename domain_traits<dim>::grid_base_object elem_type;
	typedef typename DoFDistribution::traits<elem_type>::const_iterator it_type;

	// create selector to store assemble elements (if requested)
	if (m_bRemoveGhostsFromAssemblingIterator && dim == EDGE)
	{
		SmartPtr<MultiGrid> mg = m_spApprox->domain()->grid();
		m_sel = make_sp(new Selector(1 << dim));
		m_sel->assign_grid(mg.get());
		m_sel->clear();
	}

	// select elements
	GridLevel gl(GridLevel::TOP, GridLevel::SURFACE, false);
	ConstSmartPtr<DoFDistribution> dd = m_spApprox->dof_distribution(gl, false);
	it_type it = dd->template begin<elem_type>(SurfaceView::MG_ALL);
	it_type it_end = dd->template end<elem_type>(SurfaceView::MG_ALL);
	for (; it != it_end; ++it)
		StoreElem<elem_type>(*it, this);

	if (m_bRemoveGhostsFromAssemblingIterator && dim == EDGE)
	{
		m_sel->select(m_vEdge.begin(), m_vEdge.end());
		SmartPtr<AssemblingTuner<TAlgebra> > assTuner = m_spDomDisc->ass_tuner();
		assTuner->set_selector(m_sel.get());
	}
}
#endif


} // namespace cable
} // namespace ug


