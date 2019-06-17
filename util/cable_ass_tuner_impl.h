/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2015-08-08
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

#include "cable_ass_tuner.h"
#include "lib_disc/domain_traits.h"
#include "lib_grid/tools/grid_level.h"
#include "lib_grid/grid/grid_base_objects.h"


namespace ug {
namespace cable_neuron {


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


} // namespace cable_neuron
} // namespace ug


