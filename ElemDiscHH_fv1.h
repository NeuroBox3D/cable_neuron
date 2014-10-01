/*
 * Created on: 13.06.2014
 * 		Author: Pascal Gottmann
 * ElemDiscHH_base.h
 *
 * Based on
 * convection_diffusion.h
 * from andreasvogel
 */

#ifndef __H__ElemDisc_fv1__
#define __H__ElemDisc_fv1__

// library intern headers
#include "ElemDiscHH_base.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"

namespace ug{



template<	typename TDomain>
class ElemDiscHH_FV1 : public ElemDiscHH_Base<TDomain>
{
	private:
	///	Base class type
		typedef ElemDiscHH_Base<TDomain> base_type;

	///	Own type
		typedef ElemDiscHH_FV1<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor
		ElemDiscHH_FV1(SmartPtr<ApproximationSpace<TDomain> > approx,
				const char* functions, const char* subsets);

		// set diameter for dendrit
		void set_diameter(const number d);

		// Problem is solved with IFunction position could get with vcornercoords
		void set_injection(IFunction<number>& functor) { m_Injection = &functor;}


	private:
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	finishes the loop over all elements
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the stiffness part of the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the mass part of the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local right hand side
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);


	private:
	///	abbreviation for the local solution
		static const size_t _VM_ = 0;
		static const size_t _h_ = 1;
		static const size_t _m_ = 2;
		static const size_t _n_ = 3;

		IFunction<number>* m_Injection;

		using base_type::m_spec_capa;
		using base_type::m_imSource;

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns if hanging nodes are needed
		virtual bool use_hanging() const;


	protected:
	///	current regular grid flag
		bool m_bNonRegularGrid;

	/// approximation space (needed before given by domainDisc)
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	/// dendritic radius attachment and accessor
		ANumber m_aDiameter;
		Grid::AttachmentAccessor<Edge, ANumber> m_aaDiameter;

	///	register utils
	///	\{
		void register_all_funcs(bool bHang);
		template <typename TElem, typename TFVGeom> void register_func();
	/// \}
};



} // end namespace ug


#endif /*__H__ElemDisc_fv1__*/
