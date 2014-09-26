/*
 * Created on: 13.06.2014
 * 		Author: Pascal Gottmann
 * ElemDiscHH_base.h
 *
 * Based on
 * convection_diffusion.h
 * from andreasvogel
 */

#define _USE_MATH_DEFINES
#include "ElemDiscHH_fv1.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"
#include <math.h>

namespace ug{


////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ElemDiscHH_FV1<TDomain>::
ElemDiscHH_FV1(const char* functions, const char* subsets)
 : ElemDiscHH_Base<TDomain>(functions,subsets),
   m_bNonRegularGrid(false)
{
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
void ElemDiscHH_FV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check number
	if (vLfeID.size() != 4)
		UG_THROW("ElemDiscHH_FV1: Wrong number of functions given. Need exactly "<< 4);

	if (vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("ElemDiscHH FV Scheme only implemented for 1st order.");

	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
bool ElemDiscHH_FV1<TDomain>::
use_hanging() const
{
	// As this is basically a 1D discretization,
	// there will never be any hanging nodes.
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	// set local positions
	static const int refDim = TElem::dim;
	TFVGeom& geo = GeomProvider<TFVGeom>::get();

	const MathVector<refDim>* vSCVip = geo.scv_local_ips();
	const size_t numSCVip = geo.num_scv_ips();

	m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

	try{geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("KabelDiffFV1::prep_elem: Cannot update Finite Volume Geometry.");

	// set global positions
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();

	m_imSource.				set_global_ips(vSCVip, numSCVip);
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	number element_length = 0.0;
	number pre_resistance = 0.0;

// channel kinetics
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// add length of scv to element length
		element_length += scv.volume();
		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*DIAM_CONST*DIAM_CONST);


		// gating param h
		number AlphaHh = 0.07*exp(-(u(_VM_,co)+65)/20.0);
		number BetaHh = 1.0/(exp(3.0-0.1*(u(_VM_,co)+65.0))+1.0);

		// gating param m
		number AlphaHm;
		number AlphaHm_test = exp(2.5-0.1*(u(_VM_,co)+65.0))-1.0;
		if (fabs(AlphaHm_test) > 1e-5)
			AlphaHm = (2.5 - 0.1*(u(_VM_,co)+65.0)) / AlphaHm_test;
		else
			AlphaHm = 1.0;

		number BetaHm = 4.0*exp(-(u(_VM_,co)+65)/18.0);

		// gating param n
		number AlphaHn;
		number AlphaHn_test;
		AlphaHn_test = exp(1.0-0.1*(u(_VM_,co)+65.0))-1.0;
		if (fabs(AlphaHn_test) > 1e-5)
			AlphaHn = (0.1-0.01*(u(_VM_,co)+65.0)) / AlphaHn_test;
		else
			AlphaHn = 0.1;

		number BetaHn = 0.125*exp((u(_VM_,co)+65.0)/80.0);


		number rate_h = -((AlphaHh * (1.0-u(_h_,co))) - BetaHh * u(_h_,co));
		number rate_m = -((AlphaHm * (1.0-u(_m_,co))) - BetaHm * u(_m_,co));
		number rate_n = -((AlphaHn * (1.0-u(_n_,co))) - BetaHn * u(_n_,co));

		// single channel type fluxes
		const number potassium_part_of_flux = 36.0 * pow(u(_n_,co),4) * (u(_VM_,co) + 77.0);
		const number sodium_part_of_flux =  120.0 * pow(u(_m_,co),3) * u(_h_,co) * (u(_VM_, co) - 50.0);
		const number leakage_part_of_flux = 0.3 * (u(_VM_,co) + 54.4);

		// injection flux
		number inject = 0.0;

		// time for every flux the same
		number time = this->time();

		// For Different Dimensions we need different inputs
		if (dim == 3) {
			number x = vCornerCoords[0][0];
			number y = vCornerCoords[0][1];
			number z = vCornerCoords[0][2];
			(*m_Injection)(inject, 4, time, x, y, z);
		}

		if (dim == 2) {
			number x = vCornerCoords[0][0];
			number y = vCornerCoords[0][1];
			(*m_Injection)(inject, 3, time, x, y);
		}

		if (dim == 1) {
			number x = vCornerCoords[0][0];
			(*m_Injection)(inject, 2, time, x);
		}

		const number flux =   (potassium_part_of_flux
							+ sodium_part_of_flux
							+ leakage_part_of_flux);

		//std::cout<< "Injection: " << inject << std::endl;
		// fehler in defekt normal -inject/radius
		d(_VM_, co) += scv.volume()*PI*DIAM_CONST*(flux-(inject));

		// bei - defekt gates change in false direction
		d(_h_, co) += rate_h;
		d(_m_, co) += rate_m;
		d(_n_, co) += rate_n;

		/*
		std::cout << "flux: " << flux << std::endl;
		std::cout << "defekt VM: " << d(_VM_, co) << "bei Loesung: " << u(_VM_,co) << std::endl;
		std::cout << "defekt h: " << d(_h_, co) << "bei Loesung h: " << u(_h_,co)<< std::endl;
		std::cout << "defekt m: " << d(_m_, co) << "bei Loesung m: " << u(_m_,co)<< std::endl;
		std::cout << "defekt n: " << d(_n_, co) << "bei Loesung n: " << u(_n_,co)<< std::endl;
		std::cout << "A/B-Hn: " << AlphaHn <<" " << BetaHn<< std::endl;
		std::cout << "A/B-Hm: " << AlphaHm <<" " << BetaHm<< std::endl;
		std::cout << "A/B-Hh: " << AlphaHh <<" " << BetaHh<< std::endl;
		*/
	}

// diffusion part for cable equation
	// TODO: We set constant values for resistance here; this will have to be changed
	number spec_resistance = 1.0;

	MathVector<dim> grad_c;

	for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		// compute gradient and shape at ip
		VecSet(grad_c, 0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad_c, u(_VM_,sh), scvf.global_grad(sh));

		// scalar product with normal
		number diff_flux = VecDot(grad_c, scvf.normal());

		// scale by 1/resistance and by length of element
		diff_flux *= element_length / (spec_resistance*pre_resistance);

		// add to local defect
		d(_VM_, scvf.from()) -= diff_flux;
		d(_VM_, scvf.to()  ) += diff_flux;
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// TODO: This constant will have to be set elsewhere
	number spec_capacity = 1.0;

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// gating parameters time derivative
		d(_h_, co) += u(_h_, co);
		d(_m_, co) += u(_m_, co);
		d(_n_, co) += u(_n_, co);

		// potential equation time derivative
		d(_VM_, co) += PI*DIAM_CONST*scv.volume()*u(_VM_, co)*spec_capacity;
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
/*
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// loop Sub Control Volumes (SCV)
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv( ip );

		// get associated node
		const int co = scv.node_id();

// TODO: Implementiere den auskommentierten Bereich so, dass er auf unser Problem passt!
// Maybe implement injection current here instead of defect. But not necessarily.
		// Add to local rhs
		d(_VM_, co) += m_imSource[ip] * scv.volume();
	}
*/
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	number element_length = 0.0;
	number pre_resistance = 0.0;

// channel kinetics derivatives
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// add length of scv to element length
		element_length += scv.volume();
		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*DIAM_CONST*DIAM_CONST);

	// calculate several help variables for efficient calculation of derivatives
		// gating param h
		number AlphaHh = 0.07*exp(-(u(_VM_,co)+65)/20.0);
		number BetaHh = 1.0/(exp(3.0-0.1*(u(_VM_,co)+65.0))+1.0);

		// gating param m
		number AlphaHm;
		number AlphaHm_test = exp(2.5-0.1*(u(_VM_,co)+65.0))-1.0;
		if (fabs(AlphaHm_test) > 1e-5)
			AlphaHm = (2.5 - 0.1*(u(_VM_,co)+65.0)) / AlphaHm_test;
		else
			AlphaHm = 1.0;

		number BetaHm = 4.0*exp(-(u(_VM_,co)+65)/18.0);

		// gating param n
		number AlphaHn;
		number AlphaHn_test;
		AlphaHn_test = exp(1.0-0.1*(u(_VM_,co)+65.0))-1.0;
		if (fabs(AlphaHn_test) > 1e-5)
			AlphaHn = (0.1-0.01*(u(_VM_,co)+65.0)) / AlphaHn_test;
		else
			AlphaHn = 0.1;

		number BetaHn = 0.125*exp((u(_VM_,co)+65.0)/80.0);

		// gating param h derivatives
		number dAlphaHh_dVm = 0.07/20.0 * exp(-(u(_VM_,co)+65.0)/20.0);
		number help = exp(3.0-0.1*(u(_VM_,co)+65.0));
		number dBetaHh_dVm = 0.1*help / pow(help+1.0, 2);

		// gating param m derivatives
		help = 2.5 - 0.1*(u(_VM_,co)+65.0);
		number dAlphaHm_dVm;
		if (fabs(exp(help)-1.0) > 1e-5)
			dAlphaHm_dVm = -0.1 * (((1.0-help)*exp(help)-1.0)) / pow(exp(help)-1.0, 2);
		else
			dAlphaHm_dVm = 1.0;

		number dBetaHm_dVm = -4.0/18.0 * exp(-(u(_VM_,co)+65.0)/18.0);

		// gating param n derivatives
		help = 1.0 - 0.1*(u(_VM_,co)+65.0);
		number dAlphaHn_dVm;
		if (fabs(exp(help)-1.0) > 1e-5)
			dAlphaHn_dVm = -0.01 * (((1.0-help)*exp(help)-1.0)) / pow(exp(help)-1.0, 2);
		else
			dAlphaHn_dVm = 0.005;

		number dBetaHn_dVm = 0.125/80.0 * exp((u(_VM_,co)+65.0)/80.0);

	// add to Jacobian
		// derivatives of channel states
		J(_h_, co, _h_, co) += AlphaHh + BetaHh;
		J(_m_, co, _m_, co) += AlphaHm + BetaHm;
		J(_n_, co, _n_, co) += AlphaHn + BetaHn;
		J(_h_, co, _VM_, co) += -((dAlphaHh_dVm * (1.0-u(_h_,co))) - dBetaHh_dVm * u(_h_,co));
		J(_m_, co, _VM_, co) += -((dAlphaHm_dVm * (1.0-u(_m_,co))) - dBetaHm_dVm * u(_m_,co));
		J(_n_, co, _VM_, co) += -((dAlphaHn_dVm * (1.0-u(_n_,co))) - dBetaHn_dVm * u(_n_,co));

		// derivatives of potential from HH channels
		J(_VM_, co, _h_, co) += scv.volume()*PI*DIAM_CONST * 120.0*pow(u(_m_,co),3) * (u(_VM_, co) - 50.0);
		J(_VM_, co, _m_, co) += scv.volume()*PI*DIAM_CONST * 3.0*120.0*pow(u(_m_,co),2) * (u(_VM_, co) - 50.0);
		J(_VM_, co, _n_, co) += scv.volume()*PI*DIAM_CONST * 4.0*36.0*pow(u(_n_,co),3) * (u(_VM_,co) + 77.0);
		J(_VM_, co, _VM_, co) += scv.volume()*PI*DIAM_CONST * (36.0*pow(u(_n_,co),4) + 120.0*pow(u(_m_,co),3)*u(_h_,co) + 0.3);
	}

// "diffusion" derivatives
	number spec_resistance = 1.0;

	// loop Sub Control Volume Faces (SCVF)
	for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		// loop shape functions
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			// scalar product with normal
			number d_diff_flux = element_length * VecDot(scvf.global_grad(sh), scvf.normal());

			// scale by 1/resistance and by length of element
			d_diff_flux *= element_length / (spec_resistance*pre_resistance);

			// add flux term to local matrix
			J(_VM_, scvf.from(), _VM_, sh) -= d_diff_flux;
			J(_VM_, scvf.to()  , _VM_, sh) += d_diff_flux;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// TODO: This constant will have to be set elsewhere
	number spec_capacity = 1.0;

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// gating parameters
		J(_h_, co, _h_, co) += 1.0;
		J(_m_, co, _m_, co) += 1.0;
		J(_n_, co, _n_, co) += 1.0;

		// potential equation
		J(_VM_, co, _VM_, co) += PI*DIAM_CONST*scv.volume()*spec_capacity;
	}
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
void ElemDiscHH_FV1<TDomain>::
register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ElemDiscHH_FV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ElemDiscHH_FV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ElemDiscHH_FV1<Domain3d>;
#endif

} // namespace ug

