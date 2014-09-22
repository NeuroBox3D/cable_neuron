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
   m_spConvShape(new ConvectionShapesNoUpwind<dim>),
   m_bNonRegularGrid(false)
{
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
void ElemDiscHH_FV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
//	check number
	if(vLfeID.size() != 4)
		UG_THROW("ElemDiscHH_FV1: Wrong number of functions given. "
				"Need exactly "<<4);

	if(vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("ElemDiscHH FV Scheme only implemented for 1st order.");

//	remember
	m_bNonRegularGrid = bNonRegularGrid;

//	update assemble functions
	register_all_funcs(m_bNonRegularGrid);
	std::cout << "preparesetting is working"<< std::endl;
}

template<typename TDomain>
bool ElemDiscHH_FV1<TDomain>::
use_hanging() const
{
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

//	check, that upwind has been set
	if(m_spConvShape.invalid())
		UG_THROW("ElemDiscHH_FV1::prep_elem_loop:"
						" Upwind has not been set.");

//	set local positions
	if(!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_imDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);

		//	init upwind for element type
		if(!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
			UG_THROW("ElemDiscHH_fv1::prep_elem_loop:"
						" Cannot init upwind for element type.");
	}
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

// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}UG_CATCH_THROW("ElemDiscHH::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	set local positions
	if(TFVGeom::usesHangingNodes)
	{
		const int refDim = TElem::dim;
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_imDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip);

		if(m_spConvShape.valid())
			if(!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
				UG_THROW("ElemDiscHH_FV1::prep_elem_loop:"
								" Cannot init upwind for element type.");
	}

	//	set global positions
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_imDiffusion.			set_global_ips(vSCVFip, numSCVFip);
	m_imMassScale.			set_global_ips(vSCVip, numSCVip);
	m_imSource.				set_global_ips(vSCVip, numSCVip);
}

template <class TVector>
static TVector CalculateCenter(GridObject* o, const TVector* coords)
{
	TVector v;
	VecSet(v, 0);

	size_t numCoords = 0;
	switch(o->base_object_id()){
		case VERTEX: numCoords = 1; break;
		case EDGE: numCoords = static_cast<Edge*>(o)->num_vertices(); break;
		case FACE: numCoords = static_cast<Face*>(o)->num_vertices(); break;
		case VOLUME: numCoords = static_cast<Volume*>(o)->num_vertices(); break;
		default: UG_THROW("Unknown element type."); break;
	}

	for(size_t i = 0; i < numCoords; ++i)
		VecAdd(v, v, coords[i]);

	if(numCoords > 0)
		VecScale(v, v, 1. / (number)numCoords);

	return v;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	//	Diff. Tensor times Gradient
		MathVector<dim> Dgrad, grad_c, Dgrad_c;

		double element_length = 0.0;

	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		const int co = scv.node_id();
		element_length += scv.volume();
		// 	getting local values
		const LocalVector& sol = u;
		std::cout<< "LocalVec in Jac_A working" << std::endl;

		const LocalVectorTimeSeries& vLocSol = *this->local_time_solutions();
	//	remember local solutions
		number dt = vLocSol.time(0) - vLocSol.time(1);
		//number time = vLocSol.time(0);

			J(_VM_, co, _VM_, co) += scv.volume() * ( 36*pow(sol(_n_,co),4)
								+ 120*pow(sol(_m_,co),3)*sol(_h_,co) + 0.3);
			/*double AlphaHn;
			double AlphaHn_test;
			AlphaHn_test = (exp(1-0.1*(sol(_VM_,co)))-1);
			if (fabs(AlphaHn_test) > 0.00001)
			{
				AlphaHn = (0.1-0.01*sol(_VM_,co))/(exp(1-0.1*(sol(_VM_,co)))-1);
			} else {
				AlphaHn = 0.1;
			}
;
			double BetaHn = 0.125*exp(sol(_VM_,co)/80);
			double AlphaHm = (2.5 - 0.1*sol(_VM_,co))/(exp(2.5-0.1*sol(_VM_,co))-1);
			double BetaHm = 4*exp(-1*sol(_VM_,co)/18);
			double AlphaHh = 0.07*exp(-1*sol(_VM_,co)/20);
			double BetaHh = 1/(exp(3-0.1*sol(_VM_,co))+1);*/


			double AlphaHn;
			double AlphaHn_test;
			AlphaHn_test = (exp(1-0.1*(sol(_VM_,co)+65))-1);
			if (fabs(AlphaHn_test) > 0.00001)
			{
				AlphaHn = (0.1-0.01*(sol(_VM_,co)+65))/(exp(1-0.1*(sol(_VM_,co)+65))-1);
			} else {
				AlphaHn = 0.1;
			}
			double BetaHn = 0.125*exp((sol(_VM_,co)+65)/80);
			double AlphaHm = (2.5 - 0.1*(sol(_VM_,co)+65))/(exp(2.5-0.1*(sol(_VM_,co)+65))-1);
			double BetaHm = 4*exp(-1*(sol(_VM_,co)+65)/18);
			double AlphaHh = 0.07*exp(-1*(sol(_VM_,co)+65)/20);
			double BetaHh = 1/(exp(3-0.1*(sol(_VM_,co)+65))+1);



			/*long double AlphaHn = ((0.01*(sol(_VM_,co) + 55))/(1-exp(-0.1*(sol(_VM_,co)+55))));
			long double BetaHn  = ((0.125*exp(-0.0125*(sol(_VM_, co)+65))));
			long double AlphaHm = ((0.1*(sol(_VM_,co)+40))/(1-exp(-0.1*(sol(_VM_,co)+40))));
			long double BetaHm  = (4*exp(-0.0556*(sol(_VM_,co)+65)));
			long double AlphaHh = (0.07*exp(-0.05*(sol(_VM_,co)+65)));
			long double BetaHh  = (1/(1 + exp(-0.1*(sol(_VM_,co)+35))));*/
			//Ableitungen nach VM
			/*double AlphaHn1 = (0.01*(1-exp(-0.1*sol(_VM_,co) -5.5)) - (0.01*sol(_VM_,co) +0.55)*(0.1*exp(0.1*sol(_VM_,co)-5.5)))/(pow((1-exp(-0.1*sol(_VM_,co)-5.5)),2));
			double BetaHn1 = 0.125 * -0.0125 * (exp(0.0125*sol(_VM_,co) - 0.8125));
			double AlphaHm1 = ((0.1*(1-exp(-0.1*sol(_VM_, co) + 4))) - ((0.1*sol(_VM_, co) + 4)*(0.1*exp(-0.1*sol(_VM_, co) + 4 ))))/pow((1-exp(-0.1*sol(_VM_, co) + 4)),2);
			double BetaHm1 = 4 * -0.556*exp(-0.556*sol(_VM_, co) + 3.614);
			double AlphaHh1 = -0.05*0.07 *exp(0.05*sol(_VM_, co) + 3.25);
			double BetaHh1 = -1*(1+exp(-0.1*sol(_VM_,co) -3.5))/pow((1+exp(-0.1*sol(_VM_, co) -3.5)),2);*/

			// mit + 65
			double AlphaHn2 = (-0.01/(exp(1-0.1*(sol(_VM_,co) + 65))-1)
					- (0.1 - 0.01*(sol(_VM_,co) + 65))*(-0.1*exp(1-0.1*(sol(_VM_,co) + 65)))/((exp(1-0.1*(sol(_VM_,co) + 65))-1)*(exp(1-0.1*(sol(_VM_,co) + 65))-1))
					- (-0.01/(exp(1-0.1*(sol(_VM_,co) + 65))-1) - (0.1 - 0.01*(sol(_VM_,co) + 65))*(-0.1*exp(1-0.1*(sol(_VM_,co) + 65)))/((exp(1-0.1*(sol(_VM_,co) + 65))-1)*(exp(1-0.1*(sol(_VM_,co) + 65))-1))));
			double BetaHn2 = 0.125/80.0*exp(-(sol(_VM_,co) + 65)/80);

			double AlphaHm2 = (-0.1/(exp(2.5 - 0.1*(sol(_VM_,co) + 65)) -1) + exp((sol(_VM_,co) + 65)/10 - 2.5)*(2.5 - 0.1*(sol(_VM_,co) +65))*0.1/((exp(2.5 - 0.1*(sol(_VM_,co) + 65)) -1)*(exp(2.5 - 0.1*(sol(_VM_,co) + 65)) -1)) -
					(-0.1/(exp(2.5 - 0.1*(sol(_VM_,co) + 65)) -1) + exp((sol(_VM_,co) + 65)/10 - 2.5)*(2.5 - 0.1*(sol(_VM_,co) +65))*0.1/((exp(2.5 - 0.1*(sol(_VM_,co) + 65)) -1)*(exp(2.5 - 0.1*(sol(_VM_,co) + 65)) -1))));
			double BetaHm2 = -4.0/18*exp(-1*(sol(_VM_,co)+65)/18);

			double AlphaHh2 = -0.07/20*exp(-1*(sol(_VM_,co) + 65)/20);

			double BetaHh2 = (-0.07/20*exp(-1*(sol(_VM_,co) + 65)/20) +
					 0.1*exp(3 - 0.1*(sol(_VM_,co) + 65))/((exp(3 - 0.1*(sol(_VM_,co) + 65)) + 1)*(exp(3 - 0.1*(sol(_VM_,co) + 65)) + 1)) );



			//  ohne + 65
			/*double AlphaHn2 = (-0.01/(exp(1-0.1*(sol(_VM_,co)))-1)
					- (0.1 - 0.01*(sol(_VM_,co)))*(-0.1*exp(1-0.1*(sol(_VM_,co))))/((exp(1-0.1*(sol(_VM_,co)))-1)*(exp(1-0.1*(sol(_VM_,co)))-1))
					- (-0.01/(exp(1-0.1*(sol(_VM_,co)))-1) - (0.1 - 0.01*(sol(_VM_,co)))*(-0.1*exp(1-0.1*(sol(_VM_,co))))/((exp(1-0.1*(sol(_VM_,co)))-1)*(exp(1-0.1*(sol(_VM_,co)))-1))));
			double BetaHn2 = 0.125/80.0*exp(-(sol(_VM_,co))/80);

			double AlphaHm2 = (-0.1/(exp(2.5 - 0.1*(sol(_VM_,co))) -1) + exp((sol(_VM_,co))/10 - 2.5)*(2.5 - 0.1*(sol(_VM_,co)))*0.1/((exp(2.5 - 0.1*(sol(_VM_,co))) -1)*(exp(2.5 - 0.1*(sol(_VM_,co))) -1)) -
					(-0.1/(exp(2.5 - 0.1*(sol(_VM_,co))) -1) + exp((sol(_VM_,co))/10 - 2.5)*(2.5 - 0.1*(sol(_VM_,co)))*0.1/((exp(2.5 - 0.1*(sol(_VM_,co))) -1)*(exp(2.5 - 0.1*(sol(_VM_,co))) -1))));
			double BetaHm2 = -4.0/18*exp(-1*(sol(_VM_,co))/18);

			double AlphaHh2 = -0.07/20*exp(-1*(sol(_VM_,co))/20);

			double BetaHh2 = (-0.07/20*exp(-1*(sol(_VM_,co))/20) +
					 0.1*exp(3 - 0.1*(sol(_VM_,co)))/((exp(3 - 0.1*(sol(_VM_,co))) + 1)*(exp(3 - 0.1*(sol(_VM_,co))) + 1)) );*/
			/*double AlphaHn3 = (0.01*(1-exp(-0.1*sol(_VM_,co) -5.5)) - (0.01*sol(_VM_,co) +0.55)*(0.1*exp(0.1*sol(_VM_,co)-5.5)))/(pow((1-exp(-0.1*sol(_VM_,co)-5.5)),2));
			double BetaHn3 = 0.125 * -0.0125 * (exp(0.0125*sol(_VM_,co) - 0.8125));
			double AlphaHm3 = ((0.1*(1-exp(-0.1*sol(_VM_, co) + 4))) - ((0.1*sol(_VM_, co) + 4)*(0.1*exp(-0.1*sol(_VM_, co) + 4 ))))/pow((1-exp(-0.1*sol(_VM_, co) + 4)),2);
			double BetaHm3 = 4 * -0.556*exp(-0.556*sol(_VM_, co) + 3.614);
			double AlphaHh3 = -0.05*0.07 *exp(0.05*sol(_VM_, co) + 3.25);
			double BetaHh3 = -1*(1+exp(-0.1*sol(_VM_,co) -3.5))/pow((1+exp(-0.1*sol(_VM_, co) -3.5)),2);*/

			// here we need f'
			J(_h_, co, _h_, co) += scv.volume() * (AlphaHh * (1-sol(_h_,co)) + BetaHh*(sol(_h_,co)));
			J(_m_, co, _m_, co) += scv.volume() * (AlphaHm * (1-sol(_m_,co)) + BetaHm*(sol(_m_,co)));
			J(_n_, co, _n_, co) += scv.volume() * (AlphaHn * (1-sol(_n_,co)) + BetaHn*(sol(_n_,co)));
			// kapa von 1 Wert nch volume
			//J(_VM_, co, _VM_, co) += scv.volume() * (1 * 0.36*pow(sol(_n_,co),4) + 1.20*pow(sol(_m_,co),3)*sol(_h_,co) + 0.003);

			/*J(_n_, co, _VM_, co) += scv.volume() *((AlphaHn1 - BetaHn1)*(sol(_n_, co)));
			J(_m_, co, _VM_, co) += scv.volume() *((AlphaHm1 - BetaHm1)*(sol(_m_, co)));
			J(_h_, co, _VM_, co) += scv.volume() *((AlphaHh1 - BetaHh1)*(sol(_h_, co)));*/
			J(_n_, co, _VM_, co) += scv.volume() *((AlphaHn2 - BetaHn2)*(sol(_n_, co)));
			J(_m_, co, _VM_, co) += scv.volume() *((AlphaHm2 - BetaHm2)*(sol(_m_, co)));
			J(_h_, co, _VM_, co) += scv.volume() *((AlphaHh2 - BetaHh2)*(sol(_h_, co)));
			/*double m_m = sol(_m_,co) - 1/(AlphaHm+BetaHm);
			double m_h = sol(_h_,co) - 1/(AlphaHh+BetaHh);
			double m_n = sol(_n_,co) - 1/(AlphaHn+BetaHn);*/
			J(_VM_, co, _n_, co) += scv.volume() * 36*4*pow(sol(_n_,co),3)*(sol(_VM_,co) + 77 );

			J(_VM_, co, _m_, co) += scv.volume() * 120*3*pow(sol(_m_,co),2)*sol(_h_,co)*(sol(_VM_, co) - 50);

			J(_VM_, co, _h_, co) += scv.volume() * ( 120*pow(sol(_m_,co),3)*(sol(_VM_, co)- 50));

			// here we need f' ableitungen nach einzelnen kompos
			/*J(_h_, co, _h_, co) += scv.volume() * m_h;
			J(_m_, co, _m_, co) += scv.volume() * m_m;
			J(_n_, co, _n_, co) += scv.volume() * m_n;*/
			// Ableigung nach h,m,n ohne Ableitung wird matrix singulï¿½r
			/*J(_h_, co, _h_, co) += scv.volume() * ((-AlphaHh - BetaHh));
			J(_m_, co, _m_, co) += scv.volume() * ((-AlphaHm - BetaHm));
			J(_n_, co, _n_, co) += scv.volume() * ((-AlphaHn - BetaHn));*/

			// Ableitungen nach VM von h,m,n
			/*J(_h_, co, _VM_, co) += scv.volume() * ;
			J(_m_, co, _VM_, co) += scv.volume() * ;
			J(_n_, co, _VM_, co) += scv.volume() * ;*/

			/*J(_n_, co, _n_, co) += scv.volume() * 0.36*4*pow(sol(_n_,co),3)*(sol(_VM_,co) + 77 );

			J(_m_, co, _m_, co) += scv.volume() * 1.20*3*pow(sol(_m_,co),2)*sol(_h_,co)*(sol(_VM_, co) - 50);

			J(_h_, co, _h_, co) += scv.volume() * ( 1.20*pow(sol(_m_,co),3)*(sol(_VM_, co)- 50));*/

			//std::cout << "volume: " << scv.volume() << std::endl;
	}


	// erstmal hils var fuer kapa und diam
	double spec_resistance = 1.0;
	double pre_resistance = 0;



	// fluesse hier rein
	std::cout << "add_jac_A_elem" << std::endl;



	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		// here we need diam and kapa questioning
		////////////////////////////////////////////////////
		// Diffusive Term
		////////////////////////////////////////////////////
			const typename TFVGeom::SCV& scv = geo.scv(ip);
			std::cout << "In Diffusion add jac" << std::endl;
				// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					// bei ableitung muesste ja aus u VM 1 werden
					pre_resistance += scv.volume() / (0.25*PI*DIAM_CONST*DIAM_CONST);
					VecScaleAppend(grad_c, u(_VM_,sh), scvf.global_grad(sh));

					//	scale by length of element
					VecScaleAppend(Dgrad_c, element_length, grad_c);
					double D_diff_flux = VecDot(Dgrad_c, scvf.normal());

					// 	scale by 1/resistance
					D_diff_flux /= spec_resistance*pre_resistance;

				// 	Add flux term to local matrix // HIER MATRIXINDIZES!!!
					UG_ASSERT((scvf.from() < J.num_row_dof(_VM_)) && (scvf.to() < J.num_col_dof(_VM_)),
							  "Bad local dof-index on element with object-id " << elem->base_object_id()
							  << " with center: " << CalculateCenter(elem, vCornerCoords));

					//std::cout << "D_diff_flux: " << D_diff_flux << std::endl;
					J(_VM_, scvf.from(), _VM_, sh) -= D_diff_flux;
					J(_VM_, scvf.to()  , _VM_, sh) += D_diff_flux;
			}
		}

}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// fluesse hier rein
// 	get finite volume geometry
	std::cout << "add_jac_M_elem" << std::endl;
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	//if(!m_imMassScale.data_given()) return;

	double spec_capacity = 1.0;

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

		// gating parameters
		J(_h_, co, _h_, co) += 1.0;
		J(_m_, co, _m_, co) += 1.0;
		J(_n_, co, _n_, co) += 1.0;

		// potential equation
		J(_VM_, co, _VM_, co) += PI*DIAM_CONST*scv.volume()*spec_capacity;
	}

//	m_imMass part does not explicitly depend on associated unknown function
	std::cout << "add_jac_M_elem  working" << std::endl;
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	std::cout << "add_def_A_elem" << std::endl;


// zum durchgehen der elemente
	std::vector<DoFIndex> multInd;
	//std::vector<position_type> vPos;

	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	double element_length = 0.0;
	double pre_resistance = 0.0;

	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);
		// 	get associated node
		const int co = scv.node_id();
		// 	Add to local defect

		// add length of scv to element length
		element_length += scv.volume();

		// add "pre_resistance" parts
		pre_resistance += scv.volume() / (0.25*PI*DIAM_CONST*DIAM_CONST);

		std::cout << " fehler kommt gleich" << std::endl;
		const LocalVectorTimeSeries& vLocSol = *this->local_time_solutions();
		const LocalVector& sol = vLocSol.solution(0);
		const LocalVector& solnew = vLocSol.solution(1);

	//	remember local solutions
		number dt = vLocSol.time(0) - vLocSol.time(1);
		number time = vLocSol.time(0);

		// umwandlung der positionen muss fï¿½r 3d noch geï¿½ndert werden
		std::ostringstream stream_x;
		stream_x << vCornerCoords[0];
		std::string String_x = stream_x.str();
		/*std::cout << "tests" << std::endl;
		std::cout << String_x << std::endl;*/
		size_t end = String_x.find(")",0);
		size_t beg = 1;
		std::string String_xx = String_x.substr(beg, end-1);
		//std::cout << String_xx << std::endl;
		// hier kommt fehler
		std::stringstream String_xxx;
		number x;
		String_xxx << String_xx;
		String_xxx >> x;
		//std::cout << String_xxx << std::endl;
		//std::cout << x << std::endl;

		//const LFEID& test = u.local_finite_element_id(0);
		//const LocalVector& teser = u.lo
		//InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, LFEID);
		// hier punkt und sooo...
		std::cout << "teser: " << vCornerCoords[0] << " umgewandelt: " << x << std::endl;
		//std::cout << "vPos: " << vPos << std::endl;

		//std::cout << "x-Wert: " << x << std::endl;

		//std::cout << "solution ist doof: " << sol(_VM_,co) << "bei Zeitschritt: " << dt << "mit zweiter Loesung: "<< solnew(_VM_,co)<<std::endl;

		/*long double AlphaHn = ((0.01*(sol(_VM_,co) + 55))/(1-exp(-0.1*(sol(_VM_,co)+55))));
		long double BetaHn  = ((0.125*exp(-0.0125*(sol(_VM_, co)+65))));
		long double AlphaHm = ((0.1*(sol(_VM_,co)+40))/(1-exp(-0.1*(sol(_VM_,co)+40))));
		long double BetaHm  = (4*exp(-0.0556*(sol(_VM_,co)+65)));
		long double AlphaHh = (0.07*exp(-0.05*(sol(_VM_,co)+65)));
		long double BetaHh  = (1/(1 + exp(-0.1*(sol(_VM_,co)+35))));*/

		/*double AlphaHn;
		double AlphaHn_test;
		AlphaHn_test = (exp(1-0.1*(sol(_VM_,co)))-1);
		if (fabs(AlphaHn_test) > 0.00001)
		{
			AlphaHn = (0.1-0.01*sol(_VM_,co))/(exp(1-0.1*(sol(_VM_,co)))-1);
		} else {
			AlphaHn = 0.1;
		}
;
		double BetaHn = 0.125*exp(sol(_VM_,co)/80);
		double AlphaHm = (2.5 - 0.1*sol(_VM_,co))/(exp(2.5-0.1*sol(_VM_,co))-1);
		double BetaHm = 4*exp(-1*sol(_VM_,co)/18);
		double AlphaHh = 0.07*exp(-1*sol(_VM_,co)/20);
		double BetaHh = 1/(exp(3-0.1*sol(_VM_,co))+1);*/


		double AlphaHn;
		double AlphaHn_test;
		AlphaHn_test = (exp(1-0.1*(sol(_VM_,co)+65))-1);
		if (fabs(AlphaHn_test) > 0.00001)
		{
			AlphaHn = (0.1-0.01*(sol(_VM_,co)+65))/(exp(1-0.1*(sol(_VM_,co)+65))-1);
		} else {
			AlphaHn = 0.1;
		}
		double BetaHn = 0.125*exp((sol(_VM_,co)+65)/80);
		double AlphaHm = (2.5 - 0.1*(sol(_VM_,co)+65))/(exp(2.5-0.1*(sol(_VM_,co)+65))-1);
		double BetaHm = 4*exp(-1*(sol(_VM_,co)+65)/18);
		double AlphaHh = 0.07*exp(-1*(sol(_VM_,co)+65)/20);
		double BetaHh = 1/(exp(3-0.1*(sol(_VM_,co)+65))+1);


		double flux_m = -((AlphaHm * (1-sol(_m_,co))) - BetaHm * sol(_m_,co));
		double flux_h = -((AlphaHh * (1-sol(_h_,co))) - BetaHh * sol(_h_,co));
		double flux_n = -((AlphaHn * (1-sol(_n_,co))) - BetaHn * sol(_n_,co));

		//std::cout << "m: " << m_m << "h: " << m_h << "n: " << m_n << std::endl;
		//std::cout << sol(_VM_,co) << std::endl;

		//d(_h_, co) += scv.volume() * AlphaHh * (1-sol(_h_,co)) - BetaHh*(sol(_h_,co));


		//std::cout << "H " << (scv.volume() * AlphaHm * (1-sol(_m_,co)) - BetaHm*(sol(_m_,co))) << std::endl;

		//d(_m_, co) += scv.volume() * AlphaHm * (1-sol(_m_,co)) - BetaHm*(sol(_m_,co));
		//d(_n_, co) += scv.volume() * AlphaHn * (1-sol(_n_,co)) - BetaHn*(sol(_n_,co));
		// capazitï¿½ï¿½t erstmal nicht benutzen in diesem fall nur im ersten wert
		//const number capacitive_part_of_flux = m_capacity * ( sol(_VM_, co) - oldSol(_VM_, co) ) / 0.01 ;

		const number potassium_part_of_flux = 36*pow(sol(_n_,co),4)*(sol(_VM_,co) + 77);
		const number sodium_part_of_flux =  120*pow(sol(_m_,co),3)*sol(_h_,co) * (sol(_VM_, co) - 50);
		const number leakage_part_of_flux = 0.3*(sol(_VM_,co) + 54.4);
		number inject = 0;
		// laesst sich spaeter mit for schleife bezug auf dim loesen
		/*std::cout << "Eig-XWert: " << vCornerCoords[0] << std::endl;
		std::cout << "Eig-yWert: " << vCornerCoords[1] << std::endl;
		std::cout << "Eig-zWert: " << vCornerCoords[2] << std::endl;*/
		std::cout << "X-Wert: " << x << std::endl;
		// this works only in 1D
		(*m_Injection)(inject, 2, time, x);


		const number flux =   (potassium_part_of_flux
							+ sodium_part_of_flux
							+ leakage_part_of_flux);

		//std::cout<< "Injection: " << inject << std::endl;
		// fehler in defekt normal -inject/radius
		d(_VM_, co) += scv.volume()*PI*DIAM_CONST*(flux-(inject)); // * scv.volume; //scv.volume() * flux * 0.001;
		//value1 = flux;
		//std::cout << "defekt: " << d(_VM_, co) << std::endl;
		// bei - deffekt lŠufts in die falsche richtung
		d(_h_, co) += flux_h;//*((sol(_h_,co)-sol(_h_,co)-flux_h));//-valueh);
		d(_m_, co) += flux_m;//flux_m;//*((sol(_m_,co)-sol(_m_,co)-flux_m));//-valuem);//0;
		d(_n_, co) += flux_n;//*((sol(_n_,co)-sol(_n_,co)-flux_n));//-valuen);
		//

		//std::cout << "volumen*flux: " << (scv.volume() *flux) << std::endl;
		/*number valueh = flux_h;
		number valuem = flux_m;
		number valuen = flux_n;*/
		//if (x == 0.5) {
			std::cout << "flux: " << flux << std::endl;
			std::cout << "defekt VM: " << d(_VM_, co) << "bei Loesung: " << sol(_VM_,co) << std::endl;
			std::cout << "defekt h: " << d(_h_, co) << "bei Loesung h: " << sol(_h_,co)<< std::endl;
			std::cout << "defekt m: " << d(_m_, co) << "bei Loesung m: " << sol(_m_,co)<< std::endl;
			std::cout << "defekt n: " << d(_n_, co) << "bei Loesung n: " << sol(_n_,co)<< std::endl;
			std::cout << "A/B-Hn: " << AlphaHn <<" " << BetaHn<< std::endl;
			std::cout << "A/B-Hm: " << AlphaHm <<" " << BetaHm<< std::endl;
			std::cout << "A/B-Hh: " << AlphaHh <<" " << BetaHh<< std::endl;
		//}
	}

	// erstmal hils var fuer spec_resistance und diam
	double spec_resistance = 1.0;
	double diam = 10;
	// set source data
	//scvf//m_imSource
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
	/////////////////////////////////////////////////////
	// Diffusive Term
	/////////////////////////////////////////////////////

	//	to compute D \nabla c
		MathVector<dim> Dgrad_c, grad_c;

	// 	compute gradient and shape at ip
		VecSet(grad_c, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad_c, u(_VM_,sh), scvf.global_grad(sh));

	//	scale by length of element
		VecScaleAppend(Dgrad_c, element_length, grad_c);
		double diff_flux = VecDot(Dgrad_c, scvf.normal());

	// 	scale by 1/resistance
		diff_flux /= spec_resistance*pre_resistance;

		std::cout << "diff_flux: " << diff_flux << std::endl;
	// 	Add to local defect
		d(_VM_, scvf.from()) -= diff_flux;
		d(_VM_, scvf.to()  ) += diff_flux;
		std::cout << "defekt diff fluss: " << diff_flux << std::endl;
	}



}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	std::cout << "add_def_M_elem" << std::endl;
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	double spec_capacity = 1.0;

	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

		// gating parameters
		d(_h_, co) += u(_h_, co);
		d(_m_, co) += u(_m_, co);
		d(_n_, co) += u(_n_, co);

		// potential equation
		d(_VM_, co) += PI*DIAM_CONST*scv.volume()*u(_VM_, co)*spec_capacity;

	}


}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// loop Sub Control Volumes (SCV)
	for (size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = geo.scv( ip );

		// get associated node
		const int co = scv.node_id();

/* TODO: Implementiere den auskommentierten Bereich so, dass er auf unser Problem passt!
		// Add to local rhs
		d(_VM_, co) += m_imSource[ip] * scv.volume();
*/
	}


}






//
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
ex_value(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         const LocalVector& u,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<TFVGeom::dim> vLocIP[],
         const size_t nip,
         bool bDeriv,
         std::vector<std::vector<number> > vvvDeriv[])
{
	std::cout << "ex_value" << std::endl;
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

//	FV1 SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				vValue[ip] += u(_VM_, sh) * scvf.shape(sh);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_VM_][sh] = scvf.shape(sh);
		}
	}
//	FV1 SCV ip
	else if(vLocIP == geo.scv_local_ips())
	{
	//	solution at ip
		for(size_t sh = 0; sh < numSH; ++sh)
			vValue[sh] = u(_VM_, sh);

	//	set derivatives if needed
		if(bDeriv)
			for(size_t sh = 0; sh < numSH; ++sh)
				for(size_t sh2 = 0; sh2 < numSH; ++sh2)
					vvvDeriv[sh][_VM_][sh2] = (sh==sh2) ? 1.0 : 0.0;
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		number vShape[numSH];

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
				vValue[ip] += u(_VM_, sh) * vShape[sh];

		//	compute derivative w.r.t. to unknowns iff needed
		//	\todo: maybe store shapes directly in vvvDeriv
			if(bDeriv)
				for(size_t sh = 0; sh < numSH; ++sh)
					vvvDeriv[ip][_VM_][sh] = vShape[sh];
		}
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
ex_grad(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFVGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
	std::cout << "ex_grad" << std::endl;
// 	Get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = ref_elem_type::dim;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

//	FV1 SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			VecSet(vValue[ip], 0.0);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(vValue[ip], u(_VM_, sh), scvf.global_grad(sh));

			if(bDeriv)
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_VM_][sh] = scvf.global_grad(sh);
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < numSH; ++sh)
				VecScaleAppend(locGrad, u(_VM_, sh), vLocGrad[sh]);

		//	compute global grad
			mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
			MatVecMult(vValue[ip], JTInv, locGrad);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < numSH; ++sh)
					MatVecMult(vvvDeriv[ip][_VM_][sh], JTInv, vLocGrad[sh]);
		}
	}
};

////////////////////////////////////////////////////////////////////////////////
//	upwind
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ElemDiscHH_FV1<TDomain>::
set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}





////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void ElemDiscHH_FV1<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		register_func<RegularEdge, HFV1Geometry<RegularEdge, dim> >();
	}
}
#endif

#ifdef UG_DIM_2
template<>
void ElemDiscHH_FV1<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else
	{
		register_func<RegularEdge, HFV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, HFV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, HFV1Geometry<Quadrilateral, dim> >();
	}
}
#endif

#ifdef UG_DIM_3
template<>
void ElemDiscHH_FV1<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else
	{
		register_func<RegularEdge, HFV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, HFV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, HFV1Geometry<Quadrilateral, dim> >();
		register_func<Tetrahedron, HFV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, HFV1Geometry<Prism, dim> >();
		register_func<Pyramid, HFV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, HFV1Geometry<Hexahedron, dim> >();
	}
}
#endif

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


//	exports
	m_exValue->	   template set_fct<T,refDim>(id, this, &T::template ex_value<TElem, TFVGeom>);
	m_exGrad->template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem, TFVGeom>);
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

