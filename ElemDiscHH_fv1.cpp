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
		m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imVectorSource.template 	set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imReactionRate.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReaction.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReactionRateExpl.template set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReactionExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imSourceExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMass.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);

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
		m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imVectorSource.template 	set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imReactionRate.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReaction.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReactionRateExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReactionExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceExpl.template		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imMass.template 			set_local_ips<refDim>(vSCVip,numSCVip);

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
	m_imVelocity.			set_global_ips(vSCVFip, numSCVFip);
	m_imFlux.				set_global_ips(vSCVFip, numSCVFip);
	m_imSource.				set_global_ips(vSCVip, numSCVip);
	m_imVectorSource.		set_global_ips(vSCVFip, numSCVFip);
	m_imReactionRate.		set_global_ips(vSCVip, numSCVip);
	m_imReactionRateExpl.	set_global_ips(vSCVip, numSCVip);
	m_imReactionExpl.		set_global_ips(vSCVip, numSCVip);
	m_imSourceExpl.			set_global_ips(vSCVip, numSCVip);
	m_imReaction.			set_global_ips(vSCVip, numSCVip);
	m_imMassScale.			set_global_ips(vSCVip, numSCVip);
	m_imMass.				set_global_ips(vSCVip, numSCVip);
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
		MathVector<dim> Dgrad;

	//	get conv shapes
		const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);


	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		const int co = scv.node_id();
		// 	getting local values
		const LocalVector& sol = u;
		std::cout<< "LocalVec in Jac_A working" << std::endl;

		// 	Add to local matrix
			//std::cout << "before zeugs dass ich mache!" << std::endl;
			//std::cout << "co: " << co << std::endl;
			//std::cout << "VM: " << sol(_VM_, co) << "h: " << sol(_h_, co) << "m: " << sol(_m_, co)<< "n: " << sol(_n_, co) << std::endl;
		// umwandlung der positionen muss für 3d noch geändert werden

			J(_VM_, co, _VM_, co) += scv.volume() * ( 0.36*pow(sol(_n_,co),4)
								+ 1.20*pow(sol(_m_,co),3)*sol(_h_,co) + 0.003);

			double AlphaHn = (0.1-0.01*sol(_VM_,co))/(exp(1-0.1*(sol(_VM_,co)))-1);
			double BetaHn = 0.125*exp(sol(_VM_,co)/80);
			double AlphaHm = (2.5 - 0.1*sol(_VM_,co))/(exp(2.5-0.1*sol(_VM_,co))-1);
			double BetaHm = 4*exp(-1*sol(_VM_,co)/18);
			double AlphaHh = 0.07*exp(-1*sol(_VM_,co)/20);
			double BetaHh = 1/(exp(3-0.1*sol(_VM_,co))+1);
			/*long double AlphaHn = ((0.01*(sol(_VM_,co) + 55))/(1-exp(-0.1*(sol(_VM_,co)+55))));
			long double BetaHn  = ((0.125*exp(-0.0125*(sol(_VM_, co)+65))));
			long double AlphaHm = ((0.1*(sol(_VM_,co)+40))/(1-exp(-0.1*(sol(_VM_,co)+40))));
			long double BetaHm  = (4*exp(-0.0556*(sol(_VM_,co)+65)));
			long double AlphaHh = (0.07*exp(-0.05*(sol(_VM_,co)+65)));
			long double BetaHh  = (1/(1 + exp(-0.1*(sol(_VM_,co)+35))));*/
			//Ableitungen nach VM
			double AlphaHn1 = (0.01*(1-exp(-0.1*sol(_VM_,co) -5.5)) - (0.01*sol(_VM_,co) +0.55)*(0.1*exp(0.1*sol(_VM_,co)-5.5)))/(pow((1-exp(-0.1*sol(_VM_,co)-5.5)),2));
			double BetaHn1 = 0.125 * -0.0125 * (exp(0.0125*sol(_VM_,co) - 0.8125));
			double AlphaHm1 = ((0.1*(1-exp(-0.1*sol(_VM_, co) + 4))) - ((0.1*sol(_VM_, co) + 4)*(0.1*exp(-0.1*sol(_VM_, co) + 4 ))))/pow((1-exp(-0.1*sol(_VM_, co) + 4)),2);
			double BetaHm1 = 4 * -0.556*exp(-0.556*sol(_VM_, co) + 3.614);
			double AlphaHh1 = -0.05*0.07 *exp(0.05*sol(_VM_, co) + 3.25);
			double BetaHh1 = -1*(1+exp(-0.1*sol(_VM_,co) -3.5))/pow((1+exp(-0.1*sol(_VM_, co) -3.5)),2);
			// here we need f'
			/*J(_h_, co, _h_, co) += scv.volume() * AlphaHh * (1-sol(_h_,co)) + BetaHh*(sol(_h_,co));
			J(_m_, co, _m_, co) += scv.volume() * AlphaHm * (1-sol(_m_,co)) + BetaHm*(sol(_m_,co));
			J(_n_, co, _n_, co) += scv.volume() * AlphaHn * (1-sol(_n_,co)) + BetaHn*(sol(_n_,co));*/
			// kapa von 1 Wert nch volume
			//J(_VM_, co, _VM_, co) += scv.volume() * (1 * 0.36*pow(sol(_n_,co),4) + 1.20*pow(sol(_m_,co),3)*sol(_h_,co) + 0.003);

			/*J(_n_, co, _VM_, co) += scv.volume() *(AlphaHn1*sol(_n_,co) - BetaHn1*(sol(_n_, co)));
			J(_m_, co, _VM_, co) += scv.volume() *(AlphaHm1*sol(_m_, co) - BetaHm1*(sol(_m_, co)));
			J(_h_, co, _VM_, co) += scv.volume() *(AlphaHh1*sol(_h_, co) - BetaHh1*(sol(_h_, co)));*/
			/*double m_m = sol(_m_,co) - 1/(AlphaHm+BetaHm);
			double m_h = sol(_h_,co) - 1/(AlphaHh+BetaHh);
			double m_n = sol(_n_,co) - 1/(AlphaHn+BetaHn);*/
			J(_VM_, co, _n_, co) += scv.volume() * 0.36*4*pow(sol(_n_,co),3)*(sol(_VM_,co) + 77 );

			J(_VM_, co, _m_, co) += scv.volume() * 1.20*3*pow(sol(_m_,co),2)*sol(_h_,co)*(sol(_VM_, co) - 50);

			J(_VM_, co, _h_, co) += scv.volume() * ( 1.20*pow(sol(_m_,co),3)*(sol(_VM_, co)- 50));

			// here we need f' ableitungen nach einzelnen kompos
			/*J(_h_, co, _h_, co) += scv.volume() * m_h;
			J(_m_, co, _m_, co) += scv.volume() * m_m;
			J(_n_, co, _n_, co) += scv.volume() * m_n;*/
			// Ableigung nach h,m,n ohne Ableitung wird matrix singulär
			J(_h_, co, _h_, co) += scv.volume() * ((-AlphaHh - BetaHh));
			J(_m_, co, _m_, co) += scv.volume() * ((-AlphaHm - BetaHm));
			J(_n_, co, _n_, co) += scv.volume() * ((-AlphaHn - BetaHn));

			// Ableitungen nach VM von h,m,n
			/*J(_h_, co, _VM_, co) += scv.volume() * ;
			J(_m_, co, _VM_, co) += scv.volume() * ;
			J(_n_, co, _VM_, co) += scv.volume() * ;*/

			/*J(_n_, co, _n_, co) += scv.volume() * 0.36*4*pow(sol(_n_,co),3)*(sol(_VM_,co) + 77 );

			J(_m_, co, _m_, co) += scv.volume() * 1.20*3*pow(sol(_m_,co),2)*sol(_h_,co)*(sol(_VM_, co) - 50);

			J(_h_, co, _h_, co) += scv.volume() * ( 1.20*pow(sol(_m_,co),3)*(sol(_VM_, co)- 50));*/

			std::cout << "volume: " << scv.volume() << std::endl;
	}







	// erstmal hils var fuer kapa und diam
	double kapa = 1;
	double diam = 10;


	// fluesse hier rein
	std::cout << "add_jac_A_elem" << std::endl;



//	Diffusion and Velocity Term
	//if(m_imDiffusion.data_given() || m_imVelocity.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		// here we need diam and kapa questioning
		////////////////////////////////////////////////////
		// Diffusive Term
		////////////////////////////////////////////////////

				std::cout << "In Diffusion add jac" << std::endl;
			// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					std::cout << "In Diffusion add jac  for schleife" << std::endl;
				// 	Compute cabel Equatation-Diffusion Tensor times Gradient
					double cabel = kapa*2*diam*M_PI;
					MatVecMult(Dgrad, m_imDiffusion[ip], scvf.global_grad(sh));

				//	Compute flux at IP
					const number D_diff_flux = VecDot(Dgrad, scvf.normal());

				// 	Add flux term to local matrix // HIER MATRIXINDIZES!!!
					UG_ASSERT((scvf.from() < J.num_row_dof(_VM_)) && (scvf.to() < J.num_col_dof(_VM_)),
							  "Bad local dof-index on element with object-id " << elem->base_object_id()
							  << " with center: " << CalculateCenter(elem, vCornerCoords));
					std::cout << "cabel: " << cabel << std::endl;
					std::cout << "D_diff_flux: " << D_diff_flux << std::endl;
					J(_VM_, scvf.from(), _VM_, sh) -= D_diff_flux*cabel;
					J(_VM_, scvf.to()  , _VM_, sh) += D_diff_flux*cabel;

			}

		////////////////////////////////////////////////////
		// Convective Term (not used)
		////////////////////////////////////////////////////
			/*if(m_imVelocity.data_given())
			{
				std::cout << "In Velocity add jac" << std::endl;
			//	Add Flux contribution
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					const number D_conv_flux = convShape(ip, sh);

				//	Add flux term to local matrix
					// zu testzweicken null gesetzt
					J(_VM_, scvf.from(), _VM_, sh) += D_conv_flux;
					J(_VM_, scvf.to(),   _VM_, sh) -= D_conv_flux;
				}
			}*/

			// no explicit dependency on flux import
		}
	//	std::cout << "add Jac A l√§uft" << std::endl;
	}

//	UG_LOG("Local Matrix is: \n"<<J<<"\n");

////////////////////////////////////////////////////
// Reaction Term (using lumping)
////////////////////////////////////////////////////
	// not in for loop
	/*const LocalVectorTimeSeries& vLocSol = *this->local_time_solutions();
	std::cout<< "LocalVTS working! " << std::endl;
	std::cout<< vLocSol.size() << std::endl;
	std::cout<< "Loesungen in U: "<< u << std::endl;
	// Fehler

	const LocalVector& sol = vLocSol.solution(0);
	std::cout<< "LocalVec working! " << std::endl;

// 	loop Sub Control Volume (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);
		std::cout<< "scv Working! " << std::endl;
	// 	get associated node
		const int co = scv.node_id();
		std::cout<< "co Working! " << std::endl;





	// 	Add to local matrix
		std::cout << "before zeugs dass ich mache!" << std::endl;
		double AlphaHn = ((0.01*(sol(_VM_,co) + 55))/(1-exp(-0.1*(sol(_VM_,co)+55))));
		double BetaHn  = ((0.125*exp(-0.0125*(sol(_VM_, co)+65))));
		double AlphaHm = ((0.1*(sol(_VM_,co)+40))/(1-exp(-0.1*(sol(_VM_,co)+40))));
		double BetaHm  = (4*exp(-0.0556*(sol(_VM_,co)+65)));
		double AlphaHh = (0.07*exp(-0.05*(sol(_VM_,co)+65)));
		double BetaHh  = (1/(1 + exp(-0.1*(sol(_VM_,co)+35))));
		J(_h_, co, _h_, co) += scv.volume() * AlphaHh * (1-sol(_h_,co)) - BetaHh*(sol(_h_,co));
		J(_m_, co, _m_, co) += scv.volume() * AlphaHm * (1-sol(_m_,co)) - BetaHm*(sol(_m_,co));
		J(_n_, co, _n_, co) += scv.volume() * AlphaHn * (1-sol(_n_,co)) - BetaHn*(sol(_n_,co));
		J(_VM_, co, _VM_, co) += scv.volume() * 2;
		std::cout<< "Volumen: " << scv.volume() << std::endl;*/
	//}*/
	//std::cout << "second addj is not working" << std::endl;
	std::cout << "Jacobi: " << J << std::endl;
//	reaction term does not explicitly depend on the associated unknown function
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

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();
		// need something getting values


	// 	getting local values
		const LocalVector& sol = u;


		//std::cout << "shape_vec: " << scv.shape_vector() << std::endl;
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

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);
	number value1 = 0;
	number valueh = 0;
	number valuem = 0;
	number valuen = 0;

		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
			{
			// 	get current SCV
				const typename TFVGeom::SCV& scv = geo.scv(ip);
				// 	get associated node
				const int co = scv.node_id();
				// 	Add to local defect

				std::cout << " fehler kommt gleich" << std::endl;
				const LocalVectorTimeSeries& vLocSol = *this->local_time_solutions();
				const LocalVector& sol = vLocSol.solution(0);
				const LocalVector& solnew = vLocSol.solution(1);

			//	remember local solutions
				number dt = vLocSol.time(0) - vLocSol.time(1);
				number time = vLocSol.time(0);

				// umwandlung der positionen muss für 3d noch geändert werden
				std::ostringstream stream_x;
				stream_x << vCornerCoords[0];
				std::string String_x = stream_x.str();
				std::cout << "tests" << std::endl;
				std::cout << String_x << std::endl;
				size_t end = String_x.find(")",0);
				size_t beg = 1;
				std::string String_xx = String_x.substr(beg, end-1);
				std::cout << String_xx << std::endl;
				// hier kommt fehler
				std::stringstream String_xxx;
				number x;
				String_xxx << String_xx;
				String_xxx >> x;
				std::cout << String_xxx << std::endl;
				std::cout << x << std::endl;

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
				double AlphaHn = (0.1-0.01*sol(_VM_,co))/(exp(1-0.1*(sol(_VM_,co)))-1);
				double BetaHn = 0.125*exp(sol(_VM_,co)/80);
				double AlphaHm = (2.5 - 0.1*sol(_VM_,co))/(exp(2.5-0.1*sol(_VM_,co))-1);
				double BetaHm = 4*exp(-1*sol(_VM_,co)/18);
				double AlphaHh = 0.07*exp(-1*sol(_VM_,co)/20);
				double BetaHh = 1/(exp(3-0.1*sol(_VM_,co))+1);
				// hier hab ich zuletzte dran gearbeitet...

				double flux_m = ((AlphaHm * (1-sol(_m_,co)) - BetaHm * sol(_m_,co)));
				double flux_h = ((AlphaHh * (1-sol(_h_,co)) - BetaHh * sol(_h_,co)));
				double flux_n = ((AlphaHn * (1-sol(_n_,co)) - BetaHn * sol(_n_,co)));

				//std::cout << "m: " << m_m << "h: " << m_h << "n: " << m_n << std::endl;
				//std::cout << sol(_VM_,co) << std::endl;
				std::cout << "A/B-Hn: " << AlphaHn <<" " << BetaHn<< std::endl;
				std::cout << "A/B-Hm: " << AlphaHm <<" " << BetaHm<< std::endl;
				std::cout << "A/B-Hh: " << AlphaHh <<" " << BetaHh<< std::endl;
				//d(_h_, co) += scv.volume() * AlphaHh * (1-sol(_h_,co)) - BetaHh*(sol(_h_,co));


				//std::cout << "H " << (scv.volume() * AlphaHm * (1-sol(_m_,co)) - BetaHm*(sol(_m_,co))) << std::endl;

				//d(_m_, co) += scv.volume() * AlphaHm * (1-sol(_m_,co)) - BetaHm*(sol(_m_,co));
				//d(_n_, co) += scv.volume() * AlphaHn * (1-sol(_n_,co)) - BetaHn*(sol(_n_,co));
				// capazit√§t erstmal nicht benutzen in diesem fall nur im ersten wert
				//const number capacitive_part_of_flux = m_capacity * ( sol(_VM_, co) - oldSol(_VM_, co) ) / 0.01 ;

				//const number capacitive_part_of_flux = 1 * ( solnew(_VM_, co) - sol(_VM_, co) ) / dt ;
				const number potassium_part_of_flux = 0.36*pow(sol(_n_,co),4)*(sol(_VM_,co) + 77);
				const number sodium_part_of_flux =  1.20*pow(sol(_m_,co),3)*sol(_h_,co) * (sol(_VM_, co) - 50);
				const number leakage_part_of_flux = 0.003*(sol(_VM_,co) + 54.4);
				number inject = 0;
				// laesst sich spaeter mit for schleife bezug auf dim loesen
				std::cout << "Eig-XWert: " << vCornerCoords[0] << std::endl;
				std::cout << "Eig-yWert: " << vCornerCoords[1] << std::endl;
				std::cout << "Eig-zWert: " << vCornerCoords[2] << std::endl;
				std::cout << "X-Wert: " << x << std::endl;
				(*m_Injection)(inject, 2, time, x);


				const number flux =   (//capacitive_part_of_flux +
									 potassium_part_of_flux
									+ sodium_part_of_flux
									+ leakage_part_of_flux//);
									- inject);
				std::cout<< "Injection: " << inject << std::endl;
				// fehler in defekt normal -inject/radius
				d(_VM_, co) += scv.volume()*((flux)*((sol(_VM_,co)-(sol(_VM_, co)-(flux))))-inject); // * scv.volume; //scv.volume() * flux * 0.001;
				//value1 = flux;
				//std::cout << "defekt: " << d(_VM_, co) << std::endl;
				d(_h_, co) += scv.volume()*flux_h;//*((sol(_h_,co)-sol(_h_,co)-flux_h));//-valueh);
				d(_m_, co) += scv.volume()*flux_m;//flux_m;//*((sol(_m_,co)-sol(_m_,co)-flux_m));//-valuem);//0;
				d(_n_, co) += scv.volume()*flux_n;//*((sol(_n_,co)-sol(_n_,co)-flux_n));//-valuen);
				//
				//std::cout << "volumen*flux: " << (scv.volume() *flux) << std::endl;
				number valueh = flux_h;
				number valuem = flux_m;
				number valuen = flux_n;
				std::cout << "flux: " << flux << std::endl;
				std::cout << "defekt VM: " << d(_VM_, co) << "bei Loesung: " << sol(_VM_,co) << std::endl;
				std::cout << "defekt h: " << d(_h_, co) << "bei Loesung h: " << sol(_h_,co)<< std::endl;
				std::cout << "defekt m: " << d(_m_, co) << "bei Loesung m: " << sol(_m_,co)<< std::endl;
				std::cout << "defekt n: " << d(_n_, co) << "bei Loesung n: " << sol(_n_,co)<< std::endl;
			}
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

			//	scale by diffusion tensor
				MatVecMult(Dgrad_c, m_imDiffusion[ip], grad_c);

			// 	Compute flux
				const number diff_flux = VecDot(Dgrad_c, scvf.normal());
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
add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	/*static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reaction rate
	if(m_imReactionRateExpl.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_VM_, co) += u(_VM_, co) * m_imReactionRateExpl[ip] * scv.volume();
		}
	}

//	reaction
	if(m_imReactionExpl.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_VM_, co) += m_imReactionExpl[ip] * scv.volume();
		}
	}

	if(m_imSourceExpl.data_given())
	{
		// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
			// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

			// 	get associated node
			const int co = scv.node_id();

			// 	Add to local rhs
			d(_VM_, co) -= m_imSourceExpl[ip] * scv.volume();
		}
	}*/
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	std::cout << "add_def_M_elem" << std::endl;
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	//if(!m_imMassScale.data_given() && !m_imMass.data_given()) return;

// 	loop Sub Control Volumes (SCV)


}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	std::cout << "add_rhs_elem" << std::endl;

	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	/*
	// loop Sub Control Volumes (SCV)
	// Quellterme fuer spannung
	// source daten werden nicht gegeben sondern hier berechnet!
	//if ( m_imSource.data_given() ) {*/
		for ( size_t ip = 0; ip < geo.num_scv(); ++ip ) {
			// get current SCV
			const typename TFVGeom::SCV& scv = geo.scv( ip );

			// get associated node
			const int co = scv.node_id();
			// help vars
			const LocalVectorTimeSeries& vLocSol = *this->local_time_solutions();
			const LocalVector& sol = vLocSol.solution(0);

			/*double AlphaHn = ((0.01*(sol(_VM_,co) + 55))/(1-exp(-0.1*(sol(_VM_,co)+55))));
			double BetaHn  = ((0.125*exp(-0.0125*(sol(_VM_, co)+65))));
			double AlphaHm = ((0.1*(sol(_VM_,co)+40))/(1-exp(-0.1*(sol(_VM_,co)+40))));
			double BetaHm  = (4*exp(-0.0556*(sol(_VM_,co)+65)));
			double AlphaHh = (0.07*exp(-0.05*(sol(_VM_,co)+65)));
			double BetaHh  = (1/(1 + exp(-0.1*(sol(_VM_,co)+35))));
			// Add to local rhs here you could give some source data for VM
			d(_h_, co) += scv.volume() * AlphaHh * (1-sol(_h_,co)) - BetaHh*(sol(_h_,co));
			d(_m_, co) += scv.volume() * AlphaHm * (1-sol(_m_,co)) - BetaHm*(sol(_m_,co));
			d(_n_, co) += scv.volume() * AlphaHn * (1-sol(_n_,co)) - BetaHn*(sol(_n_,co));
			std::cout<< "first source data given" << std::endl;
			std::cout << "h: "<< d(_h_, co) << std::endl;
			//d(_VM_, co) +=  * scv.volume();
		//}
	}

	// loop Sub Control Volumes (SCVF)
	//if ( m_imVectorSource.data_given() ) {
		double test = -0.1;
		ug::MathVector<dim> v;
		// 1er vector bauen
		v = test;*/


			if ( m_imVectorSource.data_given() ) {
				for ( size_t ip = 0; ip < geo.num_scvf(); ++ip ) {
					// get current SCVF
					const typename TFVGeom::SCVF& scvf = geo.scvf( ip );

					// Add to local rhs
					d(_VM_, scvf.from()) -= VecDot(m_imVectorSource[ip], scvf.normal() );
					d(_VM_, scvf.to()  ) += VecDot(m_imVectorSource[ip], scvf.normal() );
				}
			}
	}
}


////////////////////////////////////
///   error estimation (begin)   ///

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	std::cout << "prep_err_est_elem_loop" << std::endl;
	//	get the error estimator data object and check that it is of the right type
	//	we check this at this point in order to be able to dispense with this check later on
	//	(i.e. in prep_err_est_elem and compute_err_est_A_elem())
	if (this->m_spErrEstData.get() == NULL)
	{
		UG_THROW("No ErrEstData object has been given to this ElemDisc!");
	}

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (!err_est_data)
	{
		UG_THROW("Dynamic cast to SideAndElemErrEstData failed."
				<< std::endl << "Make sure you handed the correct type of ErrEstData to this discretization.");
	}


//	check that upwind has been set
	if (m_spConvShape.invalid())
		UG_THROW("ConvectionDiffusionFV1::prep_err_est_elem_loop: "
				 "Upwind has not been set.");

//	set local positions
	if (!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		// get local IPs
		size_t numSideIPs, numElemIPs;
		const MathVector<refDim>* sideIPs;
		const MathVector<refDim>* elemIPs;
		try
		{
			numSideIPs = err_est_data->num_all_side_ips(roid);
			numElemIPs = err_est_data->num_elem_ips(roid);
			sideIPs = err_est_data->template side_local_ips<refDim>(roid);
			elemIPs = err_est_data->template elem_local_ips<refDim>(roid);

			if (!sideIPs || !elemIPs) return;	// are NULL if TElem is not of the same dim as TDomain
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		// set local IPs in imports
		m_imDiffusion.template 		set_local_ips<refDim>(sideIPs, numSideIPs, false);
		m_imVelocity.template 		set_local_ips<refDim>(sideIPs, numSideIPs, false);
		m_imFlux.template 			set_local_ips<refDim>(sideIPs, numSideIPs, false);
		m_imSource.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
		m_imVectorSource.template 	set_local_ips<refDim>(sideIPs, numSideIPs, false);
		m_imReactionRate.template 	set_local_ips<refDim>(elemIPs, numElemIPs, false);
		m_imReaction.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
		m_imMassScale.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
		m_imMass.template 			set_local_ips<refDim>(elemIPs, numElemIPs, false);

		//	init upwind for element type
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		if (!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
			UG_THROW("ConvectionDiffusionFV1::prep_err_est_elem_loop: "
					 "Cannot init upwind for element type.");

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numElemIPs, numSideIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numElemIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtElemIP(ip), elemIPs[ip]);
		for (size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());
	std::cout << "prep_err_est_elem" << std::endl;
// 	update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try
	{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("ConvectionDiffusionFV1::prep_err_est_elem: Cannot update Finite Volume Geometry.");

//	roid
	ReferenceObjectID roid = elem->reference_object_id();

//	set local positions
	if (TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		size_t numSideIPs, numElemIPs;
		const MathVector<refDim>* sideIPs;
		const MathVector<refDim>* elemIPs;
		try
		{
			numSideIPs = err_est_data->num_all_side_ips(roid);
			numElemIPs = err_est_data->num_elem_ips(roid);
			sideIPs = err_est_data->template side_local_ips<refDim>(roid);
			elemIPs = err_est_data->template elem_local_ips<refDim>(roid);

			if (!sideIPs || !elemIPs) return;	// are NULL if TElem is not of the same dim as TDomain
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		m_imDiffusion.template 		set_local_ips<refDim>(sideIPs, numSideIPs);
		m_imVelocity.template 		set_local_ips<refDim>(sideIPs, numSideIPs);
		m_imFlux.template 			set_local_ips<refDim>(sideIPs, numSideIPs);
		m_imSource.template 		set_local_ips<refDim>(elemIPs, numElemIPs);
		m_imVectorSource.template 	set_local_ips<refDim>(sideIPs, numSideIPs);
		m_imReactionRate.template 	set_local_ips<refDim>(elemIPs, numElemIPs);
		m_imReaction.template 		set_local_ips<refDim>(elemIPs, numElemIPs);
		m_imMassScale.template 		set_local_ips<refDim>(elemIPs, numElemIPs);
		m_imMass.template 			set_local_ips<refDim>(elemIPs, numElemIPs);

		//	init upwind for element type
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		if (!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
			UG_THROW("ConvectionDiffusionFV1::prep_err_est_elem_loop: "
					 "Cannot init upwind for element type.");

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numElemIPs, numSideIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numElemIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtElemIP(ip), elemIPs[ip]);
		for (size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}

//	set global positions
	size_t numSideIPs, numElemIPs;
	std::vector<MathVector<dim> > sideIPs;
	std::vector<MathVector<dim> > elemIPs;

	try
	{
		numSideIPs = err_est_data->num_all_side_ips(roid);
		numElemIPs = err_est_data->num_elem_ips(roid);
		sideIPs = std::vector<MathVector<dim> >(numSideIPs);
		elemIPs = std::vector<MathVector<dim> >(numElemIPs);

		err_est_data->all_side_global_ips(&sideIPs[0], elem, vCornerCoords);
		err_est_data->elem_global_ips(&elemIPs[0], elem, vCornerCoords);
	}
	UG_CATCH_THROW("Global integration points for error estimator cannot be set.");

	m_imDiffusion.			set_global_ips(&sideIPs[0], numSideIPs);
	m_imVelocity.			set_global_ips(&sideIPs[0], numSideIPs);
	m_imFlux.				set_global_ips(&sideIPs[0], numSideIPs);
	m_imSource.				set_global_ips(&elemIPs[0], numElemIPs);
	m_imVectorSource.		set_global_ips(&sideIPs[0], numSideIPs);
	m_imReactionRate.		set_global_ips(&elemIPs[0], numElemIPs);
	m_imReaction.			set_global_ips(&elemIPs[0], numElemIPs);
	m_imMassScale.			set_global_ips(&elemIPs[0], numElemIPs);
	m_imMass.				set_global_ips(&elemIPs[0], numElemIPs);
}

//	computes the error estimator contribution (stiffness part) for one element
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	std::cout << "compute_err_est_A_elem" << std::endl;
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

//	request geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

////////////////
// SIDE TERMS //
////////////////

//	get the sides of the element
	//	We have to cast elem to a pointer of type SideAndElemErrEstData::elem_type
	//	for the SideAndElemErrEstData::operator() to work properly.
	//	This cannot generally be achieved by casting to TElem*, since this method is also registered for
	//	lower-dimensional types TElem, and must therefore be compilable, even if it is never EVER to be executed.
	//	The way we achieve this here, is by calling associated_elements_sorted() which has an implementation for
	//	all possible types. Whatever comes out of it is of course complete nonsense if (and only if)
	//	SideAndElemErrEstData::elem_type != TElem. To be on the safe side, we throw an error if the number of
	//	entries in the list is not as it should be.

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::side_type>::secure_container side_list;
	pErrEstGrid->associated_elements_sorted(side_list, (TElem*) elem);
	if (side_list.size() != (size_t) ref_elem_type::numSides)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

// 	some help variables
	MathVector<dim> fluxDensity, gradC, normal;

// calculate grad u (take grad from first scvf ip (grad u is constant on the entire element))
	if (geo.num_scvf() < 1) {UG_THROW("Element has no SCVFs!");}
	const typename TFVGeom::SCVF& scvf = geo.scvf(0);

	VecSet(gradC, 0.0);
	for (size_t j=0; j<m_shapeValues.num_sh(); j++)
		VecScaleAppend(gradC, u(_VM_,j), scvf.global_grad(j));

// calculate flux through the sides
	size_t passedIPs = 0;
	for (size_t side=0; side < (size_t) ref_elem_type::numSides; side++)
	{
		// normal on side
		SideNormal<ref_elem_type,dim>(normal, side, vCornerCoords);
		VecNormalize(normal, normal);

		try
		{
			for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
			{
				size_t ip = passedIPs + sip;

				VecSet(fluxDensity, 0.0);

			////// diffusion //////
				if (m_imDiffusion.data_given())
					MatVecScaleMultAppend(fluxDensity, -1.0, m_imDiffusion[ip], gradC);

			////// convection //////
				if (m_imVelocity.data_given())
				{
					number val = 0.0;
					for (size_t sh = 0; sh < m_shapeValues.num_sh(); sh++)
						val += u(_VM_,sh) * m_shapeValues.shapeAtSideIP(sh,sip);

					VecScaleAppend(fluxDensity, val, m_imVelocity[ip]);
				}

			////// general flux //////
				if (m_imFlux.data_given())
					VecAppend(fluxDensity, m_imFlux[ip]);

				(*err_est_data)(side_list[side],sip) += scale * VecDot(fluxDensity, normal);
			}

			passedIPs += err_est_data->num_side_ips(side_list[side]);
		}
		UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
				<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
	}

//////////////////
// VOLUME TERMS //
//////////////////

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
		{
			number total = 0.0;

		////// diffusion //////	TODO ONLY FOR (PIECEWISE) CONSTANT DIFFUSION TENSOR SO FAR!
		// div(D*grad(c)) = div(v)*u + v*grad(c)
		// nothing to do, as u is piecewise linear and div(D*grad(c)) disappears

		////// convection ////// TODO ONLY FOR CONSTANT VELOCITY FIELDS SO FAR!
		// div(v*c) = div(v)*u + v*grad(c) -- gradC has been calculated above
			if (m_imVelocity.data_given())
				total += VecDot(m_imVelocity[ip], gradC);

		////// general flux ////// TODO ONLY FOR DIVERGENCE-FREE FLUX FIELD SO FAR!
		// nothing to do

		////// reaction //////
			if (m_imReactionRate.data_given())
			{
				number val = 0.0;
				for (size_t sh = 0; sh < geo.num_sh(); sh++)
					val += u(_VM_,sh) * m_shapeValues.shapeAtElemIP(sh,ip);

				total += m_imReactionRate[ip] * val;
			}

			if (m_imReaction.data_given())
			{
				total += m_imReaction[ip];
			}

			(*err_est_data)(elem_list[0],ip) += scale * total;
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	computes the error estimator contribution (mass part) for one element
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
// note: mass parts only enter volume term
	std::cout << "compute_err_est_M_elem" << std::endl;
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

//	request geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop integration points
	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
		{
			number total = 0.0;

		////// mass scale //////
			if (m_imMassScale.data_given())
			{
				number val = 0.0;
				for (size_t sh = 0; sh < geo.num_sh(); sh++)
					val += u(_VM_,sh) * m_shapeValues.shapeAtElemIP(sh,ip);

				total += m_imMassScale[ip] * val;
			}

		////// mass //////
			if (m_imMass.data_given())
			{
				total += m_imMass[ip];
			}

			(*err_est_data)(elem_list[0],ip) += scale * total;
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	computes the error estimator contribution (rhs part) for one element
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	std::cout << "compute_err_est_rhs_elem" << std::endl;
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

////////////////
// SIDE TERMS //
////////////////
//	get the sides of the element
	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::side_type>::secure_container side_list;
	pErrEstGrid->associated_elements_sorted(side_list, (TElem*) elem);
	if (side_list.size() != (size_t) ref_elem_type::numSides)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

// loop sides
	size_t passedIPs = 0;
	for (size_t side = 0; side < (size_t) ref_elem_type::numSides; side++)
	{
		// normal on side
		MathVector<dim> normal;
		SideNormal<ref_elem_type,dim>(normal, side, vCornerCoords);
		VecNormalize(normal, normal);

		try
		{
			for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
			{
				size_t ip = passedIPs + sip;

			////// vector source //////
				if (m_imVectorSource.data_given())
					(*err_est_data)(side_list[side],sip) += scale * VecDot(m_imVectorSource[ip], normal);
			}

			passedIPs += err_est_data->num_side_ips(side_list[side]);
		}
		UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
				<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
	}

//////////////////
// VOLUME TERMS //
//////////////////
	if (!m_imSource.data_given()) return;

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

////// source //////
	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
			(*err_est_data)(elem_list[0],ip) += scale * m_imSource[ip];
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
fsh_err_est_elem_loop()
{
	std::cout << "fsh_err_est_elem_loop" << std::endl;
//	finish the element loop in the same way as the actual discretization
	this->template fsh_elem_loop<TElem, TFVGeom> ();
};

///   error estimation (end)     ///
////////////////////////////////////

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_velocity(const LocalVector& u,
                 std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                 const size_t nip)
{
	std::cout << "lin_def_velocity" << std::endl;
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	//	sum up contributions of convection shapes
		MathVector<dim> linDefect;
		VecSet(linDefect, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(linDefect, u(_VM_,sh), convShape.D_vel(ip, sh));

	//	add parts for both sides of scvf
		vvvLinDef[ip][_VM_][scvf.from()] += linDefect;
		vvvLinDef[ip][_VM_][scvf.to()] -= linDefect;
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_diffusion(const LocalVector& u,
                  std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                  const size_t nip)
{
	std::cout << "lin_def_diffusion" << std::endl;
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	// 	compute gradient at ip
		MathVector<dim> grad_u;	VecSet(grad_u, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad_u, u(_VM_,sh), scvf.global_grad(sh));

	//	compute the lin defect at this ip
		MathMatrix<dim,dim> linDefect;

	//	part coming from -\nabla u * \vec{n}
		for(size_t k=0; k < (size_t)dim; ++k)
			for(size_t j = 0; j < (size_t)dim; ++j)
				linDefect(j,k) = (scvf.normal())[j] * grad_u[k];

	//	add contribution from convection shapes
		if(convShape.non_zero_deriv_diffusion())
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				MatAdd(linDefect, convShape.D_diffusion(ip, sh), u(_VM_, sh));

	//	add contributions
		vvvLinDef[ip][_VM_][scvf.from()] -= linDefect;
		vvvLinDef[ip][_VM_][scvf.to()  ] += linDefect;
	}
}

template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_flux(const LocalVector& u,
             std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
             const size_t nip)
{
	std::cout << "lin_def_flux" << std::endl;
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	//	add parts for both sides of scvf
		vvvLinDef[ip][_VM_][scvf.from()] += scvf.normal();
		vvvLinDef[ip][_VM_][scvf.to()] -= scvf.normal();
	}
}
//	computes the linearized defect w.r.t to the reaction rate
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_reaction_rate(const LocalVector& u,
                      std::vector<std::vector<number> > vvvLinDef[],
                      const size_t nip)
{
	std::cout << "lin_def_reaction_rate" << std::endl;
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_VM_][co] = u(_VM_, co) * scv.volume();
	}
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_reaction(const LocalVector& u,
                 std::vector<std::vector<number> > vvvLinDef[],
                 const size_t nip)
{
	std::cout << "lin_def_reaction" << std::endl;
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_VM_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_source(const LocalVector& u,
               std::vector<std::vector<number> > vvvLinDef[],
               const size_t nip)
{
	std::cout << "lin_def_source" << std::endl;
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_VM_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the vector source
//	(in analogy to velocity)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_vector_source(const LocalVector& u,
                      std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                      const size_t nip)
{
	std::cout << "lin_def_vector_source" << std::endl;
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

	// loop Sub Control Volumes Faces (SCVF)
	for ( size_t ip = 0; ip < geo.num_scvf(); ++ip ) {
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf( ip );

		// add parts for both sides of scvf
		vvvLinDef[ip][_VM_][scvf.from()] -= scvf.normal();
		vvvLinDef[ip][_VM_][scvf.to()] += scvf.normal();
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_mass_scale(const LocalVector& u,
                   std::vector<std::vector<number> > vvvLinDef[],
                   const size_t nip)
{
	std::cout << "lin_def_mass_scale" << std::endl;
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t co = 0; co < geo.num_scv(); ++co)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(co);

	// 	Check associated node
		UG_ASSERT(co == scv.node_id(), "Only one shape per SCV");

	// 	set lin defect
		vvvLinDef[co][_VM_][co] = u(_VM_, co) * scv.volume();
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ElemDiscHH_FV1<TDomain>::
lin_def_mass(const LocalVector& u,
             std::vector<std::vector<number> > vvvLinDef[],
             const size_t nip)
{
	std::cout << "lin_def_mass" << std::endl;
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t co = 0; co < geo.num_scv(); ++co)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(co);

	// 	Check associated node
		UG_ASSERT(co == scv.node_id(), "Only one shape per SCV");

	// 	set lin defect
		vvvLinDef[co][_VM_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the velocity
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

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
const typename ElemDiscHH_FV1<TDomain>::conv_shape_type&
ElemDiscHH_FV1<TDomain>::
get_updated_conv_shapes(const FVGeometryBase& geo)
{
	std::cout << "get_updated_conv_shapes" << std::endl;
//	compute upwind shapes for transport equation
//	\todo: we should move this computation into the preparation part of the
//			disc, to only compute the shapes once, reusing them several times.
	if(m_imVelocity.data_given())
	{
	//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		if(m_imDiffusion.data_given()) vDiffusion = m_imDiffusion.values();

	//	update convection shapes
		if(!m_spConvShape->update(&geo, m_imVelocity.values(), vDiffusion, true))
		{
			UG_LOG("ERROR in 'ElemDiscHHFV1::add_jac_A_elem': "
					"Cannot compute convection shapes.\n");
		}
	}

//	return a const (!!) reference to the upwind
	return *const_cast<const IConvectionShapes<dim>*>(m_spConvShape.get());
}



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
	this->set_add_def_A_expl_elem_fct(id, &T::template add_def_A_expl_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TFVGeom>);

// error estimator parts
	this->set_prep_err_est_elem_loop(id, &T::template prep_err_est_elem_loop<TElem, TFVGeom>);
	this->set_prep_err_est_elem(id, &T::template prep_err_est_elem<TElem, TFVGeom>);
	this->set_compute_err_est_A_elem(id, &T::template compute_err_est_A_elem<TElem, TFVGeom>);
	this->set_compute_err_est_M_elem(id, &T::template compute_err_est_M_elem<TElem, TFVGeom>);
	this->set_compute_err_est_rhs_elem(id, &T::template compute_err_est_rhs_elem<TElem, TFVGeom>);
	this->set_fsh_err_est_elem_loop(id, &T::template fsh_err_est_elem_loop<TElem, TFVGeom>);

//	set computation of linearized defect w.r.t velocity
	m_imDiffusion.set_fct(id, this, &T::template lin_def_diffusion<TElem, TFVGeom>);
	m_imVelocity. set_fct(id, this, &T::template lin_def_velocity<TElem, TFVGeom>);
	m_imFlux.set_fct(id, this, &T::template lin_def_flux<TElem, TFVGeom>);
	m_imReactionRate. set_fct(id, this, &T::template lin_def_reaction_rate<TElem, TFVGeom>);
	m_imReaction. set_fct(id, this, &T::template lin_def_reaction<TElem, TFVGeom>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source<TElem, TFVGeom>);
	m_imVectorSource.set_fct(id, this, &T::template lin_def_vector_source<TElem, TFVGeom>);
	m_imMassScale.set_fct(id, this, &T::template lin_def_mass_scale<TElem, TFVGeom>);
	m_imMass.	set_fct(id, this, &T::template lin_def_mass<TElem, TFVGeom>);

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

