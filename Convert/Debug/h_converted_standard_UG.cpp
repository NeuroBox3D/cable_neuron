#include "h_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::alpt(double v) 
{ 
  return  exp(0.0378*zetat*(v-vhalft)) ; 
}
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::bett(double v) 
{ 
  return  exp(0.0378*zetat*gmt*(v-vhalft)) ; 
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::approx_space_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getehd() 
{ 
return ehd; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getghdbar() 
{ 
return ghdbar; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getvhalfl() 
{ 
return vhalfl; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getkl() 
{ 
return kl; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getvhalft() 
{ 
return vhalft; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::geta0t() 
{ 
return a0t; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getzetat() 
{ 
return zetat; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getgmt() 
{ 
return gmt; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getq10() 
{ 
return q10; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getqtl() 
{ 
return qtl; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setehd(double val) 
{ 
ehd = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setghdbar(double val) 
{ 
ghdbar = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setvhalfl(double val) 
{ 
vhalfl = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setkl(double val) 
{ 
kl = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setvhalft(double val) 
{ 
vhalft = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::seta0t(double val) 
{ 
a0t = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setzetat(double val) 
{ 
zetat = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setgmt(double val) 
{ 
gmt = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setq10(double val) 
{ 
q10 = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setqtl(double val) 
{ 
qtl = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->lGate)) 
UG_THROW("Attachment necessary (lGate) for h_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->lGate); 
this->aalGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->lGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> h_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
{ 
	 //var for output 
	 std::vector<number> GatingAccesors; 
 
	 typedef ug::MathVector<TDomain::dim> position_type; 
 
	 position_type coord; 
 
	 if (coord.size()==1) 
	 	 coord[0]=x; 
	 if (coord.size()==2) 
	 { 
	 	 coord[0] = x;
	 	 coord[1] = y;
	 } 
	 if (coord.size()==3) 
	 { 
	 	 coord[0] = x;
	 	 coord[1] = y;
	 	 coord[2] = z;
	 } 
	 //accesors 
	 typedef Attachment<position_type> position_attachment_type; 
	 typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accesor_type; 
 
	 // Definitions for Iteration over all Elements 
	 typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; 
	 SubsetGroup ssGrp; 
	 try { ssGrp = SubsetGroup(m_pVMDisc->approx_space()->domain()->subset_handler(), this->m_vSubset);} 
	 UG_CATCH_THROW("Subset group creation failed."); 
 
	 itType iter; 
	 number bestDistSq, distSq; 
	 Vertex* bestVrt; 
 
	 // Iterate only if there is one Gtting needed 
	 if (m_log_lGate == true )
	 { 
	 	 // iterating over all elements 
	 	 for (size_t si=0; si < ssGrp.size(); si++) 
	 	 { 
	 	 	 itType iterBegin = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]); 
	 	 	 itType iterEnd = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
	 	 	 const position_accesor_type& aaPos = m_pVMDisc->approx_space()->domain()->position_accessor(); 
	 	 	 if (si==0) 
	 	 	 { 
	 	 	 	 bestVrt = *iterBegin; 
	 	 	 	 bestDistSq = VecDistanceSq(coord, aaPos[bestVrt]); 
	 	 	 } 
	 	 	 iter = iterBegin; 
	 	 	 iter++; 
	 	 	 while(iter != iterEnd) 
	 	 	 { 
	 	 	 	 distSq = VecDistanceSq(coord, aaPos[*iter]); 
	 	 	 	 { 
	 	 	 	 	 bestDistSq = distSq; 
	 	 	 	 	 bestVrt = *iter; 
	 	 	 	 } 
	 	 	 	 ++iter; 
	 	 	 } 
	 	 } 
	 	 if (m_log_lGate == true) 
	 	 	 GatingAccesors.push_back(this->aalGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void h_converted_standard_UG<TDomain>::set_log_lGate(bool bLoglGate) { m_log_lGate = bLoglGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
//get celsius and time
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
number dt = m_pVMDisc->time(); 
// make preparing vor getting values of every edge 
number v = vrt_values[VMDisc<TDomain>::_v_]; 

 
double          a,qt; 
qt= pow(q10 , ((celsius-33)/10)); 
        a = alpt(v); 
        linf = 1/(1 + exp(-(v-vhalfl)/kl)); 
;//       linf = 1/(1+ alpl(v))
double         taul = bett(v)/(qtl*qt*a0t*(1+a)); 
aalGate[vrt] =linf; 
}  
 
 
 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
 number FARADAY = m_pVMDisc->F; 
 number dt = newTime - m_pVMDisc->time(); 
number v = vrt_values[VMDisc<TDomain>::_v_]; 

 
double l = aalGate[vrt]; 

 
 
double          a,qt; 
qt= pow(q10 , ((celsius-33)/10)); 
        a = alpt(v); 
        linf = 1/(1 + exp(-(v-vhalfl)/kl)); 
//       linf = 1/(1+ alpl(v))
double         taul = bett(v)/(qtl*qt*a0t*(1+a)); 
        l  +=   (linf - l)/taul*dt; 
; 
 

 
 
aalGate[vrt] = l; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number l = aalGate[ver]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 
number ghd = ghdbar*l; 

 
 
outCurrentValues.push_back( ghd*(v-ehd)); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class h_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class h_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class h_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


