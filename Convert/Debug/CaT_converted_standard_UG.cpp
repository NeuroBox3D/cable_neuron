#include "CaT_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
 	F = m_pVMDisc->F; 
 R = m_pVMDisc->R; 
 K = m_pVMDisc->temperature(); 
 celsius = m_pVMDisc->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getgbar() 
{ 
return gbar; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvshift() 
{ 
return vshift; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getcao() 
{ 
return cao; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getcai() 
{ 
return cai; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvmin() 
{ 
return vmin; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvmax() 
{ 
return vmax; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getv12m() 
{ 
return v12m; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getv12h() 
{ 
return v12h; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvwm() 
{ 
return vwm; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvwh() 
{ 
return vwh; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getam() 
{ 
return am; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getah() 
{ 
return ah; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvm1() 
{ 
return vm1; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvm2() 
{ 
return vm2; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvh1() 
{ 
return vh1; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getvh2() 
{ 
return vh2; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getwm1() 
{ 
return wm1; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getwm2() 
{ 
return wm2; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getwh1() 
{ 
return wh1; 
} 
template<typename TDomain> 
double CaT_converted_standard_UG<TDomain>::getwh2() 
{ 
return wh2; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setgbar(double val) 
{ 
gbar = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvshift(double val) 
{ 
vshift = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setcao(double val) 
{ 
cao = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setcai(double val) 
{ 
cai = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvmin(double val) 
{ 
vmin = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvmax(double val) 
{ 
vmax = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setv12m(double val) 
{ 
v12m = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setv12h(double val) 
{ 
v12h = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvwm(double val) 
{ 
vwm = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvwh(double val) 
{ 
vwh = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setam(double val) 
{ 
am = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setah(double val) 
{ 
ah = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvm1(double val) 
{ 
vm1 = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvm2(double val) 
{ 
vm2 = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvh1(double val) 
{ 
vh1 = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setvh2(double val) 
{ 
vh2 = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setwm1(double val) 
{ 
wm1 = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setwm2(double val) 
{ 
wm2 = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setwh1(double val) 
{ 
wh1 = val; 
} 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::setwh2(double val) 
{ 
wh2 = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for CaT_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for CaT_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> CaT_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 if (m_log_mGate || m_log_hGate )
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
	 	 if (m_log_mGate == true) 
	 	 	 GatingAccesors.push_back(this->aamGate[bestVrt]); 
	 	 if (m_log_hGate == true) 
	 	 	 GatingAccesors.push_back(this->aahGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void CaT_converted_standard_UG<TDomain>::set_log_mGate(bool bLogmGate) { m_log_mGate = bLogmGate; }
template<typename TDomain> void CaT_converted_standard_UG<TDomain>::set_log_hGate(bool bLoghGate) { m_log_hGate = bLoghGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
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
number v = vrt_values[CableEquation<TDomain>::_v_]; 
number ca = vrt_values[CableEquation<TDomain>::_ca_]; 

 
double          tinc; 
double           a, b; 
double v_ = v; 
double 	minf = 1.0 / ( 1 + exp(-(v_+v12m)/vwm) ); 
double 	hinf = 1.0 / ( 1 + exp((v_+v12h)/vwh) ); 
double 	mtau = ( am + 1.0 / ( exp((v_+vm1)/wm1) + exp(-(v_+vm2)/wm2) ) ) ; 
double 	htau = ( ah + 1.0 / ( exp((v_+vh1)/wh1) + exp(-(v_+vh2)/wh2) ) ) ; 
        tinc = -dt ; 
double         mexp = 1 - exp(tinc/mtau); 
double         hexp = 1 - exp(tinc/htau); 
aamGate[vrt] = minf; 
aahGate[vrt] = hinf; 
}  
 
 
 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
 number FARADAY = m_pVMDisc->F; 
 number dt = newTime - m_pVMDisc->time(); 
number v = vrt_values[CableEquation<TDomain>::_v_]; 
number ca = vrt_values[CableEquation<TDomain>::_ca_]; 

 
double m = aamGate[vrt]; 
double h = aahGate[vrt]; 

 
 

 
 
v = v+vshift; 
double v_ = v+vshift; 
v_ = v; 
double 	minf = 1.0 / ( 1 + exp(-(v_+v12m)/vwm) ); 
double 	hinf = 1.0 / ( 1 + exp((v_+v12h)/vwh) ); 
double 	mtau = ( am + 1.0 / ( exp((v_+vm1)/wm1) + exp(-(v_+vm2)/wm2) ) ) ; 
double 	htau = ( ah + 1.0 / ( exp((v_+vh1)/wh1) + exp(-(v_+vh2)/wh2) ) ) ; 
double         tinc = -dt ; 
double         mexp = 1 - exp(tinc/mtau); 
double         hexp = 1 - exp(tinc/htau); 
        m = m + mexp*(minf-m); 
        h = h + hexp*(hinf-h); 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number ca = vrt_values[m_pVMDisc->_ca_]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 

 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
number eca; 
if (m_pVMDisc->rev_pot_ca() == 0) 
{ 
	  eca = helpV*(log(m_pVMDisc->ca_out()/ca)); 
} 
else 
{ 
	  eca = m_pVMDisc->rev_pot_ca(); 
} 
 
 
number gca = gbar*m*m*h; 

 
 
outCurrentValues.push_back( gca * (v - eca)); 
} 
 
 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class CaT_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class CaT_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class CaT_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


