#include "HH2_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::vtrap(double x, double y) 
{ 
double vtrap; 
	if (fabs(x/y) < 1e-6) {
		return  y*(1 - x/y/2); 
	}else{
		return  x/(Exp(x/y)-1); 
	}
}
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::Exp(double x) 
{ 
double Exp; 
	if (x < -100) {
		return  0; 
	}else{
		return  exp(x); 
	}
} 

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
 	F = m_pVMDisc->F; 
 R = m_pVMDisc->R; 
 K = m_pVMDisc->temperature(); 
 celsius = m_pVMDisc->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getgnabar() 
{ 
return gnabar; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getgkbar() 
{ 
return gkbar; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getena() 
{ 
return ena; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getek() 
{ 
return ek; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getvtraub() 
{ 
return vtraub; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setgnabar(double val) 
{ 
gnabar = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setgkbar(double val) 
{ 
gkbar = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setena(double val) 
{ 
ena = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setek(double val) 
{ 
ek = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setvtraub(double val) 
{ 
vtraub = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->SGate)) 
UG_THROW("Attachment necessary (SGate) for HH2_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->SGate); 
this->aaSGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->SGate); 
 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for HH2_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for HH2_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
if (spGrid->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for HH2_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->nGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> HH2_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 if (m_log_SGate || m_log_mGate || m_log_hGate || m_log_nGate )
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
	 	 if (m_log_SGate == true) 
	 	 	 GatingAccesors.push_back(this->aaSGate[bestVrt]); 
	 	 if (m_log_mGate == true) 
	 	 	 GatingAccesors.push_back(this->aamGate[bestVrt]); 
	 	 if (m_log_hGate == true) 
	 	 	 GatingAccesors.push_back(this->aahGate[bestVrt]); 
	 	 if (m_log_nGate == true) 
	 	 	 GatingAccesors.push_back(this->aanGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void HH2_converted_standard_UG<TDomain>::set_log_SGate(bool bLogSGate) { m_log_SGate = bLogSGate; }
template<typename TDomain> void HH2_converted_standard_UG<TDomain>::set_log_mGate(bool bLogmGate) { m_log_mGate = bLogmGate; }
template<typename TDomain> void HH2_converted_standard_UG<TDomain>::set_log_hGate(bool bLoghGate) { m_log_hGate = bLoghGate; }
template<typename TDomain> void HH2_converted_standard_UG<TDomain>::set_log_nGate(bool bLognGate) { m_log_nGate = bLognGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
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
number na = vrt_values[CableEquation<TDomain>::_na_]; 
number k = vrt_values[CableEquation<TDomain>::_k_]; 

 
}  
 
 
 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
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
number na = vrt_values[CableEquation<TDomain>::_na_]; 
number k = vrt_values[CableEquation<TDomain>::_k_]; 

 
double S = aaSGate[vrt]; 
double m = aamGate[vrt]; 
double h = aahGate[vrt]; 
double n = aanGate[vrt]; 

 
 

 
 
v = v; 
double 	v2 = v - vtraub ;// convert to traub convention; 
//--//       a = 0.32 * (13-v2) / ( Exp((13-v2)/4) - 1); 
double 	a = 0.32 * vtrap(13-v2, 4); 
//--//       b = 0.28 * (v2-40) / ( Exp((v2-40)/5) - 1); 
double 	b = 0.28 * vtrap(v2-40, 5); 
double 	tau_m = 1 / (a + b) / tadj; 
double 	m_inf = a / (a + b); 
	a = 0.128 * Exp((17-v2)/18); 
	b = 4 / ( 1 + Exp((40-v2)/5) ); 
double 	tau_h = 1 / (a + b) / tadj; 
double 	h_inf = a / (a + b); 
//--//       a = 0.032 * (15-v2) / ( Exp((15-v2)/5) - 1); 
	a = 0.032 * vtrap(15-v2, 5); 
	b = 0.5 * Exp((10-v2)/40); 
double 	tau_n = 1 / (a + b) / tadj; 
double 	n_inf = a / (a + b); 
double 	m_exp = 1 - Exp(-dt/tau_m); 
double 	h_exp = 1 - Exp(-dt/tau_h); 
double 	n_exp = 1 - Exp(-dt/tau_n); 
	m = m + m_exp * (m_inf - m); 
	h = h + h_exp * (h_inf - h); 
	n = n + n_exp * (n_inf - n); 
aaSGate[vrt] = S; 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
aanGate[vrt] = n; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number S = aaSGate[ver]; 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number n = aanGate[ver]; 
number na = vrt_values[m_pVMDisc->_na_]; 
number k = vrt_values[m_pVMDisc->_k_]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 

 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 

 
 
outCurrentValues.push_back( gnabar * m*m*m*h * (v - ena) +  gkbar * n*n*n*n * (v - ek)); 
} 
 
 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class HH2_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class HH2_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class HH2_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


