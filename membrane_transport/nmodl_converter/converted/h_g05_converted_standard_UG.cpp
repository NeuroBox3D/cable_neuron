#include "h_g05_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable_neuron { 
 
 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::alpha(double v) 
{ 
	return  a0*exp((0.001)*(-ab)*zeta*(v-vhalf)*farad/(gas*(273.16+celsius))) ; 
}
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::beta(double v) 
{ 
	return  a0*exp((0.001)*(1-ab)*zeta*(v-vhalf)*farad/(gas*(273.16+celsius))); 
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::ce_obj_available()  
{  
	init_attachments();  
 	F = m_pCE->F; 
 R = m_pCE->R; 
 K = m_pCE->temperature(); 
 celsius = m_pCE->temperature_celsius(); 
}  
 

 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::getgbar() 
{ 
return gbar; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::geterevh() 
{ 
return erevh; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::getvhalf() 
{ 
return vhalf; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::geta0() 
{ 
return a0; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::getzeta() 
{ 
return zeta; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::getab() 
{ 
return ab; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::getqten() 
{ 
return qten; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::gettemp() 
{ 
return temp; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::getgas() 
{ 
return gas; 
} 
template<typename TDomain> 
double h_g05_converted_standard_UG<TDomain>::getfarad() 
{ 
return farad; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::setgbar(double val) 
{ 
gbar = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::seterevh(double val) 
{ 
erevh = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::setvhalf(double val) 
{ 
vhalf = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::seta0(double val) 
{ 
a0 = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::setzeta(double val) 
{ 
zeta = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::setab(double val) 
{ 
ab = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::setqten(double val) 
{ 
qten = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::settemp(double val) 
{ 
temp = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::setgas(double val) 
{ 
gas = val; 
} 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::setfarad(double val) 
{ 
farad = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->SGate)) 
UG_THROW("Attachment necessary (SGate) for h_g05_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->SGate); 
this->aaSGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->SGate); 
 
if (spGrid->has_vertex_attachment(this->hhGate)) 
UG_THROW("Attachment necessary (hhGate) for h_g05_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hhGate); 
this->aahhGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hhGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> h_g05_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 try { ssGrp = SubsetGroup(m_pCE->approx_space()->domain()->subset_handler(), this->m_vSubset);} 
	 UG_CATCH_THROW("Subset group creation failed."); 
 
	 itType iter; 
	 number bestDistSq, distSq; 
	 Vertex* bestVrt; 
 
	 // Iterate only if there is one Gtting needed 
	 if (m_log_SGate || m_log_hhGate )
	 { 
	 	 // iterating over all elements 
	 	 for (size_t si=0; si < ssGrp.size(); si++) 
	 	 { 
	 	 	 itType iterBegin = m_pCE->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]); 
	 	 	 itType iterEnd = m_pCE->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
	 	 	 const position_accesor_type& aaPos = m_pCE->approx_space()->domain()->position_accessor(); 
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
	 	 if (m_log_hhGate == true) 
	 	 	 GatingAccesors.push_back(this->aahhGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void h_g05_converted_standard_UG<TDomain>::set_log_SGate(bool bLogSGate) { m_log_SGate = bLogSGate; }
template<typename TDomain> void h_g05_converted_standard_UG<TDomain>::set_log_hhGate(bool bLoghhGate) { m_log_hhGate = bLoghhGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
//get celsius and time
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number celsius = m_pCE->temperature_celsius(); 
number dt = m_pCE->time(); 
// make preparing vor getting values of every edge 
number v = vrt_values[CableEquation<TDomain>::_v_]; 

 
double      a, b, q10; 
q10= pow(qten , ((celsius-temp)/10));
    a = q10*alpha(v); 
    b = q10*beta(v); 
    inf = a/(a+b); 
double     tau = 1/(a+b); 
    if (tau<2) {tau=2;};; 
aahhGate[vrt] =inf; 
}  
 
 
 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number celsius = m_pCE->temperature_celsius(); 
 number FARADAY = m_pCE->F; 
 number dt = newTime - m_pCE->time(); 
number v = vrt_values[CableEquation<TDomain>::_v_]; 

 
double S = aaSGate[vrt]; 
double hh = aahhGate[vrt]; 

 
 
double      a, b, q10; 
q10= pow(qten , ((celsius-temp)/10));
    a = q10*alpha(v); 
    b = q10*beta(v); 
    inf = a/(a+b); 
double     tau = 1/(a+b); 
    if (tau<2) {tau=2;};; 
    hh  +=  (inf-hh)/tau*dt; 
; 
 

 
 
aaSGate[vrt] = S; 
aahhGate[vrt] = hh; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number S = aaSGate[ver]; 
number hh = aahhGate[ver]; 
number v =  vrt_values[m_pCE->_v_]; 
 
 
number t = m_pCE->time(); 
 
 

 
 
const number helpV = 1e3*(m_pCE->R*m_pCE->temperature())/m_pCE->F; 
 
 

 
 
number Ih = (0.0001)*gbar*hh*(v-erevh); 
outCurrentValues.push_back(0 + Ih); 
} 
 
 
template<typename TDomain> 
void h_g05_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class h_g05_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class h_g05_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class h_g05_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable_neuron
} // namespace ug


