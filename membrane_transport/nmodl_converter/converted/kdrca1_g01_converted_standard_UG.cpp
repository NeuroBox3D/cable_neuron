#include "kdrca1_g01_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable_neuron { 
 
 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::alpn(double v) 
{ 
  return  exp(1.e-3*zetan*(v-vhalfn)*9.648e4(degC/mV)/(8.315*(273.16+celsius))) ; 
}
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::betn(double v) 
{ 
  return  exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4(degC/mV)/(8.315*(273.16+celsius))) ; 
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::ce_obj_available()  
{  
	init_attachments();  
 	F = m_pCE->F; 
 R = m_pCE->R; 
 K = m_pCE->temperature(); 
 celsius = m_pCE->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getek() 
{ 
return ek; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getgkdrbar() 
{ 
return gkdrbar; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getikmax() 
{ 
return ikmax; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getvhalfn() 
{ 
return vhalfn; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::geta0n() 
{ 
return a0n; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getzetan() 
{ 
return zetan; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getgmn() 
{ 
return gmn; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getnmax() 
{ 
return nmax; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getq10() 
{ 
return q10; 
} 
template<typename TDomain> 
double kdrca1_g01_converted_standard_UG<TDomain>::getnscale() 
{ 
return nscale; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setek(double val) 
{ 
ek = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setgkdrbar(double val) 
{ 
gkdrbar = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setikmax(double val) 
{ 
ikmax = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setvhalfn(double val) 
{ 
vhalfn = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::seta0n(double val) 
{ 
a0n = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setzetan(double val) 
{ 
zetan = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setgmn(double val) 
{ 
gmn = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setnmax(double val) 
{ 
nmax = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setq10(double val) 
{ 
q10 = val; 
} 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::setnscale(double val) 
{ 
nscale = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->SGate)) 
UG_THROW("Attachment necessary (SGate) for kdrca1_g01_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->SGate); 
this->aaSGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->SGate); 
 
if (spGrid->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for kdrca1_g01_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->nGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> kdrca1_g01_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 if (m_log_SGate || m_log_nGate )
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
	 	 if (m_log_nGate == true) 
	 	 	 GatingAccesors.push_back(this->aanGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void kdrca1_g01_converted_standard_UG<TDomain>::set_log_SGate(bool bLogSGate) { m_log_SGate = bLogSGate; }
template<typename TDomain> void kdrca1_g01_converted_standard_UG<TDomain>::set_log_nGate(bool bLognGate) { m_log_nGate = bLognGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
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
number k = vrt_values[CableEquation<TDomain>::_k_]; 

 
double          a,qt; 
qt= pow(q10 , ((celsius-24)/10(degC))); 
        a = alpn(v); 
double         ninf = 1/(1+a); 
double         taun = betn(v)/(qt*a0n*(1+a)); 
	if (taun<nmax) {taun=nmax;};; 
        taun=taun/nscale; 
double         facn = (1 - exp(-dt/taun)); 
aanGate[vrt] =ninf; 
double gkdr =  gkdrbar*n; 
double ik =  gkdr*(v-ek); 
}  
 
 
 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
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
number k = vrt_values[CableEquation<TDomain>::_k_]; 

 
double S = aaSGate[vrt]; 
double n = aanGate[vrt]; 

 
 
double          a,qt; 
qt= pow(q10 , ((celsius-24)/10(degC))); 
        a = alpn(v); 
double         ninf = 1/(1+a); 
double         taun = betn(v)/(qt*a0n*(1+a)); 
	if (taun<nmax) {taun=nmax;};; 
        taun=taun/nscale; 
double         facn = (1 - exp(-dt/taun)); 
        n  +=  (ninf-n)/taun*dt; 
; 
 

 
 
aaSGate[vrt] = S; 
aanGate[vrt] = n; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number S = aaSGate[ver]; 
number n = aanGate[ver]; 
number k = vrt_values[m_pCE->_k_]; 
number v =  vrt_values[m_pCE->_v_]; 
 
 
number t = m_pCE->time(); 
 
 

 
 
const number helpV = 1e3*(m_pCE->R*m_pCE->temperature())/m_pCE->F; 
 
 
	SOLVE states METHOD cnexp; 
	gkdr = gkdrbar*n; 
	ik = gkdr*(v-ek); 
double i; 
 
if (ik>ikmax)
{ 
 ik=ikmax ; 
} 
}
outCurrentValues.push_back(0); 
} 
 
 
template<typename TDomain> 
void kdrca1_g01_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class kdrca1_g01_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class kdrca1_g01_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class kdrca1_g01_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable_neuron
} // namespace ug


