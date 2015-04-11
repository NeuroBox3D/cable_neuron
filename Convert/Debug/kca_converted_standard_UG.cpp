#include "kca_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::getgbar() 
{ 
return gbar; 
} 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::getcai() 
{ 
return cai; 
} 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::getcaix() 
{ 
return caix; 
} 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::getRa() 
{ 
return Ra; 
} 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::getRb() 
{ 
return Rb; 
} 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::gettemp() 
{ 
return temp; 
} 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::getq10() 
{ 
return q10; 
} 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::getvmin() 
{ 
return vmin; 
} 
template<typename TDomain> 
double kca_converted_standard_UG<TDomain>::getvmax() 
{ 
return vmax; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::setgbar(double val) 
{ 
gbar = val; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::setcai(double val) 
{ 
cai = val; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::setcaix(double val) 
{ 
caix = val; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::setRa(double val) 
{ 
Ra = val; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::setRb(double val) 
{ 
Rb = val; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::settemp(double val) 
{ 
temp = val; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::setq10(double val) 
{ 
q10 = val; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::setvmin(double val) 
{ 
vmin = val; 
} 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::setvmax(double val) 
{ 
vmax = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for kca_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->nGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
{ 
//get celsius and time
number celsius = m_pVMDisc->celsius; 
number dt = m_pVMDisc->time(); 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t size_l = 0; size_l< vl.size(); size_l++) 
{ 
	 Vertex* vrt = vl[size_l]; 
 
 
number v = u(m_pVMDisc->_v_, size_l); 
number k = u(m_pVMDisc->_k_, size_l); 
number ca = u(m_pVMDisc->_ca_, size_l); 

 
double cai = ca; 
 
double a=Ra* pow(cai , caix); 
double         b = Rb; 
tadj= pow(q10 , ((celsius-temp)/10)); 
double         ntau = 1/tadj/(a+b); 
double 	ninf = a/(a+b); 
;//        tinc = -dt * tadj
;//        nexp = 1 - exp(tinc/ntau)
aanGate[vrt] = ninf; 
}  
}  
 
 
 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
{ 
number celsius = m_pVMDisc->celsius; 
 number FARADAY = m_F; 
 // make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t size_l = 0; size_l< vl.size(); size_l++) 
{ 
	 Vertex* vrt = vl[size_l]; 
 
 
number dt = newTime - m_pVMDisc->m_aaTime[vrt]; 
number v = u(m_pVMDisc->_v_, size_l); 
number k = u(m_pVMDisc->_k_, size_l); 
number ca = u(m_pVMDisc->_ca_, size_l); 

 
double n = aanGate[vrt]; 

 
 
double cai = ca; 
 
double a=Ra* pow(cai , caix); 
double         b = Rb; 
tadj= pow(q10 , ((celsius-temp)/10)); 
double         ntau = 1/tadj/(a+b); 
double 	ninf = a/(a+b); 
//        tinc = -dt * tadj
//        nexp = 1 - exp(tinc/ntau)
        n  +=   (ninf-n)/ntau*dt; 
; 
 

 
 
aanGate[vrt] = n; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void kca_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number n = aanGate[ver]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number ek; 
if (m_pVMDisc->get_ek() == 0) 
{ 
	  ek = helpV*(log(m_pVMDisc->k_out/k)); 
} 
else 
{ 
	  ek = m_pVMDisc->get_ek(); 
} 
 
 

 
 
number gk = tadj*gbar*n; 
outCurrentValues.push_back( (1e-4) * gk * (v - ek)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class kca_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class kca_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class kca_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
