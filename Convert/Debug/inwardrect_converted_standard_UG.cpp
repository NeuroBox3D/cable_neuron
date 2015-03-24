#include "inwardrect_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void inwardrect_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
 // creating Method for attachments 
template<typename TDomain> 
void inwardrect_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for inwardrect_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->nGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void inwardrect_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
{ 
//get celsius 
number celsius = m_pVMDisc->celsius; 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t l = 0; l< vl.size(); l++) 
{ 
	 Vertex* vrt = vl[l]; 
 
 
number v = u(m_pVMDisc->_v_, l); 
number k = u(m_pVMDisc->_k_, l); 

 
double          tinc; 
double v = v; 
                      //Call once from HOC to initialize inf at resting v.
double         b = Ra * (v - tha) / (1 - exp(-(v - tha)/qa)); 
double         a = -Rb * (v - tha) / (1 - exp((v - tha)/qa)); 
double         ntau = 1/(a+b); 
double 	ninf = a*ntau; 
                      //Call once from HOC to initialize inf at resting v.
tadj= pow(q10 , ((celsius-temp)/10)); 
        tinc = -dt * tadj; 
double         nexp = 1 - exp(tinc/ntau); 
aanGate[vrt] = ninf; 
}  
}  
 
 
 
template<typename TDomain> 
void inwardrect_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
{ 
number celsius = m_pVMDisc->celsius; 
 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t l = 0; l< vl.size(); l++) 
{ 
	 Vertex* vrt = vl[l]; 
 
 
number dt = newTime - m_pVMDisc->m_aaTime[vrt]; 
number v = u(m_pVMDisc->_v_, l); 
number k = u(m_pVMDisc->_k_, l); 

 
double n = aanGate[vrt]; 

 
 

 
 
v = v; 
v = v; 
                      //Call once from HOC to initialize inf at resting v.
double         b = Ra * (v - tha) / (1 - exp(-(v - tha)/qa)); 
double         a = -Rb * (v - tha) / (1 - exp((v - tha)/qa)); 
double         ntau = 1/(a+b); 
double 	ninf = a*ntau; 
double tadj= pow(q10 , ((celsius-temp)/10)); 
double         tinc = -dt * tadj; 
double         nexp = 1 - exp(tinc/ntau); 
        n = n + nexp*(ninf-n); 
aanGate[vrt] = n; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void inwardrect_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number n = aanGate[ver]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number ek = helpV*(log(m_pVMDisc->k_out/k)); 
 
 

 
 
outCurrentValues.push_back( (1e-4) * gk * (v - ek)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class inwardrect_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class inwardrect_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class inwardrect_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
