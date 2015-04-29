#include "myseclamp_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double myseclamp_converted_standard_UG<TDomain>::getrs() 
{ 
return rs; 
} 
template<typename TDomain> 
double myseclamp_converted_standard_UG<TDomain>::getdur1() 
{ 
return dur1; 
} 
template<typename TDomain> 
double myseclamp_converted_standard_UG<TDomain>::getdur2() 
{ 
return dur2; 
} 
template<typename TDomain> 
double myseclamp_converted_standard_UG<TDomain>::getdur3() 
{ 
return dur3; 
} 
template<typename TDomain> 
double myseclamp_converted_standard_UG<TDomain>::getdur4() 
{ 
return dur4; 
} 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::setrs(double val) 
{ 
rs = val; 
} 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::setdur1(double val) 
{ 
dur1 = val; 
} 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::setdur2(double val) 
{ 
dur2 = val; 
} 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::setdur3(double val) 
{ 
dur3 = val; 
} 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::setdur4(double val) 
{ 
dur4 = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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

 
double tc2 =  dur1 + dur2; 
double tc3 =  tc2 + dur3; 
double tc4 =  tc3 + dur4; 
double on =  0; 
}  
}  
 
 
 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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

 

 
 

 
 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void myseclamp_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
 
 

 
 
number vstim(); 
outCurrentValues.push_back(0); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class myseclamp_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class myseclamp_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class myseclamp_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
