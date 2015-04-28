#include "release_BMK_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void release_BMK_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double release_BMK_converted_standard_UG<TDomain>::getdel() 
{ 
return del; 
} 
template<typename TDomain> 
double release_BMK_converted_standard_UG<TDomain>::getdur() 
{ 
return dur; 
} 
template<typename TDomain> 
double release_BMK_converted_standard_UG<TDomain>::getamp() 
{ 
return amp; 
} 
template<typename TDomain> 
void release_BMK_converted_standard_UG<TDomain>::setdel(double val) 
{ 
del = val; 
} 
template<typename TDomain> 
void release_BMK_converted_standard_UG<TDomain>::setdur(double val) 
{ 
dur = val; 
} 
template<typename TDomain> 
void release_BMK_converted_standard_UG<TDomain>::setamp(double val) 
{ 
amp = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void release_BMK_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void release_BMK_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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

 
double T =  0; 
}  
}  
 
 
 
template<typename TDomain> 
void release_BMK_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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
void release_BMK_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
 
 
; 
; 
; 
double T; 
 
if (t < del + dur && t > del)
{ 
		T = amp; 
} 
else 
{ 
		T = 0; 
}
outCurrentValues.push_back(0);
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class release_BMK_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class release_BMK_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class release_BMK_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
