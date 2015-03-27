#include "kir_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::taumkir(double v)
{  
return v; 
} 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getgkbar() 
{ 
return gkbar; 
} 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getmvhalf() 
{ 
return mvhalf; 
} 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getmslope() 
{ 
return mslope; 
} 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getmshift() 
{ 
return mshift; 
} 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getqfact() 
{ 
return qfact; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setgkbar(double val) 
{ 
gkbar = val; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setmvhalf(double val) 
{ 
mvhalf = val; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setmslope(double val) 
{ 
mslope = val; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setmshift(double val) 
{ 
mshift = val; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setqfact(double val) 
{ 
qfact = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for kir_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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

 
double 			minf = 1  /  ( 1 + exp( (v - mvhalf + mshift) / mslope) ); 
aamGate[vrt] = minf; 
}  
}  
 
 
 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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

 
double m = aamGate[vrt]; 

 
 
double 			minf = 1  /  ( 1 + exp( (v - mvhalf + mshift) / mslope) ); 
        m  +=  (minf - m) / ( taumkir(v)/qfact )*dt; 

 
 
aamGate[vrt] = m; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number m = aamGate[ver]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number ek = helpV*(log(m_pVMDisc->k_out/k)); 
 
 
number gk = gkbar * m; 

 
 
outCurrentValues.push_back( gk * ( v - ek )); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class kir_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class kir_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class kir_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
