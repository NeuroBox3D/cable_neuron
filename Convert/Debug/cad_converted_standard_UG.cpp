#include "cad_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double cad_converted_standard_UG<TDomain>::getdepth() 
{ 
return depth; 
} 
template<typename TDomain> 
double cad_converted_standard_UG<TDomain>::gettaur() 
{ 
return taur; 
} 
template<typename TDomain> 
double cad_converted_standard_UG<TDomain>::getcainf() 
{ 
return cainf; 
} 
template<typename TDomain> 
double cad_converted_standard_UG<TDomain>::getcai() 
{ 
return cai; 
} 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::setdepth(double val) 
{ 
depth = val; 
} 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::settaur(double val) 
{ 
taur = val; 
} 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::setcainf(double val) 
{ 
cainf = val; 
} 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::setcai(double val) 
{ 
cai = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->caSGate)) 
UG_THROW("Attachment necessary (caSGate) for cad_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->caSGate); 
this->aacaSGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->caSGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
{ 
//get celsius and time
number celsius = m_pVMDisc->celsius; 
number dt = m_pVMDisc->time(); 
number ica = m_pVMDisc->get_flux_ca(); 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t size_l = 0; size_l< vl.size(); size_l++) 
{ 
	 Vertex* vrt = vl[size_l]; 
 
 
number v = u(m_pVMDisc->_v_, size_l); 
number ca = u(m_pVMDisc->_ca_, size_l); 

 
aacaSGate[vrt] =  cainf; 
cai =  ca; 
}  
}  
 
 
 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
{ 
number celsius = m_pVMDisc->celsius; 
 number FARADAY = m_F; 
 number ica = m_pVMDisc->get_flux_ca(); 
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
number ca = u(m_pVMDisc->_ca_, size_l); 

 
double caS = aacaSGate[vrt]; 

 
 
double 	drive_channel =  - (10000) * ica / (2 * FARADAY * depth); 
; 
 
if (drive_channel <= 0.)
{ 
 drive_channel = 0. ; 
} 
caS +=  drive_channel + (cainf-caS)/taur*dt; 
; 
 
	cai = ca; 
; 
 

 
 
aacaSGate[vrt] = caS; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number ica = m_pVMDisc->get_flux_ca(); 
number caS = aacaSGate[ver]; 
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
 
 

 
 
outCurrentValues.push_back(   caS); 
outCurrentValues.push_back(   caS ); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class cad_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class cad_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class cad_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
