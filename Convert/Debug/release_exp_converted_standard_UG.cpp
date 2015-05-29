#include "release_exp_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double release_exp_converted_standard_UG<TDomain>::gettau1() 
{ 
return tau1; 
} 
template<typename TDomain> 
double release_exp_converted_standard_UG<TDomain>::gettau2() 
{ 
return tau2; 
} 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::settau1(double val) 
{ 
tau1 = val; 
} 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::settau2(double val) 
{ 
tau2 = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->AGate)) 
UG_THROW("Attachment necessary (AGate) for release_exp_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->AGate); 
this->aaAGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->AGate); 
 
if (spGrid->has_vertex_attachment(this->BGate)) 
UG_THROW("Attachment necessary (BGate) for release_exp_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->BGate); 
this->aaBGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->BGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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

 
total =  0; 
tau1 =  .9999*tau2; 
aaAGate[vrt] = 0; 
aaBGate[vrt] = 0; 
double tp =  (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1); 
double factor =  -exp(-tp/tau1) + exp(-tp/tau2); 
factor =  1/factor;
}  
}  
 
 
 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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

 
double A = aaAGate[vrt]; 
double B = aaBGate[vrt]; 

 
 
    B  +=  -B/tau2*dt; 
; 
 

 
 
aaAGate[vrt] = A; 
aaBGate[vrt] = B; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number A = aaAGate[ver]; 
number B = aaBGate[ver]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
 
 

 
 
number T = B - A; 
outCurrentValues.push_back(0); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class release_exp_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class release_exp_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class release_exp_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


