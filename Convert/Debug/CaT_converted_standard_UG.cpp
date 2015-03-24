#include "CaT_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
 // creating Method for attachments 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for CaT_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for CaT_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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
number ca = u(m_pVMDisc->_ca_, l); 

 
double          tinc; 
double           a, b; 
double v_ = v; 
double 	minf = 1.0 / ( 1 + exp(-(v_+v12m)/vwm) ); 
double 	hinf = 1.0 / ( 1 + exp((v_+v12h)/vwh) ); 
double 	mtau = ( am + 1.0 / ( exp((v_+vm1)/wm1) + exp(-(v_+vm2)/wm2) ) ) ; 
double 	htau = ( ah + 1.0 / ( exp((v_+vh1)/wh1) + exp(-(v_+vh2)/wh2) ) ) ; 
        tinc = -dt ; 
double         mexp = 1 - exp(tinc/mtau); 
double         hexp = 1 - exp(tinc/htau); 
aamGate[vrt] = minf; 
aahGate[vrt] = hinf; 
}  
}  
 
 
 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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
number ca = u(m_pVMDisc->_ca_, l); 

 
double m = aamGate[vrt]; 
double h = aahGate[vrt]; 

 
 

 
 
v = v+vshift; 
double v_ = v+vshift; 
double v_ = v; 
double 	minf = 1.0 / ( 1 + exp(-(v_+v12m)/vwm) ); 
double 	hinf = 1.0 / ( 1 + exp((v_+v12h)/vwh) ); 
double 	mtau = ( am + 1.0 / ( exp((v_+vm1)/wm1) + exp(-(v_+vm2)/wm2) ) ) ; 
double 	htau = ( ah + 1.0 / ( exp((v_+vh1)/wh1) + exp(-(v_+vh2)/wh2) ) ) ; 
double                       ; 
double         tinc = -dt ; 
double         mexp = 1 - exp(tinc/mtau); 
double         hexp = 1 - exp(tinc/htau); 
        m = m + mexp*(minf-m); 
        h = h + hexp*(hinf-h); 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void CaT_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number eca = helpV*(log(m_pVMDisc->ca_out/ca)); 
 
 

 
 
outCurrentValues.push_back( gca * (v - eca)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class CaT_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class CaT_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class CaT_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
