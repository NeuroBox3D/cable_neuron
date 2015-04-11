#include "h_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::alpt(double v) 
{ 
  return  exp(0.0378*zetat*(v-vhalft)) ; 
}
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::bett(double v) 
{ 
  return  exp(0.0378*zetat*gmt*(v-vhalft)) ; 
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getehd() 
{ 
return ehd; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getghdbar() 
{ 
return ghdbar; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getvhalfl() 
{ 
return vhalfl; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getkl() 
{ 
return kl; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getvhalft() 
{ 
return vhalft; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::geta0t() 
{ 
return a0t; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getzetat() 
{ 
return zetat; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getgmt() 
{ 
return gmt; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getq10() 
{ 
return q10; 
} 
template<typename TDomain> 
double h_converted_standard_UG<TDomain>::getqtl() 
{ 
return qtl; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setehd(double val) 
{ 
ehd = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setghdbar(double val) 
{ 
ghdbar = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setvhalfl(double val) 
{ 
vhalfl = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setkl(double val) 
{ 
kl = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setvhalft(double val) 
{ 
vhalft = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::seta0t(double val) 
{ 
a0t = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setzetat(double val) 
{ 
zetat = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setgmt(double val) 
{ 
gmt = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setq10(double val) 
{ 
q10 = val; 
} 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::setqtl(double val) 
{ 
qtl = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->lGate)) 
UG_THROW("Attachment necessary (lGate) for h_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->lGate); 
this->aalGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->lGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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

 
double          a,qt; 
qt= pow(q10 , ((celsius-33)/10)); 
        a = alpt(v); 
        linf = 1/(1 + exp(-(v-vhalfl)/kl)); 
;//       linf = 1/(1+ alpl(v))
double         taul = bett(v)/(qtl*qt*a0t*(1+a)); 
aalGate[vrt] =linf; 
}  
}  
 
 
 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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

 
double l = aalGate[vrt]; 

 
 
double          a,qt; 
qt= pow(q10 , ((celsius-33)/10)); 
        a = alpt(v); 
        linf = 1/(1 + exp(-(v-vhalfl)/kl)); 
//       linf = 1/(1+ alpl(v))
double         taul = bett(v)/(qtl*qt*a0t*(1+a)); 
        l  +=   (linf - l)/taul*dt; 
; 
 

 
 
aalGate[vrt] = l; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void h_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number l = aalGate[ver]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
 
 
number ghd = ghdbar*l; 

 
 
outCurrentValues.push_back( ghd*(v-ehd)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class h_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class h_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class h_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
