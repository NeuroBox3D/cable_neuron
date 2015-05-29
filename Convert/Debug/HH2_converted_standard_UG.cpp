#include "HH2_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::vtrap(double x, double y) 
{ 
	if (fabs(x/y) < 1e-6) {
		return  y*(1 - x/y/2); 
	}else{
		return  x/(Exp(x/y)-1); 
	}
}
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::Exp(double x) 
{ 
	if (x < -100) {
		return  0; 
	}else{
		return  exp(x); 
	}
} 

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getgnabar() 
{ 
return gnabar; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getgkbar() 
{ 
return gkbar; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getena() 
{ 
return ena; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getek() 
{ 
return ek; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getcelsius() 
{ 
return celsius; 
} 
template<typename TDomain> 
double HH2_converted_standard_UG<TDomain>::getvtraub() 
{ 
return vtraub; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setgnabar(double val) 
{ 
gnabar = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setgkbar(double val) 
{ 
gkbar = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setena(double val) 
{ 
ena = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setek(double val) 
{ 
ek = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setcelsius(double val) 
{ 
celsius = val; 
} 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::setvtraub(double val) 
{ 
vtraub = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for HH2_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for HH2_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
if (spGrid->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for HH2_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->nGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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
number na = u(m_pVMDisc->_na_, size_l); 
number k = u(m_pVMDisc->_k_, size_l); 

 
tadj =  pow(3.0 , ((celsius-36)/10)); 
aamGate[vrt] = 0; 
aahGate[vrt] = 0; 
aanGate[vrt] = 0; 
}  
}  
 
 
 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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
number na = u(m_pVMDisc->_na_, size_l); 
number k = u(m_pVMDisc->_k_, size_l); 

 
double m = aamGate[vrt]; 
double h = aahGate[vrt]; 
double n = aanGate[vrt]; 

 
 

 
 
v = v; 
double 	v2 = v - vtraub ;// convert to traub convention; 
//       a = 0.32 * (13-v2) / ( Exp((13-v2)/4) - 1)
double 	a = 0.32 * vtrap(13-v2, 4); 
//       b = 0.28 * (v2-40) / ( Exp((v2-40)/5) - 1)
double 	b = 0.28 * vtrap(v2-40, 5); 
double 	tau_m = 1 / (a + b) / tadj; 
double 	m_inf = a / (a + b); 
	a = 0.128 * Exp((17-v2)/18); 
	b = 4 / ( 1 + Exp((40-v2)/5) ); 
double 	tau_h = 1 / (a + b) / tadj; 
double 	h_inf = a / (a + b); 
//       a = 0.032 * (15-v2) / ( Exp((15-v2)/5) - 1)
	a = 0.032 * vtrap(15-v2, 5); 
	b = 0.5 * Exp((10-v2)/40); 
double 	tau_n = 1 / (a + b) / tadj; 
double 	n_inf = a / (a + b); 
double 	m_exp = 1 - Exp(-dt/tau_m); 
double 	h_exp = 1 - Exp(-dt/tau_h); 
double 	n_exp = 1 - Exp(-dt/tau_n); 
	m = m + m_exp * (m_inf - m); 
	h = h + h_exp * (h_inf - h); 
	n = n + n_exp * (n_inf - n); 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
aanGate[vrt] = n; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void HH2_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number n = aanGate[ver]; 
number na = vrt_values[VMDisc<TDomain>::_na_]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
 
 

 
 
outCurrentValues.push_back( gnabar * m*m*m*h * (v - ena) +  gkbar * n*n*n*n * (v - ek)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class HH2_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class HH2_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class HH2_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


