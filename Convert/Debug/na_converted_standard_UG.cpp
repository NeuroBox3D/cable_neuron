#include "na_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::trap0(double v, double th, double a, double q) 
{ 
	if (fabs((v-th)/q) > 1e-6) {
	        return  a * (v - th) / (1 - exp(-(v - th)/q)); 
	} else {
	        return  a * q; 
 	}
}	

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
 // creating Method for attachments 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for na_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for na_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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
number na = u(m_pVMDisc->_na_, l); 

 
double           a, b; 
double vm = v; 
	a = trap0(vm,tha,Ra,qa); 
	b = trap0(-vm,-tha,Rb,qa); 
tadj= pow(q10 , ((celsius-temp)/10)); 
double 	mtau = 1/tadj/(a+b); 
double 	minf = a/(a+b); 
		//"h" inactivation 
	a = trap0(vm,thi1,Rd,qi); 
	b = trap0(-vm,-thi2,Rg,qi); 
double 	htau = 1/tadj/(a+b); 
double 	hinf = 1/(1+exp((vm-thinf)/qinf)); 
//        tinc = -dt * tadj
//        mexp = 1 - exp(tinc/mtau)
//        hexp = 1 - exp(tinc/htau)
aamGate[vrt] = minf; 
aahGate[vrt] = hinf; 
}  
}  
 
 
 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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
number na = u(m_pVMDisc->_na_, l); 

 
double m = aamGate[vrt]; 
double h = aahGate[vrt]; 

 
 
double           a, b; 
double vm = v; 
	a = trap0(vm,tha,Ra,qa); 
	b = trap0(-vm,-tha,Rb,qa); 
tadj= pow(q10 , ((celsius-temp)/10)); 
double 	mtau = 1/tadj/(a+b); 
double 	minf = a/(a+b); 
		//"h" inactivation 
	a = trap0(vm,thi1,Rd,qi); 
	b = trap0(-vm,-thi2,Rg,qi); 
double 	htau = 1/tadj/(a+b); 
double 	hinf = 1/(1+exp((vm-thinf)/qinf)); 
//        tinc = -dt * tadj
//        mexp = 1 - exp(tinc/mtau)
//        hexp = 1 - exp(tinc/htau)
        m  +=   (minf-m)/mtau*dt; 
        h  +=   (hinf-h)/htau*dt; 

 
 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number na = vrt_values[VMDisc<TDomain>::_na_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number ena = helpV*(log(m_pVMDisc->na_out/na)); 
 
 

 
 
outCurrentValues.push_back( (1e-4) * gna * (v - ena)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class na_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class na_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class na_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
