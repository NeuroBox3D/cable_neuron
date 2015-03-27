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
 
 
 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getgbar() 
{ 
return gbar; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getvshift() 
{ 
return vshift; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::gettha() 
{ 
return tha; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getqa() 
{ 
return qa; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getRa() 
{ 
return Ra; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getRb() 
{ 
return Rb; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getthi1() 
{ 
return thi1; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getthi2() 
{ 
return thi2; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getqi() 
{ 
return qi; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getthinf() 
{ 
return thinf; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getqinf() 
{ 
return qinf; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getRg() 
{ 
return Rg; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getRd() 
{ 
return Rd; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::gettemp() 
{ 
return temp; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getq10() 
{ 
return q10; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getvmin() 
{ 
return vmin; 
} 
template<typename TDomain> 
double na_converted_standard_UG<TDomain>::getvmax() 
{ 
return vmax; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setgbar(double val) 
{ 
gbar = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setvshift(double val) 
{ 
vshift = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::settha(double val) 
{ 
tha = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setqa(double val) 
{ 
qa = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setRa(double val) 
{ 
Ra = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setRb(double val) 
{ 
Rb = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setthi1(double val) 
{ 
thi1 = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setthi2(double val) 
{ 
thi2 = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setqi(double val) 
{ 
qi = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setthinf(double val) 
{ 
thinf = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setqinf(double val) 
{ 
qinf = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setRg(double val) 
{ 
Rg = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setRd(double val) 
{ 
Rd = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::settemp(double val) 
{ 
temp = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setq10(double val) 
{ 
q10 = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setvmin(double val) 
{ 
vmin = val; 
} 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::setvmax(double val) 
{ 
vmax = val; 
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

 
double           a, b; 
double vm = v; 
	a = trap0(vm,tha,Ra,qa); 
	b = trap0(-vm,-tha,Rb,qa); 
tadj= pow(q10 , ((celsius-temp)/10)); 
double 	mtau = 1/tadj/(a+b); 
double 	minf = a/(a+b); 
		;//"h" inactivation 
	a = trap0(vm,thi1,Rd,qi); 
	b = trap0(-vm,-thi2,Rg,qi); 
double 	htau = 1/tadj/(a+b); 
double 	hinf = 1/(1+exp((vm-thinf)/qinf)); 
;//        tinc = -dt * tadj
;//        mexp = 1 - exp(tinc/mtau)
;//        hexp = 1 - exp(tinc/htau)
aamGate[vrt] = minf; 
aahGate[vrt] = hinf; 
}  
}  
 
 
 
template<typename TDomain> 
void na_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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
 
 

 
 
number gna = tadj*gbar*m*m*m*h; 
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
  
  
