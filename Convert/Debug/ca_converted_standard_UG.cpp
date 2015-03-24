#include "ca_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
template<typename TDomain> 
double ca_converted_standard_UG<TDomain>::efun(double z) 
{ 
	if (fabs(z) < 1e-4) {
		return  1 - z/2; 
	}else{
		return  z/(exp(z) - 1); 
	}
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void ca_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
 // creating Method for attachments 
template<typename TDomain> 
void ca_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for ca_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for ca_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void ca_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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

 
double           a, b; 
double vm = v; 
tadj= pow(q10 , ((celsius-temp)/10)); 
	a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1); 
	b = 0.94*exp((-75-vm)/17); 
double 	mtau = 1/tadj/(a+b); 
double 	minf = a/(a+b); 
		//"h" inactivation 
	a = 0.000457*exp((-13-vm)/50); 
	b = 0.0065/(exp((-vm-15)/28) + 1); 
double 	htau = 1/tadj/(a+b); 
double 	hinf = a/(a+b); 
//        tinc = -dt * tadj
//        mexp = 1 - exp(tinc/mtau)
//        hexp = 1 - exp(tinc/htau)
aamGate[vrt] = minf; 
aahGate[vrt] = hinf; 
}  
}  
 
 
 
template<typename TDomain> 
void ca_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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

 
 
double           a, b; 
double vm = v; 
tadj= pow(q10 , ((celsius-temp)/10)); 
	a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1); 
	b = 0.94*exp((-75-vm)/17); 
double 	mtau = 1/tadj/(a+b); 
double 	minf = a/(a+b); 
		//"h" inactivation 
	a = 0.000457*exp((-13-vm)/50); 
	b = 0.0065/(exp((-vm-15)/28) + 1); 
double 	htau = 1/tadj/(a+b); 
double 	hinf = a/(a+b); 
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
void ca_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number eca = helpV*(log(m_pVMDisc->ca_out/ca)); 
 
 

 
 
outCurrentValues.push_back( (1e-4) * gca * (v - eca)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class ca_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class ca_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class ca_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
