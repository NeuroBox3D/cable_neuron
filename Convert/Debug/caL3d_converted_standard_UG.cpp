#include "caL3d_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::ghk(double v, double  ci, double  co) 
{ 
	double z; 

	z = (0.001)*2*F*v/(R*(celsius+273.15)); 
	return  (.001)*2*F*(ci*efun(-z) - co*efun(z)); 
}
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::efun(double z) 
{ 
	if (fabs(z) < 1e-4) {
		return  1 - z/2; 
	}else{
		return  z/(exp(z) - 1); 
	}
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getp() 
{ 
return p; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getth() 
{ 
return th; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getq() 
{ 
return q; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getRa() 
{ 
return Ra; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getRb() 
{ 
return Rb; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::gettemp() 
{ 
return temp; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getq10() 
{ 
return q10; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setp(double val) 
{ 
p = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setth(double val) 
{ 
th = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setq(double val) 
{ 
q = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setRa(double val) 
{ 
Ra = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setRb(double val) 
{ 
Rb = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::settemp(double val) 
{ 
temp = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setq10(double val) 
{ 
q10 = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->CGate)) 
UG_THROW("Attachment necessary (CGate) for caL3d_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->CGate); 
this->aaCGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->CGate); 
 
if (spGrid->has_vertex_attachment(this->OGate)) 
UG_THROW("Attachment necessary (OGate) for caL3d_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->OGate); 
this->aaOGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->OGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
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
number ca = u(m_pVMDisc->_ca_, size_l); 

 
aaCGate[vrt] = 1 ; 
}  
}  
 
 
 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
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
number ca = u(m_pVMDisc->_ca_, size_l); 

 
double C = aaCGate[vrt]; 
double O = aaOGate[vrt]; 

 
 

 
 
tadj= pow(q10 , ((celsius-temp)/10(degC))); 
	a = Ra / (1 + exp(-(v-th)/q)) * tadj; 
	b = Rb / (1 + exp((v-th)/q)) * tadj; 
 
 
 
C+=(-C*a+O*b)*dt; 
O+=(C*a+-O*b)*dt; 
 
 
 
aaCGate[vrt] = C; 
aaOGate[vrt] = O; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number C = aaCGate[ver]; 
number O = aaOGate[ver]; 
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number cai =  ca; 
 
 

 
 
number rates(v); 
number cao = m_pVMDisc->ca_out; 

 
 
outCurrentValues.push_back( O * p * ghk(v,cai,cao)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class caL3d_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class caL3d_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class caL3d_converted_standard_UG<Domain3d>; 
#endif 
 
 
}  
  
  
