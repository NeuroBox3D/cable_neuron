#include "pump_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
 	F = m_pVMDisc->F; 
 R = m_pVMDisc->R; 
 K = m_pVMDisc->temperature(); 
 celsius = m_pVMDisc->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double pump_converted_standard_UG<TDomain>::getnai() 
{ 
return nai; 
} 
template<typename TDomain> 
double pump_converted_standard_UG<TDomain>::getipumpmax() 
{ 
return ipumpmax; 
} 
template<typename TDomain> 
double pump_converted_standard_UG<TDomain>::getkm() 
{ 
return km; 
} 
template<typename TDomain> 
double pump_converted_standard_UG<TDomain>::getn() 
{ 
return n; 
} 
template<typename TDomain> 
double pump_converted_standard_UG<TDomain>::getnainit() 
{ 
return nainit; 
} 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::setnai(double val) 
{ 
nai = val; 
} 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::setipumpmax(double val) 
{ 
ipumpmax = val; 
} 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::setkm(double val) 
{ 
km = val; 
} 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::setn(double val) 
{ 
n = val; 
} 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::setnainit(double val) 
{ 
nainit = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
} 
 
 
 
template<typename TDomain> 
std::vector<number> pump_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
{ 
	 //var for output 
	 std::vector<number> GatingAccesors; 
 
	 typedef ug::MathVector<TDomain::dim> position_type; 
 
	 position_type coord; 
 
	 if (coord.size()==1) 
	 	 coord[0]=x; 
	 if (coord.size()==2) 
	 { 
	 	 coord[0] = x;
	 	 coord[1] = y;
	 } 
	 if (coord.size()==3) 
	 { 
	 	 coord[0] = x;
	 	 coord[1] = y;
	 	 coord[2] = z;
	 } 
	 //accesors 
	 typedef Attachment<position_type> position_attachment_type; 
	 typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accesor_type; 
 
	 // Definitions for Iteration over all Elements 
	 typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; 
	 SubsetGroup ssGrp; 
	 try { ssGrp = SubsetGroup(m_pVMDisc->approx_space()->domain()->subset_handler(), this->m_vSubset);} 
	 UG_CATCH_THROW("Subset group creation failed."); 
 
	 itType iter; 
	 number bestDistSq, distSq; 
	 Vertex* bestVrt; 
 
	 // Iterate only if there is one Gtting needed 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
 // Init Method for using gatings 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
//get celsius and time
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
number dt = m_pVMDisc->time(); 
// make preparing vor getting values of every edge 
number v = vrt_values[CableEquation<TDomain>::_v_]; 
number na = vrt_values[CableEquation<TDomain>::_na_]; 
number k = vrt_values[CableEquation<TDomain>::_k_]; 

 
}  
 
 
 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
 number FARADAY = m_pVMDisc->F; 
 number dt = newTime - m_pVMDisc->time(); 
number v = vrt_values[CableEquation<TDomain>::_v_]; 
number na = vrt_values[CableEquation<TDomain>::_na_]; 
number k = vrt_values[CableEquation<TDomain>::_k_]; 

 

 
 

 
 
 
 
 
} 
 
 
 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number na = vrt_values[m_pVMDisc->_na_]; 
number k = vrt_values[m_pVMDisc->_k_]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
number nai = na;

 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 
number inapump = ipumpmax*(1/(1 + pow(km/nai,n))); 

 
 
outCurrentValues.push_back( 3.0*inapump +  -2.0*inapump); 
} 
 
 
template<typename TDomain> 
void pump_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class pump_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class pump_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class pump_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


