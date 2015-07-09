#include "passive_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void passive_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double passive_converted_standard_UG<TDomain>::getg() 
{ 
return g; 
} 
template<typename TDomain> 
double passive_converted_standard_UG<TDomain>::gete() 
{ 
return e; 
} 
template<typename TDomain> 
void passive_converted_standard_UG<TDomain>::setg(double val) 
{ 
g = val; 
} 
template<typename TDomain> 
void passive_converted_standard_UG<TDomain>::sete(double val) 
{ 
e = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void passive_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
} 
 
 
 
template<typename TDomain> 
std::vector<number> passive_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
void passive_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
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
number v = vrt_values[VMDisc<TDomain>::_v_]; 

 
}  
 
 
 
template<typename TDomain> 
void passive_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
 number FARADAY = m_pVMDisc->F; 
 number dt = newTime - m_pVMDisc->time(); 
number v = vrt_values[VMDisc<TDomain>::_v_]; 

 

 
 

 
 
 
 
 
} 
 
 
 
template<typename TDomain> 
void passive_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 

 
 
outCurrentValues.push_back( g*(v - e)); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class passive_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class passive_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class passive_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


