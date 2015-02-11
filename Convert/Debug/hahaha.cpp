#include " hahaha.h"	
#include " lib_grid/lg_base.h" 
#include " lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include " lib_disc/function_spaces/grid_function.h" 
#include " lib_disc/function_spaces/local_transfer_interface.h" 
namespace ug { 
 
 
 
{ 

 
 // Init Method for using gatings 
template<typename TDomain, typename TAlgebra> 
void hahaha<TDomain, TAlgebra>::init(number time, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) 
{ 
// attach attachments 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->vGate)) 
UG_THROW("Attachment necessary (vGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->vGate); 
this->aavGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->vGate); 
 
// creates multiindeces 
std::vector<DoFIndex> multIndv; 
 
typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;  
SubsetGroup ssGrp;  
try { ssGrp = SubsetGroup(spGridFct->domain()->subset_handler(), this->m_vSubset);}  
UG_CATCH_THROW("Subset group creation failed."); 
  
for (std::size_t si = 0; si < ssGrp.size(); si++)  
{  
itType iterBegin = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);  
itType iterEnd = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
  
for (itType iter = iterBegin; iter != iterEnd; ++iter)  
{  
 
size_t vm = spGridFct->fct_id_by_name("v");  

 
spGridFct->dof_indices(*iter, vm, multIndv);  
 
 for (size_t i=0; i<multIndv.size(); i++)  
{  
aavGate[*iter] = DoFRef(*spGridFct, multIndv[i]);  
}  
 
 
}  
}  
}  
 
 
 
template<typename TDomain, typename TAlgebra> 
void hahaha<TDomain, TAlgebra>::update_gating(number newTime, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) 
{ 
typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; 
SubsetGroup ssGrp; 
try { ssGrp = SubsetGroup(spGridFct->domain()->subset_handler(), this->m_vSubset);} 
UG_CATCH_THROW("Subset group creation failed."); 
 
for (std::size_t si = 0; si < ssGrp.size(); si++) 
{ 
 
itType iterBegin = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);  
itType iterEnd = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
for (itType iter = iterBegin; iter != iterEnd; ++iter) 
{ 
 
Vertex* vrt = *iter; 
// needed konzentration has to be set 
size_t vm = spGridFct->fct_id_by_name("v");  

 
spGridFct->dof_indices(*iter, vm, multIndv);  
 
 for (size_t i=0; i<multIndv.size(); i++)  
{  
aavGate[*iter] = DoFRef(*spGridFct, multIndv[i]);  
}  
 
 

 
 

 
 
 
 
 
} 
} 
} 
 
 
 
template<typename TDomain, typename TAlgebra> 
void hahaha<TDomain, TAlgebra>::ionic_current(Vertex* ver, std::vector<number>& outCurrentValues) 
{ 
 
double v = m_aav[ver]; 
 
 
const number helpV = (m_R*m_T)/m_F; 
 
 

 
 
outCurrentValues.push_back( +  g*(v - e)/m_F ) 
  
}  
  
  
