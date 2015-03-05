#include "hh_converted_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
template<typename TDomain, typename TAlgebra> 
double hh_converted_UG<TDomain, TAlgebra>::vtrap(double x, double y)
{ 
        if (fabs(x/y) < 1e-6) {
                return  y*(1 - x/y/2); 
        }else{
                return  x/(exp(x/y) - 1); 
        }
}

 
 // Init Method for using gatings 
template<typename TDomain, typename TAlgebra> 
void hh_converted_UG<TDomain, TAlgebra>::init(number time, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) 
{ 
//get celsius 
number celsius = m_vmDisc->celsius; 
// attach attachments 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->mGate); 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->hGate); 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->nGate); 
 
// creates multiindeces 
std::vector<DoFIndex> multIndna; 
std::vector<DoFIndex> multIndk; 
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
size_t fct_id_na = spGridFct->fct_id_by_name("na"); 
size_t fct_id_k = spGridFct->fct_id_by_name("k"); 

 
spGridFct->dof_indices(*iter, fct_id_na, multIndna); 
spGridFct->dof_indices(*iter, fct_id_k, multIndk); 
spGridFct->dof_indices(*iter, vm, multIndv);  
 
 UG_ASSERT(multIndv.size() == 1, "multi-index has size !=1 "); 
number v = DoFRef(*spGridFct, multIndv[0]);  
number na = DoFRef(*spGridFct, multIndna[0]);  
number k = DoFRef(*spGridFct, multIndk[0]);  
  
 
 
double           alpha, beta, sum, q10; 
                      //Call once from HOC to initialize inf at resting v.
q10= pow(3 , ((celsius-6.3)/10)); 
                //"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10); 
        beta =  4 * exp(-(v+65)/18); 
        sum = alpha + beta; 
	mtau = 1/(q10*sum); 
        minf = alpha/sum; 
                //"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(q10*sum); 
        hinf = alpha/sum; 
                //"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) ; 
        beta = .125*exp(-(v+65)/80); 
	sum = alpha + beta; 
        ntau = 1/(q10*sum); 
        ninf = alpha/sum; 
aamGate[*iter] = minf; 
aahGate[*iter] = hinf; 
aanGate[*iter] = ninf; 
}  
}  
}  
 
 
 
template<typename TDomain, typename TAlgebra> 
void hh_converted_UG<TDomain, TAlgebra>::update_gating(number dt, ConstSmartPtr<typename TAlgebra::vector_type> uOld) 
{ 
number celsius = m_vmDisc->celsius; 
 
// creates multiindeces 
std::vector<DoFIndex> multIndna; 
std::vector<DoFIndex> multIndk; 
std::vector<DoFIndex> multIndv; 
 
typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; 
SubsetGroup ssGrp; 
try { ssGrp = SubsetGroup(m_vmDisc->approx_space()->subset_handler(), this->m_vSubset);} 
UG_CATCH_THROW("Subset group creation failed."); 
 
for (std::size_t si = 0; si < ssGrp.size(); si++) 
{ 
 
ConstSmartPtr<DoFDistribution> dd = m_vmDisc->approx_space()->dof_distribution(GridLevel::TOP); 
itType iterBegin = m_vmDisc->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);  
itType iterEnd = m_vmDisc->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
for (itType iter = iterBegin; iter != iterEnd; ++iter) 
{ 
 
Vertex* vrt = *iter; 
// needed konzentration has to be set 
const size_t vm = dd->fct_id_by_name("v");  
const size_t fct_id_na = dd->fct_id_by_name("na"); 
const size_t fct_id_k = dd->fct_id_by_name("k"); 

 
dd->dof_indices(*iter, fct_id_na, multIndna); 
dd->dof_indices(*iter, fct_id_k, multIndk); 
dd->dof_indices(*iter, vm, multIndv);  
 
 UG_ASSERT(multIndv.size() == 1, "multi-index has size !=1 "); 
number v = DoFRef(*uOld, multIndv[0]);  
number na = DoFRef(*uOld, multIndna[0]);  
number k = DoFRef(*uOld, multIndk[0]);  
  
 
 
double m = aamGate[*iter]; 
double h = aahGate[*iter]; 
double n = aanGate[*iter]; 

 
 
double           alpha, beta, sum, q10; 
                      //Call once from HOC to initialize inf at resting v.
q10= pow(3 , ((celsius-6.3)/10)); 
                //"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10); 
        beta =  4 * exp(-(v+65)/18); 
        sum = alpha + beta; 
	mtau = 1/(q10*sum); 
        minf = alpha/sum; 
                //"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(q10*sum); 
        hinf = alpha/sum; 
                //"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) ; 
        beta = .125*exp(-(v+65)/80); 
	sum = alpha + beta; 
        ntau = 1/(q10*sum); 
        ninf = alpha/sum; 
        m +=  (minf-m)/mtau*dt; 
        h += (hinf-h)/htau*dt; 
        n += (ninf-n)/ntau*dt; 

 
 
aamGate[*iter] = m; 
aahGate[*iter] = h; 
aanGate[*iter] = n; 
 
 
 
} 
} 
} 
 
 
 
template<typename TDomain, typename TAlgebra> 
void hh_converted_UG<TDomain, TAlgebra>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number n = aanGate[ver]; 
number na = vrt_values[VMDisc<TDomain, TAlgebra>::_na_]; 
number k = vrt_values[VMDisc<TDomain, TAlgebra>::_k_]; 
number v =  vrt_values[VMDisc<TDomain, TAlgebra>::_v_]; 
 
 
const number helpV = (m_R*m_T)/m_F; 
number ena = helpV*(log(m_vmDisc->na_out/na)); 
number ek = helpV*(log(m_vmDisc->k_out/k)); 
 
 
number gna = gnabar*m*m*m*h; 
number gk = gkbar*n*n*n*n; 

 
 
outCurrentValues.push_back( gna*(v - ena) +  gk*(v - ek)       +  gl*(v - el)); 
outCurrentValues.push_back( gna*(v - ena)/m_F ); 
outCurrentValues.push_back( gk*(v - ek)      /m_F ); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
#ifdef UG_CPU_1 
template class hh_converted_UG<Domain1d, CPUAlgebra>; 
#endif 
#ifdef UG_CPU_2 
template class hh_converted_UG<Domain1d, CPUBlockAlgebra<2> >; 
#endif 
#ifdef UG_CPU_3 
template class hh_converted_UG<Domain1d, CPUBlockAlgebra<3> >; 
#endif 
#ifdef UG_CPU_4 
template class hh_converted_UG<Domain1d, CPUBlockAlgebra<4> >; 
#endif 
#ifdef UG_CPU_VAR 
template class hh_converted_UG<Domain1d, CPUVariableBlockAlgebra >; 
#endif 
#endif 
 
 
#ifdef UG_DIM_2 
#ifdef UG_CPU_1 
template class hh_converted_UG<Domain2d, CPUAlgebra>; 
#endif 
#ifdef UG_CPU_2 
template class hh_converted_UG<Domain2d, CPUBlockAlgebra<2> >; 
#endif 
#ifdef UG_CPU_3 
template class hh_converted_UG<Domain2d, CPUBlockAlgebra<3> >; 
#endif 
#ifdef UG_CPU_4 
template class hh_converted_UG<Domain2d, CPUBlockAlgebra<4> >; 
#endif 
#ifdef UG_CPU_VAR 
template class hh_converted_UG<Domain2d, CPUVariableBlockAlgebra >; 
#endif 
#endif 
 
 
#ifdef UG_DIM_3 
#ifdef UG_CPU_1 
template class hh_converted_UG<Domain3d, CPUAlgebra>; 
#endif 
#ifdef UG_CPU_2 
template class hh_converted_UG<Domain3d, CPUBlockAlgebra<2> >; 
#endif 
#ifdef UG_CPU_3 
template class hh_converted_UG<Domain3d, CPUBlockAlgebra<3> >; 
#endif 
#ifdef UG_CPU_4 
template class hh_converted_UG<Domain3d, CPUBlockAlgebra<4> >; 
#endif 
#ifdef UG_CPU_VAR 
template class hh_converted_UG<Domain3d, CPUVariableBlockAlgebra >; 
#endif 
#endif 
 
 
}  
  
  
