#ifndef hahaha_H_
#define hahaha_H_
#include "lib_grid/lg_base.h" 
#include "lib_grid/grid/grid_base_objects.h" 

#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h" 
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h" 
#include "lib_disc/spatial_disc/disc_util/geom_provider.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include "lib_disc/common/local_algebra.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/interpolate.h" 

#include "bridge/bridge.h" 
#include "bridge/util.h" 
#include "bridge/util_domain_algebra_dependent.h" 
#include "bridge/util_domain_dependent.h" 

#include "common/util/smart_pointer.h" 
#include "common/util/vector_util.h" 

#include <vector> 
#include <stdio.h> 
#include "bindings/lua/lua_user_data.h" 
namespace ug { 
template <typename TDomain, typename TAlgebra> 
class hahaha
    : public IChannel<TDomain, TAlgebra> 
{ 
    public: 
 
 

 
 
const static 	g = .001	; 
const static 	e = -70	; 
std::vector<number> m_diff; 



/// constructor 
 
ChannelHHNernst(const char* functions, const char* subsets) 
: IChannel<TDomain, TAlgebra>(functions, subsets), 
m_bNonRegularGrid(false), m_R(8314), m_T(279.45), m_F(96485.0) 
{ 
register_all_funcs(m_bNonRegularGrid); 
}; 
 
/// destructor 
 
virtual ~ChannelHHNernst() {}; 
// inherited from IChannel 
 
virtual const std::vector<number> get_diff(); 
virtual void set_diff(const std::vector<number>& diff); 
virtual void init(number time, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct); 
virtual void update_gating(number newTime, SmartPtr<GridFunction<TDomain, TAlgebra> > sgGridFct); 
virtual void ionic_current(Vertex* v, std::vector<number>& outCurrentValues); 
/// assembles the local right hand side 
 
template<typename TElem, typename TFVGeom> 
void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]); 

 
protected: 
bool m_bNonRegularGrid; 
 
void register_all_funcs(bool bHang); 
template <typename TElem, typename TFVGeom> 
void register_func(); 
 
 
private: 
 
ADouble v; 
Grid::AttachmentAccessor<Vertex, ADouble> aav; 
 
//nernst const values 
number m_R 
number m_T 
number m_F 
 
/// Base type 
typedef IChannel<TDomain, TAlgebra> base_type; 
///	Own type 
typedef ChannelHH<TDomain, TAlgebra> this_type; 
/// GridFunction type 
typedef GridFunction<TDomain, TAlgebra> TGridFunction; 
 
}; 
 
} // namespace ug 
 
 
#endif // hahaha_H_
