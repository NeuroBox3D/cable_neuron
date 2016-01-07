/*
 * leakage.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann, mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__Na_K_Pump_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__Na_K_Pump_H__


#include "channel_interface.h"


namespace ug {
namespace cable {


template <typename TDomain>
class Na_K_Pump
	: public ICableMembraneTransport<TDomain>
{
	public:
		using ICableMembraneTransport<TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*)
		Na_K_Pump(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&)
		Na_K_Pump(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~Na_K_Pump() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		void set_K_K(number K);
		void set_K_Na(number Na);

		void set_IMAX_P(number IMAX);

		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		number K_K;
		number K_Na;
		number IMAX_P;				// mol*s^-1

};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__LEAKAGE_H__
