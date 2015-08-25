/*
 * ion_leakage.h
 *
 *  Created on: 24.08.2015
 *      Author: pgottmann, mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__ION_LEAKAGE_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__ION_LEAKAGE_H__


#include "channel_interface.h"


namespace ug {
namespace cable {


template <typename TDomain>
class IonLeakage
	: public IChannel<TDomain>
{
	public:
		using IChannel<TDomain>::m_pVMDisc;

		/// @copydoc IChannel<TDomain>::IChannel(const char*)
		IonLeakage(const char* functions, const char* subsets);

		/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&)
		IonLeakage(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~IonLeakage() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// set leaking substance
		void set_leaking_quantity(const std::string& lq);

		/// set permeability
		void set_perm(number flux_at_rest, number conc_in_rest, number conc_out_rest, number vm_rest);

		// inherited from IChannel
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void vm_disc_available();

	private:
		virtual void specify_write_function_indices();

	private:
		// membrane conductivities
		number m_perm;		// C / (m^2 * mV * ms)

		// leaking function name
		std::string m_leaking_fct;

		// leaking function index in VMDisc
		size_t m_lfInd;

		number m_flux_at_rest;
		number m_conc_in_rest;
		number m_conc_out_rest;
		number m_vm_rest;
};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__ION_LEAKAGE_H__
