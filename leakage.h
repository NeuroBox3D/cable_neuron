/*
 * leakage.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__LEAKAGE_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__LEAKAGE_H__


#include "channel_interface.h"


namespace ug {
namespace cable {


template <typename TDomain>
class ChannelLeak
	: public IChannel<TDomain>
{
	public:
		using IChannel<TDomain>::m_pVMDisc;

		/// @copydoc IChannel<TDomain>::IChannel(const char*)
		ChannelLeak(const char* functions, const char* subsets);

		/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&)
		ChannelLeak(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~ChannelLeak() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// set leakage conductivity
		void set_cond(number g);

		/// set leakage equilibrium potential
		void set_rev_pot(number e);

		// inherited from IChannel
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void vm_disc_available();
		virtual std::vector<number> state_values(number x, number y, number z);
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		// membrane conductivities
		number m_g;		// C / (m^2 * mV * ms)

		// equilibrium potential
		number m_E;
};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__LEAKAGE_H__
