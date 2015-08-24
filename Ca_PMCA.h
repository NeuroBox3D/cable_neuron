/*
 * leakage.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann, mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__Ca_PMCA_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__Ca_PMCA_H__


#include "channel_interface.h"


namespace ug {
namespace cable {


template <typename TDomain>
class Ca_PMCA
	: public IChannel<TDomain>
{
	public:
		using IChannel<TDomain>::m_pVMDisc;

		/// @copydoc IChannel<TDomain>::IChannel(const char*)
		Ca_PMCA(const char* functions, const char* subsets);

		/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&)
		Ca_PMCA(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~Ca_PMCA() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		void set_KD_P(number KD);

		void set_IMAX_P(number IMAX);

		// inherited from IChannel
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void vm_disc_available();
		virtual std::vector<number> state_values(number x, number y, number z);
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		number KD_P;					// mol*m^-3 (Elwess et al.)
		//const number KD_P = 3.4e-04;		// mol*m^-3 (Graupner)
		number IMAX_P;				// mol*s^-1


		ADouble gatingFactorGate;
		Grid::AttachmentAccessor<Vertex, ADouble> aagatingFactorGate;
};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__LEAKAGE_H__
