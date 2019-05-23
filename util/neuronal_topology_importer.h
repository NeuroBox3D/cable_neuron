/*!
 * \file neuronal_topology_importer.h
 * \brief provides geometry importing and grid generating routines
 *        (swc, hoc, ngx, txt, xml to ugx)
 *
 *  Created on:      2014-10-14
 *  Original author: Stephan Grein
 *
 *  Restructured and augmented 2019/04 by mbreit.
 */

#ifndef UG__PLUGINS__CABLE_NEURON__UTIL__NEURONAL_TOPOLOGY_IMPORTER_H
#define UG__PLUGINS__CABLE_NEURON__UTIL__NEURONAL_TOPOLOGY_IMPORTER_H

#include "common/math/ugmath_types.h"  // for vector3
#include "common/util/smart_pointer.h"  // for SmartPtr
#include "lib_grid/common_attachments.h"  // for ANumber
#include "lib_grid/grid/grid.h"  // for AttachmentAccessors
#include "lib_grid/multi_grid.h"  // for MultiGrid
#include "lib_grid/tools/subset_handler_multi_grid.h"  // for MGSubsetHandler

#include "../synapse_handling/synapse_info_io_traits.h"
#include "../synapse_handling/synapses/base_synapse.h"  // for IBaseSynapse

#include <map>
#include <string>
#include <vector>

#include <boost/property_tree/ptree.hpp>


/// debug ids
extern ug::DebugID NETIGeometry;
extern ug::DebugID NETIGrid;



namespace ug {
namespace cable_neuron {
namespace neuronal_topology_importer {


struct PointInfo
{
	int index;            // for hoc and ngx files: index of the section the point belongs to),
	                      // for swc, txt, nml: type of point (soma, axon, dend, etc.);
	                      // used to assign subsets
	vector3 coordinates;  // point 3d coordinates
	number diameter;      // diameter
	size_t connectsTo;    // connected to point with this index
	bool isDuplicate;     // whether this point is a connection and exists twice
	                      // (used only in txt import)
};


struct SectionInfo
{
	size_t index;  // section index
	size_t start;  // index of the first point of this section
	size_t end;    // index of the last point in this section
};


struct SynapseInfo
{
	SynapseType type;          // Exp2Syn or AlphaSyn
	std::string secName_pre;   // pre-synaptic section name
	std::string secName_post;  // post-synaptic section name
	number secLocCoord_pre;    // section-local coordinates of pre-synapse
	number secLocCoord_post;   // section-local coordinates of post-synapse
	number tau1;               // time constant [s]
	number tau2;               // time constant [s] (only exp2)
	number e;                  // reversal potential [V]
	number threshold;          // activation threshold potential [V] (only exp2)
	number onset;              // onset of activation [s] (only alpha)
	number delay;              // activation delay [s] (not used)
	number gMax;               // peak conductance [S]
};


/// NeuronalTopologyImporter
/*!
 * \brief imports swc and hoc geometries and generates ugx grids from them
 */
class NeuronalTopologyImporter
{
	public:
		NeuronalTopologyImporter();

		virtual ~NeuronalTopologyImporter();

		/// set all scaling factors
		void set_scaling_factors(number length, number time, number potential, number conductance);

		/// set length scaling factor
		void set_length_scaling(number lengthScale);

		/// set time scaling factor
		void set_time_scaling(number timeScale);

		/// set potential scaling factor
		void set_potential_scaling(number potentialScale);

		/// set conductance scaling factor
		void set_conductance_scaling(number conductanceScale);


		/*!
		 * \brief add and remove joining criteria
		 */
		void add_joining_criterion(const char* str);
		void rm_joining_criterion(const char* str);

		/*!
		 * \brief imports and generates the grid
		 * \param[in] file path to geometry file
		 */
		bool import_and_generate_grid(const std::string& file);

		/*!
		 * \brief imports and generates the grid
		 * \param[in] file path to geometry file
		 * \param[in] type geometry file type
		 */
		bool import_and_generate_grid(const std::string& file, const std::string& type);


	protected:
		/*!
		 * \brief imports the geometry
		 * \param[in] file path to geometry file
		 * \param[in] type geometry file type
		 * \return \c bool
		 */
		bool import_geometry(const std::string& file);
		bool import_geometry_from_SWC();
		bool import_geometry_from_HOC();
		bool import_geometry_from_NGX();
		bool import_geometry_from_TXT();
		bool import_geometry_from_NML();

		void read_identifier_file_txt(const std::string& fid);

		void extract_points_and_sections_from_hoc(const std::string& text);
		void extract_points_and_sections_from_ngx(const boost::property_tree::ptree& tree);
		void extract_points_and_sections_from_txt(const std::string& fsec);

		void extract_connections_from_hoc(const std::string& text);
		void extract_connections_from_ngx(const boost::property_tree::ptree& tree);
		void extract_connections_from_txt(const std::string& fconn);

		void extract_synapses_from_hoc(const std::string& text);
		void extract_synapses_from_ngx(const boost::property_tree::ptree& tree);
		void extract_synapses_from_txt(const std::string& fsynapses);

		/*!
		 * \brief generates the grid, using the supplied grid as initial grid, out of the imported geometry
		 * \param[in] grid supplied initial grid
		 */
		bool generate_grid();

		void subset_processing();
		void subset_processing_SWC();
		void subset_processing_HOC();
		void subset_processing_NGX();
		void subset_processing_TXT();
		void subset_processing_NML();


	private:
		static const number REMOVE_DOUBLES_THRESHOLD;

		std::string m_filename;
		std::string m_extension;

		/// scale length (typically from um -> m)
		number m_scaleLength;

		/// scale time (typically from ms -> s)
		number m_scaleTime;

		/// scale potential (typically from mV -> V)
		number m_scalePotential;

		/// scale conductance (typically from uS -> S)
		number m_scaleConductance;

		std::vector<std::string> m_joiningCriteria;

		std::vector<PointInfo> m_vPtInfo;
		std::map<std::string, SectionInfo> m_mSecInfo;
		std::vector<SynapseInfo> m_vSynInfo;

		SmartPtr<MultiGrid> m_grid;
		SmartPtr<MGSubsetHandler> m_sh;

		ANumber m_aDiameterAttachment;
		Grid::VertexAttachmentAccessor<ANumber> m_aaDiameter;

		typedef Attachment<std::vector<synapse_handler::IBaseSynapse*> > ASynapse;
		ASynapse m_aSynapse;
		Grid::EdgeAttachmentAccessor<ASynapse> m_aaSynapse;

		// txt-specific variables
		int m_networkID;  // 0 for neocortex, 1 for hippocampus
		size_t m_numCellTypes;
		bool m_bWithCellType;
};


}  // namespace neuronal_topology_importer
}  // namespace cable_neuron
}  // namespace ug

#endif  // UG__PLUGINS__CABLE_NEURON__UTIL__NEURONAL_TOPOLOGY_IMPORTER_H
