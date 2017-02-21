/*
 * synapse_distributor.cpp
 *
 *  Created on: Mar 24, 2016
 *      Author: lreinhardt
 */

#include "synapse_distributor.h"
#include "lib_grid/global_attachments.h"
#include "../synapse_handling/synapse_info_io_traits.h"
#include "synapse_dealer.h"
#include "../util/functions.h"  // subset_length

// boost includes
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

SynapseDistributor::SynapseDistributor(const std::string& infile)
{
	m_bDomBased = false;

	m_InputFile = infile;
	m_LastMessage = "";

	// grid setup
	m_spGrid = make_sp(new MultiGrid);
	m_spSubsetHandler = make_sp(new MGSubsetHandler(*m_spGrid));
	if (!LoadGridFromFile(*m_spGrid, *m_spSubsetHandler, infile.c_str()))
		UG_THROW("File '" << infile << "' could not be loaded.");

	// initialize synapse attachment and attach to grid, if not already present
	m_aSSyn = GlobalAttachments::attachment<AVSSynapse>("synapses");
	if (!m_spGrid->has_edge_attachment(m_aSSyn))
		m_spGrid->attach_to_edges(m_aSSyn);

	// initialize attachment accessors
	m_aaSSyn = Grid::EdgeAttachmentAccessor<AVSSynapse>(*m_spGrid, m_aSSyn);
	m_aaPosition = Grid::VertexAttachmentAccessor<APosition>(*m_spGrid, aPosition);

    // get maximal synapse ID in grid
    m_nextID = max_id() + 1;

	m_LastMessage = "SynapseDistributor: SynapseDistributor object successfully instantiated.";
}

SynapseDistributor::SynapseDistributor(SmartPtr<Domain3d> dom)
{
	m_bDomBased = true;

	m_InputFile = "";
	m_LastMessage = "";

	// grid setup
	m_spGrid = dom->grid();
	m_spSubsetHandler = dom->subset_handler();

	// initialize attachment and attach to grid, if not already present
	m_aSSyn = GlobalAttachments::attachment<AVSSynapse>("synapses");
	if (!m_spGrid->has_edge_attachment(m_aSSyn))
		m_spGrid->attach_to_edges(m_aSSyn);

	// initialize attachment accessors
	m_aaSSyn = Grid::EdgeAttachmentAccessor<AVSSynapse>(*m_spGrid, m_aSSyn);
	m_aaPosition = Grid::VertexAttachmentAccessor<APosition>(*m_spGrid, aPosition);

    // get maximal synapse ID in grid
    m_nextID = max_id() + 1;

	m_LastMessage = "SynapseDistributor: SynapseDistributor object successfully instantiated.";
}


SynapseDistributor::~SynapseDistributor()
{}



void SynapseDistributor::CopySynapsesFromParentToChild(Edge* parent)
{
	number localBaseCoord = 0.0;
	size_t childEdgeIndex = 0;

//	Iterate over all synapses of the parent edge
	for(size_t i = 0; i < m_aaSSyn[parent].size(); ++i)
	{
	//	Determine index of child edge, which inherits the synapse, and get new local base coords for the child edge
		if(m_aaSSyn[parent][i]->location() < 0.5)
		{
			localBaseCoord = 2*m_aaSSyn[parent][i]->location();
			childEdgeIndex = 0;
		}
		else{
			localBaseCoord = 2*m_aaSSyn[parent][i]->location() - 1;
			childEdgeIndex = 1;
		}

	//	Get correct child edge and the parent's synapse
		Edge* child = m_spGrid->get_child_edge(parent, childEdgeIndex);
		IBaseSynapse* syn = m_aaSSyn[parent][i];

	//	Copy synapse to child
		m_aaSSyn[child].push_back(syn);

	//	Modify the local base coord and associated edge of synapse
		m_aaSSyn[child].back()->set_location(localBaseCoord);
	}
}



void SynapseDistributor::CopySynapsesToAllLevels()
{
//	Iterate over all grid levels
	for(size_t i = 0; i < m_spGrid->num_levels(); ++i)
	{

#ifdef UG_PARALLEL
	// 	copy from vertical masters to vertical slaves
		typedef GridLayoutMap::Types<Edge>::Layout layout_type;
		DistributedGridManager& dgm = *m_spGrid->distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
		pcl::InterfaceCommunicator<layout_type> icom;

		ComPol_CopyAttachment<layout_type, AVSSynapse> compolCopySynapses(*m_spGrid, m_aSSyn);
		icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopySynapses);
		icom.communicate();
#endif

	//	Iterate over all edges of level i
		for(EdgeIterator eIter = m_spGrid->begin<Edge>(i); eIter != m_spGrid->end<Edge>(i); ++eIter)
		{
			Edge* e = *eIter;

			if(m_aaSSyn[e].size() > 0 && m_spGrid->has_children(e))
				CopySynapsesFromParentToChild(e);
		}
	}
}


void SynapseDistributor::remove_synapse(Edge* e)
{
	if (m_aaSSyn[e].size() > 1)
	{
		IBaseSynapse* s_post = m_aaSSyn[e].back(); 	//last pointer
		size_t id = s_post->id();
		delete s_post;								//free memory before pop_back()
		m_aaSSyn[e].pop_back();

		IBaseSynapse* s_pre = m_aaSSyn[e].back();
		size_t id_pre = s_pre->id();
		UG_COND_THROW(id_pre != id, "IDs of consecutive synapses do not match. "
		    "Does this grid contain non-primary synapses?");

		delete s_pre;
		m_aaSSyn[e].pop_back();
	}
}

void SynapseDistributor::remove_all_synapses(Edge* e)
{
	for (size_t i = 0; i < m_aaSSyn[e].size(); ++i)
		delete m_aaSSyn[e][i];

	m_aaSSyn[e].clear();
}

void SynapseDistributor::clear()
{
	for (size_t i = 0; i < m_spGrid->num_levels(); ++i)
	{
        EdgeIterator it = m_spGrid->begin<Edge>(0);
        EdgeIterator it_end = m_spGrid->end<Edge>(0);
		for (; it != it_end; ++it)
		    remove_all_synapses(*it);
	}
}

void SynapseDistributor::clear(int subsetIndex)
{
	for(size_t i = 0; i < m_spGrid->num_levels(); ++i)
	{
		for(EdgeIterator eIter = m_spSubsetHandler->begin<Edge>(subsetIndex, 0); eIter != m_spSubsetHandler->end<Edge>(subsetIndex, 0); eIter++)
		{
			Edge* e = *eIter;
			remove_all_synapses(e);
		}
	}
}


void SynapseDistributor::create_synapse_pair(const std::string& type, IPreSynapse** preOut, IPostSynapse** postOut)
{
    SynapseDealer* dealer = SynapseDealer::instance();

    IBaseSynapse* pre = dealer->deal(std::string("OnsetPreSynapse"));
    pre->set_id(m_nextID);

    IBaseSynapse* post = dealer->deal(type);
    post->set_id(m_nextID);

    *preOut = dynamic_cast<IPreSynapse*>(pre);
    *postOut = dynamic_cast<IPostSynapse*>(post);
    UG_COND_THROW(!*preOut, "Created synapse should be pre-synapse, but is not.");
    UG_COND_THROW(!*postOut, "Created synapse should be post-synapse, but is not.");

    ++m_nextID;
}



void SynapseDistributor::place_synapse_at_coords
(
    const std::vector<number>& coords,
    IPreSynapse* pre,
    IPostSynapse* post
)
{
    // check that grid is set
    UG_COND_THROW(!m_spGrid.valid(), "Grid is not set.");

    // convert coords to vec3
    vector3 c(0);
    for (size_t i = 0; i < coords.size() && i < 3; ++i)
        c[i] = coords[i];

    // find edge with minimal distance to coords
    number minDistSq = std::numeric_limits<number>::infinity();
    Edge* minEdge = NULL;
    number minLocCoord;

    EdgeIterator it = m_spGrid->begin<Edge>(0);
    EdgeIterator it_end = m_spGrid->end<Edge>(0);
    for (; it != it_end; ++it)
    {
        Edge* e = *it;

        // get vertex coords
        vector3 v0 = m_aaPosition[e->vertex(0)];
        vector3 v1 = m_aaPosition[e->vertex(1)];

        number distSq;

        // check whether there is a true distance minimum along the edge
        vector3 edgeVec, c2v0;
        VecSubtract(edgeVec, v1, v0);
        VecSubtract(c2v0, v0, c);
        number lambda = - VecProd(c2v0, edgeVec) / VecNormSquared(edgeVec);
        if (lambda >= 0 && lambda <= 1)
        {
            // determine dist
            VecScaleAdd(c2v0, 1.0, c2v0, lambda, edgeVec); // recycle c2v0
            distSq = VecNormSquared(c2v0);
        }
        // otherwise we have to check both vertices
        else
        {
            VecSubtract(edgeVec, v1, c); // recycle edgeVec
            number dist0Sq = VecNormSquared(c2v0);
            number dist1Sq = VecNormSquared(edgeVec);

            if (dist0Sq < dist1Sq)
            {
                distSq = dist0Sq;
                lambda = 0.0;
            }
            else
            {
                distSq = dist1Sq;
                lambda = 1.0;
            }
        }

        // check if distance is less than current minimum
        if (distSq < minDistSq)
        {
            minDistSq = distSq;
            minEdge = e;
            minLocCoord = lambda;
        }
    }

    UG_COND_THROW(!minEdge, "No location found near given coordinates " << c << " to put synapse to.");

    // create pre- and post-synapses and set id and location
    pre->set_id(m_nextID);
    post->set_id(m_nextID);
    pre->set_location(minLocCoord);
    post->set_location(minLocCoord);

    ++m_nextID;

    m_aaSSyn[minEdge].push_back(pre);
    m_aaSSyn[minEdge].push_back(post);
}


void SynapseDistributor::place_synapse(Edge* e, const std::string& type)
{
    IPreSynapse* pre = NULL;
    IPostSynapse* post = NULL;
    create_synapse_pair(type, &pre, &post);

	number localCoord = static_cast<number>(rand()) / RAND_MAX; //localCoord factor 0 means e[0], 1 means e[1]
	pre->set_location(localCoord);
	post->set_location(localCoord);

	m_aaSSyn[e].push_back(pre);
	m_aaSSyn[e].push_back(post);
}


void SynapseDistributor::place_synapses_uniform(const std::vector<Edge*>& vEdges, size_t numSyn, const std::string& type)
{
//	Determine discrete random distribution of given edges by their lengths
	std::vector<number> probs;
	number totLength = 0.0;

	for(size_t i = 0; i < vEdges.size(); ++i)
	{
		Edge* e = vEdges[i];
		number edgeLength = EdgeLength(e, m_aaPosition);

	//	Push back current edge length as probability (still to be divided by total length of all edges)
		probs.push_back(edgeLength);

		totLength += edgeLength;
	}

//	Scale probabilities
	VecScale(probs, probs, 1./totLength);

//	Define random int generator for random access of edges vector
	boost::random::mt19937 gen;
#ifdef UG_PARALLEL
	gen.seed(pcl::ProcRank());
#endif
	boost::random::discrete_distribution<> dist(probs.begin(), probs.end());
	boost::variate_generator<boost::random::mt19937, boost::random::discrete_distribution<> > randomIndex(gen, dist);

//	Distribute specified number of synapses along the coarse grid edges randomly (uniformly s.t. individual edge lengths)
	for (size_t i = 0; i < numSyn; ++i)
	    place_synapse(vEdges[randomIndex()], type);

//	Handle multigrid transfer
	if (m_bDomBased)
		CopySynapsesToAllLevels();
}


void SynapseDistributor::place_synapses_uniform(size_t numSyn, const std::string& type)
{
//  Save coarse grid edges in a vector
    std::vector<Edge*> vEdges(m_spGrid->begin<Edge>(0), m_spGrid->end<Edge>(0));

//  Call wrapper
    place_synapses_uniform(vEdges, numSyn, type);
}

void SynapseDistributor::place_synapses_uniform(int si, size_t numSyn, const std::string& type)
{
//  Save subset coarse grid edges in a vector
    std::vector<Edge*> vEdges(m_spSubsetHandler->begin<Edge>(si, 0),
                              m_spSubsetHandler->end<Edge>(si, 0));

//  Call wrapper
    place_synapses_uniform(vEdges, numSyn, type);
}

void SynapseDistributor::place_synapses_uniform(const char* subset, size_t numSyn, const std::string& type)
{
    int si = m_spSubsetHandler->get_subset_index(subset);
    place_synapses_uniform(si, numSyn, type);
}

void SynapseDistributor::place_synapses_uniform(int si, number density, const std::string& type)
{
//  Get total length of given subset
    number length = subset_length(si, m_spSubsetHandler);

//  Get number of synapses to place from density value
    size_t numSynapses = (size_t)(length * density);

//  Save subset coarse grid edges in a vector
    std::vector<Edge*> vEdges(m_spSubsetHandler->begin<Edge>(si, 0),
                              m_spSubsetHandler->end<Edge>(si, 0));

//  Call wrapper
    place_synapses_uniform(vEdges, numSynapses, type);
}

void SynapseDistributor::place_synapses_uniform(const char* subset, number density, const std::string& type)
{
    int si = m_spSubsetHandler->get_subset_index(subset);
    place_synapses_uniform(si, density, type);
}


void SynapseDistributor::place_synapses
(
    const std::vector<number>& distr,
    size_t numSynapses,
    const std::string& type
)
{
    // check validity of input (sum of probabilities has to be 1, number of probs needs to match number of subsets)
    number s = 0;
    for (size_t i = 0; i < distr.size(); ++i)
        s += distr[i];

    UG_COND_THROW(fabs(s - 1.0) > 1e-5, "Specified distribution incorrect. Probabilities don't sum up to 1.");
    UG_COND_THROW((int) distr.size() != m_spSubsetHandler->num_subsets(),
        "Specified distribution incorrect. Number of subsets is " << m_spSubsetHandler->num_subsets()
        << ", but number of given probabilities is " << distr.size() << ".");

    // placing synapses
    for (size_t i = 0; i < distr.size(); ++i)
    {
        size_t numSynPerSubset = numSynapses*distr[i];
        std::vector<Edge*> vEdges(m_spSubsetHandler->begin<Edge>(i, 0), m_spSubsetHandler->end<Edge>(i, 0));

        place_synapses_uniform(vEdges, numSynPerSubset, type);
    }

    // handle multi-grid transfer
    if (m_bDomBased)
        CopySynapsesToAllLevels();
}


void SynapseDistributor::place_synapses_uniform
(
    number density,
    number x,
    number y,
    number z,
    number radius,
    const std::string& type
)
{
    // convert center coords to vector
    vector3 center(x, y, z);

    // find all edges within ball and add their lengths
    number length = 0.0;
    std::vector<Edge*> vEdges;
    EdgeIterator it = m_spGrid->begin<Edge>(0);
    EdgeIterator it_end = m_spGrid->end<Edge>(0);
    for (; it != it_end; ++it)
    {
        Edge* e = *it;

        vector3 a = m_aaPosition[e->vertex(0)];
        vector3 b = m_aaPosition[e->vertex(1)];
        if (VecDistanceSq(a, center) < radius*radius && VecDistanceSq(b, center) < radius*radius)
        {
            length += EdgeLength(e, m_aaPosition);
            vEdges.push_back(e);
        }
    }

    // place synapses on in-ball edges
    size_t numSynapses = (size_t) (length * density);
    place_synapses_uniform(vEdges, numSynapses, type);
}


void SynapseDistributor::degenerate_uniform(const std::vector<Edge*>& vEdges, size_t numSynapses)
{
//	Determine number of synapses in the given vector of edges
	if(numSynapses > num_synapses(vEdges))
		UG_THROW("SynapseDistributor::degenerate_uniform(vector<Edge*> vEdges, size_t numSynapses): Cannot degenerate more synapses than placed.");

//	Determine discrete random distribution of given edges by their lengths
	std::vector<number> probs;
	number totLength = 0.0;

	for(size_t i = 0; i < vEdges.size(); ++i)
	{
		Edge* e = vEdges[i];
		number edgeLength = EdgeLength(e, m_aaPosition);

	//	Push back current edge length as probability (still to be divided by total length of all edges)
		probs.push_back(edgeLength);

		totLength += edgeLength;
	}

//	Scale probabilities
	VecScale(probs, probs, 1./totLength);

//	Define random int generator for random access of edges vector
	boost::random::mt19937 gen;
	boost::random::discrete_distribution<> dist(probs.begin(), probs.end());
	boost::variate_generator<boost::random::mt19937, boost::random::discrete_distribution<> > randomIndex(gen, dist);

//	Degenerate specified number of synapses along the coarse grid edges randomly (uniformly s.t. individual edge lengths)
	size_t i = 0;
	while (i < numSynapses)
	{
		size_t rID = randomIndex();

		// TODO: There surely is a way to do this without so many misses!
		if (m_aaSSyn[vEdges[rID]].size() > 1)
		{
		    remove_synapse(vEdges[rID]);
            ++i;
		}
	}

//	Handle multigrid transfer
	if(m_bDomBased)
		CopySynapsesToAllLevels();
}



void SynapseDistributor::degenerate_uniform(number p)
{
	if(p < 0 || p > 1)
		UG_THROW("SynapseDistributor::degenerate_uniform(number p, size_t numSynapses): Specified percentage incorrect.");

//	Save coarse grid edges in a vector
	std::vector<Edge*> vEdges = std::vector<Edge*>(m_spGrid->begin<Edge>(0), m_spGrid->end<Edge>(0));

//	Determine number of synapses to be removed from subset
	size_t numSynToRemove = p * num_synapses(vEdges);

	degenerate_uniform(vEdges, numSynToRemove);
}

void SynapseDistributor::degenerate_uniform(number p, int si)
{
	if(p < 0 || p > 1)
		UG_THROW("SynapseDistributor::degenerate_uniform(number p, int subsetIndex): Specified percentage incorrect.");

//	Save subset coarse grid edges in a vector
	std::vector<Edge*> vEdges = std::vector<Edge*>(m_spSubsetHandler->begin<Edge>(si, 0), m_spSubsetHandler->end<Edge>(si, 0));

//	Determine number of synapses to be removed from subset
	size_t numSynToRemove = p * num_synapses(vEdges);

	degenerate_uniform(vEdges, numSynToRemove);
}

void SynapseDistributor::degenerate_uniform(number p, const char* subset)
{
	if(p < 0 || p > 1)
		UG_THROW("SynapseDistributor::degenerate_uniform(number p, const char* subset): Specified percentage incorrect.");

//	Get subset index from subset name
	int si = m_spSubsetHandler->get_subset_index(subset);

	degenerate_uniform(p, si);
}

size_t SynapseDistributor::num_synapses(const std::vector<Edge*>& vEdges) const
{
    // determine number of synapses in the given vector of edges
	size_t numSynapses = 0;
	for (size_t i = 0; i < vEdges.size(); ++i)
	    numSynapses += m_aaSSyn[vEdges[i]].size();

	return numSynapses / 2; // we have counted pre- AND post-synapses
}

size_t SynapseDistributor::num_synapses(int si) const
{
//	Save subset coarse grid edges in a vector
	std::vector<Edge*> vEdges;
	if (m_bDomBased)
	    {
	        ConstEdgeIterator it = m_spGrid->begin<Edge>();
	        ConstEdgeIterator it_end = m_spGrid->end<Edge>();
	        for (; it != it_end; ++it)
	        {
	            Edge* e = *it;

	            // only surface edges; as we only use distributor in serial, this is easily verified
	            // only from correct subset
	            if (!m_spGrid->has_children(e) && m_spSubsetHandler->get_subset_index(e) == si)
	                vEdges.push_back(e);
	        }
	    }
	    else
	    {
	        // if we have created the grid ourselves, it has only 1 level
	        std::vector<Edge*> vE(m_spSubsetHandler->begin<Edge>(si, 0), m_spSubsetHandler->end<Edge>(si, 0));
	        vEdges.swap(vE);
	    }

	return num_synapses(vEdges);
}

size_t SynapseDistributor::num_synapses(const char* subset) const
{
//	Get subset index from subset name
	int si = m_spSubsetHandler->get_subset_index(subset);

	return num_synapses(si);
}

size_t SynapseDistributor::num_synapses() const
{
    // save surface grid edges in a vector
    std::vector<Edge*> vEdges;
    if (m_bDomBased)
    {
        ConstEdgeIterator it = m_spGrid->begin<Edge>();
        ConstEdgeIterator it_end = m_spGrid->end<Edge>();
        for (; it != it_end; ++it)
        {
            // only surface edges
            // as we only use distributor in serial, this is easily verified
            if (!m_spGrid->has_children(*it))
                vEdges.push_back(*it);
        }
    }
    else
    {
        // if we have created the grid ourselves, it has only 1 level
        std::vector<Edge*> vE(m_spGrid->begin<Edge>(0), m_spGrid->end<Edge>(0));
        vEdges.swap(vE);
    }

	return num_synapses(vEdges);
}


size_t SynapseDistributor::max_id() const
{
    int maxID = -1;
    if (m_bDomBased)
    {
        ConstEdgeIterator it = m_spGrid->begin<Edge>();
        ConstEdgeIterator it_end = m_spGrid->end<Edge>();
        for (; it != it_end; ++it)
        {
           // only surface edges; as we only use distributor in serial, this is easily verified
           if (m_spGrid->has_children(*it)) continue;

           const std::vector<IBaseSynapse*>& syns = m_aaSSyn[*it];
           size_t sz = syns.size();
           for (size_t i = 0; i < sz; ++i)
               maxID = std::max((int) syns[i]->id(), maxID);
        }
    }
    else
    {
        ConstEdgeIterator it = m_spGrid->begin<Edge>(0);
        ConstEdgeIterator it_end = m_spGrid->end<Edge>(0);
        for (; it != it_end; ++it)
        {
            const std::vector<IBaseSynapse*>& syns = m_aaSSyn[*it];
            size_t sz = syns.size();
            for (size_t i = 0; i < sz; ++i)
               maxID = std::max((int) syns[i]->id(), maxID);
        }
    }

    return (size_t) maxID;
}


void SynapseDistributor::print_status()
{
	// Iterates over all subsets and prints the number of edges in each subset
	UG_LOGN("#################################### Status ####################################");
	for (int i=0; i<m_spSubsetHandler->num_subsets(); ++i) {
		UG_LOGN("Subset '" << m_spSubsetHandler->subset_info(i).name << "(" << i << ")':\t\tEdges: "
			<< m_spSubsetHandler->num<Edge>(i) << "\tSynapses: " << num_synapses(i));
	}
	UG_LOGN(std::endl << "Total:\t\t\tEdges: " << m_spGrid->num_edges()
	    << "\tSynapses: " << num_synapses() << std::endl
		<< "################################################################################");
}


bool SynapseDistributor::export_grid(const std::string& outfile)
{
	m_LastMessage = "export_grid(): Saving to "+outfile+"...";

	if(SaveGridToFile(*m_spGrid, *m_spSubsetHandler, outfile.c_str())) {
		m_LastMessage.append("done.");
		return true;
	}
	else {
		m_LastMessage.append("SaveGridToFile(*pm_Grid, *m_spSubsetHandler, outfile.c_str()) failed.");
		return false;
	}
}


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
