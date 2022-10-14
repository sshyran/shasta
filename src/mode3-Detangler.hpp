#ifndef SHASTA_MODE3_DETANGLER_HPP
#define SHASTA_MODE3_DETANGLER_HPP

// Shasta.
#include "ReadId.hpp"

// Standard library.
#include <map>
#include "vector.hpp"

/*******************************************************************************

Class mode3::Detangler contains data structures and code used to detangle the
mode3::AssemblyGraph.

*******************************************************************************/

namespace shasta {
    namespace mode3 {
        class Detangler;

        class AssemblyGraph;
    }
}



class shasta::mode3::Detangler {
public:

    // Type used to describe the sequence of segment ids
    // visited by an oriented read.
    using Journey = vector<uint64_t>;

    // The journey of each oriented read.
    // Obtained from the AssemblyGraph::assemblyGraphJourneys.
    // Indexed by OrientedReadId::getValue().
    // This does not change during detangle.
     vector<Journey> journeys;

    // Type used to identify an entry in a journey.
    class JourneyEntryId {
    public:
        OrientedReadId orientedReadId;

        // The position of this entry in the journey of this oriented read.
        uint64_t position;
    };

    // A cluster is a set of journey entries all corresponding to the same
    // segment id.
    // At the beginning, there is a cluster for each segment.
    // During detangling, clusters are split, so there can be more
    // than one cluster for each segment.
    class Cluster {
    public:
        uint64_t segmentId;
        vector<JourneyEntryId> journeyEntryIds;
    };

    // Store the clusters keyed by segmentId.
    std::map<uint64_t, vector<Cluster> > clusters;

    // The clusterTable gives the cluster that each journey entry corresponds to.
    // It changes during detangling.
    // It is indexed just like the journeys,
    // that is [orientedReadId.getValue()][position].
    vector< vector<const Cluster*> > clusterTable;

    Detangler(const AssemblyGraph&);
private:
    void createJourneys(const AssemblyGraph&);
    void createInitialClusters();
};



#endif

