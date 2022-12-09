#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP

// See comments in mode3a.hpp.

// Shasta.
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <list>
#include "tuple.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3a {

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            AssemblyGraphVertex, AssemblyGraphEdge>;

        class JourneyEntry;
        class Transition;

        class PackedMarkerGraph;
    }
}



// An entry in the journey of an oriented read.
class shasta::mode3a::JourneyEntry {
public:
    OrientedReadId orientedReadId;

    // The position in the journey for this oriented read.
    uint64_t position;
};



// Each vertex represents a segment replica.
class shasta::mode3a::AssemblyGraphVertex {
public:

    // The PackedMarkerGraph segment that this vertex corresponds to.
    uint64_t segmentId;

    // Replica index among all vertices with the same segmentId.
    uint64_t segmentReplicaIndex;

    AssemblyGraphVertex(uint64_t segmentId, uint64_t segmentReplicaIndex=0) :
        segmentId(segmentId), segmentReplicaIndex(segmentReplicaIndex) {}

    string stringId() const
    {
        string s = to_string(segmentId);
        if(segmentReplicaIndex != 0) {
            s += "." + to_string(segmentReplicaIndex);
        }
        return s;
    }

    // The journey entries that go through this vertex.
    vector<JourneyEntry> journeyEntries;
};



// A transition of an oriented read from a segment to the next.
class shasta::mode3a::Transition {
public:
    uint64_t position0;
    uint64_t position1;
    OrientedReadId orientedReadId;

    bool operator<(const Transition& that) const
    {
        return
            tie(orientedReadId, position0, position1) <
            tie(that.orientedReadId, that.position0, that.position1);
    }
    bool operator==(const Transition& that) const
    {
        return
            tie(orientedReadId, position0, position1) ==
            tie(that.orientedReadId, that.position0, that.position1);
    }
};



// Each edge represents a link.
class shasta::mode3a::AssemblyGraphEdge {
public:
};



class shasta::mode3a::AssemblyGraph : public AssemblyGraphBaseClass {
public:
    AssemblyGraph(const PackedMarkerGraph&);
    void write(const string& name) const;

    // The sequence of vertices visited by each oriented read.
    // It has one entry for each entry of the corresponding journey in the PackedMarkerGraph.
    // An entry can be null_vertex(), if the corresponding vertex was removed
    // during detangling.
    // Indexed by OrientedReadId::getValue().
    vector< vector<vertex_descriptor> > journeys;

    // Each segment in the AssemblyGraph corresponds to a segment in
    // this PackedMarkerGraph.
    const PackedMarkerGraph& packedMarkerGraph;
private:

    // The vertices created so far for each segmentId.
    // Some of these may have been deleted and set to null_vertex().
    vector< vector<vertex_descriptor> > verticesBySegment;

    void createSegmentsAndJourneys();
    void createLinks();

    // Get the stringId for a given vertex_descriptor, or "None" if v is null_vertex().
    string vertexStringId(vertex_descriptor v) const;

    // Get the transitions for an edge.
    void getEdgeTransitions(edge_descriptor, vector<Transition>&) const;
    uint64_t edgeCoverage(edge_descriptor) const;

    void writeGfa(const string& name, uint64_t minLinkCoverage) const;
    void writeLinkCoverageHistogram(const string& name) const;
    void writeJourneys(const string& name) const;

    // Simple detangling, one vertex at a time, looking only
    // at immediate parent and children.
public:
    void simpleDetangle(
        uint64_t minLinkCoverage,
        uint64_t minTangleCoverage);
private:
    void simpleDetangle(
        vertex_descriptor,
        uint64_t minLinkCoverage,
        uint64_t minTangleCoverage);

    // Find the previous and next vertex for each JourneyEntry in a given vertex.
    // On return, adjacentVertices contains a pair of vertex descriptors for
    // each JourneyEntry in vertex v, in the same order.
    // Those vertex descriptors are the previous and next vertex visited
    // by the oriented read for that JourneyEntry, and can be null_vertex()
    // if v is at the beginning or end of the journey of an oriented read.
    void findAdjacentVertices(
        vertex_descriptor v,
        vector< pair<vertex_descriptor, vertex_descriptor> >& adjacentVertices
    ) const;


    // Compute paths by following oriented reads in a starting vertex.
    // For now the paths are not stored anywhere.
    void computePaths(uint64_t threadCount);
    void computePathsThreadFunction(uint64_t threadId);
    void computePath(
        vertex_descriptor,
        uint64_t minSegmentCoverage,
        uint64_t minLinkCoverage
        );
};

#endif

